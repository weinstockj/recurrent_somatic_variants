# global reference (will be initialized in .onLoad)
snp <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
    snp <<- reticulate::import("sparse", delay_load = TRUE)
}

vaf_location = function() {
    "/net/topmed2/working/jweinstk/count_singletons/new_drivers/trimmed_vaf_47k_2020_09_06.npz"
}

sample_meta_location = function() {
    "/net/topmed2/working/jweinstk/count_singletons/new_drivers/output/output_sample_statistics_hs_2020_09_05.tsv"
}

variant_meta_location = function() {
    "/net/topmed2/working/jweinstk/count_singletons/new_drivers/output/output_variant_summary_statistics_hs_2020_09_05.tsv"
}

fai_location = function() {
    "/net/topmed2/working/jweinstk/reference_comparison/hs38DH.fa.fai"
}

subset_vaf_table = function(rows, cols) {
    stopifnot(all(rows > 0) & all(cols > 0))

    #convert to python 0-based indices
    rows = rows - 1
    cols = cols - 1 

    py_rows_iterable = paste(rows, collapse = ", ")
    py_cols_iterable = paste(cols, collapse = ", ")
    py_cmd = glue::glue("vaf_table[[{py_rows_iterable}], :][:, [{py_cols_iterable}]].todense()")
    vaf = reticulate::py_eval(py_cmd)
    return(vaf)
}

convert_to_dgCMatrix = function() {
    # converts to Matrix dgCMatrix format
    vaf = reticulate::py_to_r(reticulate::py_eval("vaf_table.tocsc()"))
    return(vaf)
}

read_in_vaf = function(vaf_path) {
    # vaf_table = snp$load_npz(vaf_path)
    py_cmd = glue::glue("import sparse as snp; vaf_table = snp.load_npz('{vaf_path}')")
    reticulate::py_run_string(py_cmd, local = FALSE, convert = FALSE)
    n_samples = reticulate::py_eval("vaf_table.shape[0]")
    n_variants = reticulate::py_eval("vaf_table.shape[1]")
    futile.logger::flog.info(glue::glue("vaf matrix has {n_samples} samples and {n_variants} variants"))
}

recode_vaf_by_threshold = function(threshold = .26) {
    # remove variants with VAF values above threshold
    py_cmd = glue::glue("vaf_table.tocsc(); vaf_table[vaf_table > {threshold}] = 0.0; vaf_table.eliminate_zeros();vaf_mat = snp.as_coo(vaf_mat)")
}

#' @export
#' @example \dontrun{
#' high_level_helper(sample_ped_path = "../prepare_ped_files/blood_phenotypes_2020_11_25.ped", outcome_column= "hemoglobin_mcnc_bld", covariate_columns = "AgeAtBloodDraw,INFERRED_SEX,STUDY,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,sPC1,sPC2,sPC3,sPC4,sPC5,sPC6,sPC7,sPC8,sPC9,sPC10,VB_DEPTH,FREEMIX", output_prefix = "hemoglobin_mcnc_bld", n_workers = 200)
#' }
high_level_helper = function(
    vaf_path = vaf_location(), 
    variant_meta_path = variant_meta_location(), 
    sample_meta_path = sample_meta_location(), 
    sample_ped_path = test_ped(),
    ref_fai_path = fai_location(),
    region_size = 1e7,
    n_workers = 20,
    run_locally = FALSE,
    slurm_template = "batchtools.slurm.tmpl",
    sample_id_col = "Sample",
    analysis_type = "quantitative",
    inverse_normalize = TRUE,
    covariate_columns = "",
    outcome_column = "",
    output_prefix = "output") {

    futile.logger::flog.info("R library site:")
    futile.logger::flog.info(.libPaths())

    acceptable_chroms = c(glue::glue("chr{1:22}"), "chrX")

    fai = vroom::vroom(ref_fai_path, col_names = c("contig", "size", "location", "bases", "bytes")) %>%
        dplyr::select(contig, size) %>%
        dplyr::filter(contig %in% acceptable_chroms)

    chunks = purrr::map2_dfr(fai$contig, fai$size, ~{
        sequence = seq(1, .y - region_size, by = region_size)
        tibble::tibble(chrom = .x, start = sequence) %>%
            dplyr::mutate(
                end = dplyr::lead(start) - 1,
                end = dplyr::coalesce(end, .y) # last end coordinate is size of chrom
            )
    })

    futile.logger::flog.info(glue::glue("now running {nrow(chunks)} separate chunks using {n_workers} workers (cores or size of slurm job array)"))

    if(run_locally) {
        # future::plan(future::sequential)
        future::plan(future::multiprocess, workers = n_workers)
    } else {
        future::plan(
            future.batchtools::batchtools_slurm, 
            template = slurm_template,
            resources = list(
                partition = "topmed", array_throttle = n_workers, ncpus = 2, memory = 5000, walltime = 1200 * 60, omp.threads = 1, blas.threads = 1)
        )
    }

    print("debugging info")
    mi <- list(rver = getRversion(), libs = .libPaths())
    print(mi)

    wi %<-% list(rver = getRversion(), libs = .libPaths())
    print(wi)

    results = future.apply::future_lapply(
    # results = lapply(
        1:nrow(chunks),
        function(i) {
            chunk_row = chunks %>% dplyr::slice(i)

            result = main_assoc(
                vaf_path = vaf_path,
                variant_meta_path = variant_meta_path,
                sample_meta_path = sample_meta_path, 
                sample_ped_path = sample_ped_path,
                sample_id_col = sample_id_col,
                analysis_type = analysis_type,
                inverse_normalize = inverse_normalize,
                covariate_columns = covariate_columns,
                outcome_column = outcome_column,
                output_prefix = output_prefix,
                chrom = chunk_row[["chrom"]],
                start = chunk_row[["start"]],
                end = chunk_row[["end"]]
            )
            result
        },
        future.packages = "somaticWas"
    )

    future::plan(future::sequential)

    results_bind = results %>%
                    purrr::compact(.) %>%
                    dplyr::bind_rows(.) %>%
                    dplyr::mutate(POS = as.integer(POS))

    results_bind_variants = results_bind %>%
                    dplyr::filter(term == tested_variant)

    arrow::write_feather(results_bind, file.path(output_prefix, "association_statistics_full.feather"))
    arrow::write_feather(
        results_bind_variants, 
        file.path(output_prefix, "association_statistics.feather")
    )
    tsv_fname = file.path(output_prefix, "association_statistics.tsv")
    readr::write_tsv(results_bind_variants, tsv_fname)
    system(glue::glue("bgzip -c {tsv_fname} > {tsv_fname}.gz"))
    system(glue::glue("tabix -S1 -s1 -b2 -e2 {tsv_fname}.gz"))
    cleanup(output_prefix)

    tsv_fname = file.path(output_prefix, "association_statistics_top_100.tsv")
    readr::write_tsv(results_bind_variants %>% dplyr::arrange(p.value) %>% dplyr::slice(1:100), tsv_fname)

    manhattan(results_bind_variants, output_prefix)
    qqunif_plot_save(results_bind_variants %>% tidyr::drop_na(p.value), output_prefix)

    futile.logger::flog.info("now done.")

    return(results_bind)
}

cleanup = function(output_prefix) {
    intermediate_feather = list.files(output_prefix, pattern = glob2rx("association_statistics_chr*.feather"), full.names = TRUE)
    purrr::walk(intermediate_feather, file.remove)
    
    system2("cat", 
        args = file.path(output_prefix, "association_statistics_chr*.log"), 
        stdout = file.path(output_prefix, "intermediate_steps.log")
    )

    intermediate_logs = list.files(output_prefix, pattern = glob2rx("association_statistics_chr*.log"), full.names = TRUE)
    purrr::walk(intermediate_logs, file.remove)
}

main_assoc = function(
    vaf_path = vaf_location(), 
    variant_meta_path = variant_meta_location(), 
    sample_meta_path = sample_meta_location(), 
    sample_ped_path = test_ped(),
    sample_id_col = "Sample",
    analysis_type = "quantitative",
    inverse_normalize = TRUE,
    covariate_columns = "",
    outcome_column = "",
    output_prefix = "output",
    chrom = "chr1",
    start = 1,
    end = 1e5) {

    stopifnot(file.exists(vaf_path))
    stopifnot(file.exists(variant_meta_path))
    stopifnot(file.exists(sample_meta_path))

    analysis_type = match.arg(analysis_type, c("quantitative", "binary"))
    IS_QUANT = analysis_type == "quantitative"

    Sys.setenv("OPENBLAS_NUM_THREADS" = 1)

    if(is.null(output_prefix) | output_prefix == "" | missing(output_prefix)) {
        output_prefix = "output/"
    }
    if(!dir.exists(output_prefix)) {
        dir.create(output_prefix)
    }

    logger_fname = file.path(output_prefix, glue::glue("association_statistics_{chrom}_{as.integer(start)}_{as.integer(end)}.log"))
    if(file.exists(logger_fname)) file.remove(logger_fname)
    appender = futile.logger::appender.file(logger_fname)
    futile.logger::flog.appender(appender, name = "logger")

    read_in_vaf(vaf_path)
    recode_vaf_by_threshold()

    variant_meta = vroom::vroom(variant_meta_path) %>%
        tibble::rowid_to_column(var = "col_index") %>%
        dplyr::select(CHROM, POS, REF, ALT, col_index)

    sample_meta = vroom::vroom(sample_meta_path) %>%
        tibble::rowid_to_column(var = "row_index") %>%
        dplyr::select(Sample, row_index)

    assert_not_empty(variant_meta)
    assert_not_empty(sample_meta)

    futile.logger::flog.info(glue::glue("variant meta file has {nrow(variant_meta)} rows"), name = "logger")

    futile.logger::flog.info(glue::glue("subsetting to {chrom}:{start}-{end}"), name = "logger")

    region = variant_meta %>%
        dplyr::filter(CHROM == {{chrom}} & POS >= {{start}} & POS <= {{end}}) %>%
        dplyr::mutate(variant_label = glue::glue("{CHROM}-{POS}-{REF}-{ALT}"))

    if(nrow(region) == 0) {
        warning("VAF matrix has no variants in this region")
        return(NULL)
    }

    futile.logger::flog.info(glue::glue("variant meta file has {nrow(region)} rows after subsetting to region"), name = "logger")

    futile.logger::flog.info(glue::glue("sample meta file has {nrow(sample_meta)} rows"), name = "logger")

    ped = vroom::vroom(sample_ped_path)
    futile.logger::flog.info(glue::glue("ped file has {nrow(ped)} rows"), appender = appender)

    ped = dplyr::inner_join(sample_meta, ped, by = sample_id_col)
    futile.logger::flog.info(glue::glue("ped file has {nrow(ped)} rows after merging with sample meta"), name = "logger")

    covariate_columns = stringr::str_split(covariate_columns, ",")[[1]]
    stopifnot(all(covariate_columns %in% names(ped)))

    ped_drop = tidyr::drop_na(ped, dplyr::all_of(covariate_columns), dplyr::all_of(outcome_column))
    futile.logger::flog.info(glue::glue("ped file has {nrow(ped_drop)} rows after dropping missing data"), name = "logger")

    if(inverse_normalize & IS_QUANT) {
        futile.logger::flog.info(glue::glue("now applying inverse normal transform to {outcome_column}"), name = "logger")
        ped_drop[[outcome_column]] = inverse_normalize(ped_drop[[outcome_column]])
    }

    count_unique_values = ped_drop %>%
        dplyr::summarize(dplyr::across(where(is.character), ~length(unique(.x)))) %>%
        tidyr::gather(key = column, value = n_unique)

    constant_columns = count_unique_values %>%
                        dplyr::filter(n_unique == 1) %>%
                        dplyr::pull(column)

    if(length(constant_columns) >= 1) {
        futile.logger::flog.info(
            glue::glue("dropping the following columns because they are constants {paste(constant_columns, collapse = ', ')}")
        )
    }

    covariate_columns = setdiff(covariate_columns, constant_columns)
    covariate_string = paste(covariate_columns, collapse = " + ")

    local_vaf = subset_vaf_table(ped_drop$row_index, region$col_index) %>%
                    tibble::as_tibble(.) %>% 
                    setNames(region$variant_label)

    futile.logger::flog.info(glue::glue("region has {ncol(local_vaf)} variants"), name = "logger")

    if(IS_QUANT) {
        results = purrr::map_dfr(region$variant_label, ~{
            variant = local_vaf[[.x]]     
            VAF_SUM = sum(variant)
            AC = sum(variant > 0)
            N_SAMPLES = length(variant)
            ped_drop[[.x]] = variant
            formula = as.formula(glue::glue("{outcome_column} ~ {covariate_string} + `{.x}`"))
            model = lm(formula, data = ped_drop)
            model_summary = broom::tidy(model, conf.int = TRUE)
            model_anova = car::Anova(model, type = "II") %>%
                            broom::tidy(.) %>%
                            dplyr::mutate(
                                term = stringr::str_replace_all(term, "`", ""),
                                partial.r.squared = sumsq / sum(sumsq)
                            ) %>%
                            dplyr::select(term, partial.r.squared)

            model_glance = broom::glance(model)
            model_summary %>%
                dplyr::mutate(
                    term     = stringr::str_replace_all(term, "`", ""),
                    r.squared = model_glance$r.squared,
                    df       = model_glance$df,
                    AIC      = model_glance$AIC,
                    tested_variant = .x,
                    AC       = AC,
                    VAF_SUM  = VAF_SUM,
                    N_SAMPLES = N_SAMPLES
                ) %>%
                dplyr::left_join(model_anova, by = "term") # left instead of inner because categorical covariates need special handling as STUDY-> STUDYaric, STUDYchrs, etc in some output, but remain collapsed in Anova output
        })
    } else {
        results = purrr::map_dfr(region$variant_label, ~{
            variant = local_vaf[[.x]]     
            VAF_SUM = sum(variant)
            AC = sum(variant > 0)
            N_SAMPLES = length(variant)
            ped_drop[[.x]] = variant
            formula = as.formula(glue::glue("{outcome_column} ~ {covariate_string} + `{.x}`"))
            model = glm(formula, data = ped_drop, family = binomial())
            # model = speedglm::speedglm(formula, data = ped_drop, family = binomial())
            model_summary = broom::tidy(model, conf.int = TRUE)

            model_glance = broom::glance(model)
            model_summary %>%
                dplyr::mutate(
                    term     = stringr::str_replace_all(term, "`", ""),
                    AIC      = model_glance$AIC,
                    tested_variant = .x,
                    AC       = AC,
                    VAF_SUM  = VAF_SUM,
                    N_SAMPLES = N_SAMPLES
                )
        })
    }

    results = dplyr::inner_join(results, region, by = c("tested_variant" = "variant_label")) %>%
        dplyr::relocate(dplyr::any_of(c("CHROM", "POS", "REF", "ALT"))) # puts these first

    if(!is.data.frame(results)) {
        return(NULL)
    }

    fname = file.path(output_prefix, glue::glue("association_statistics_{chrom}_{as.integer(start)}_{as.integer(end)}.feather"))
    arrow::write_feather(results, fname)

    return(results)
}

calculate_sample_pcs = function(vaf_path = vaf_location()) {
    start_time = Sys.time()
    futile.logger::flog.info("now coverting to dgCMatrix")
    vaf_mat = convert_to_dgCMatrix()
    n_samples = nrow(vaf_mat)
    n_variants = ncol(vaf_mat)
    n_comps = 50

    n1 = Matrix::Matrix(data = 1, nrow = n_samples, ncol = 1)
    x_row_means = Matrix::rowMeans(vaf_mat)
    x_col_means = Matrix::colMeans(vaf_mat)
    x_col_means_row_vec = Matrix::Matrix(x_col_means, nrow = 1, ncol = n_variants)
    x2_col_means = Matrix::colMeans(vaf_mat * vaf_mat)
    x_col_means2 = x_col_means ^ 2
    x_col_variance = x2_col_means - x_col_means2
    x_col_sd = sqrt(x_col_variance)
    x_scaled = vaf_mat %*% Matrix::diag(1 / x_col_sd)

    futile.logger::flog.info("now computing xxt")

    xxt = Matrix::tcrossprod(x_scaled)
    x_col_means = Matrix::Matrix(Matrix::colMeans(x_scaled), nrow = n_variants, ncol = 1)
    centering_mat = as.numeric(Matrix::crossprod(x_col_means)) * n1 %*% Matrix::t(n1) - x_scaled %*% x_col_means %*% Matrix::t(n1) - n1 %*% Matrix::t(x_col_means) %*% Matrix::t(x_scaled)

    futile.logger::flog.info("now computing svd")
    # centered_xxt = (xxt - nrow(vaf_mat) * Matrix::tcrossprod(x_means)) / (nrow(vaf_mat) - 1)
    pcs = irlba::irlba(
        xxt - centering_mat, 
        nu = n_comps, 
        nv = n_comps, 
        # center = centering_vec, 
        scale = rep(n_variants, n_samples)
        # scale = diag(centered_xxt)
    )

    end_time = Sys.time()
    # 30 minutes?
    futile.logger::flog.info(glue::glue("pc calculation took {end_time - start_time}"))
    return(pcs)
}
