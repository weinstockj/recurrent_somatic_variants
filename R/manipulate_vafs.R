identify_most_common_mutations = function(sample_meta_path, variant_meta_path, allele_threshold = 3000) {
    vaf_mat = convert_to_dgCMatrix() 
    vaf_mat = (vaf_mat > 0) * 1L # convert to matrix of 0 and 1 when value > 0 (i.e. sample i with variant j has the mutation
    variant_meta = vroom::vroom(variant_meta_path) %>%
        tibble::rowid_to_column(var = "col_index") %>%
        dplyr::select(CHROM, POS, REF, ALT, col_index) %>%
        dplyr::mutate(variant_label = glue::glue("{CHROM}-{POS}-{REF}-{ALT}"))

    assert_not_empty(variant_meta)

    sample_meta = vroom::vroom(sample_meta_path) %>%
        tibble::rowid_to_column(var = "row_index") %>%
        dplyr::select(Sample, row_index)

    assert_not_empty(sample_meta)

    col_sums = Matrix::colSums(vaf_mat)
    indices = which(col_sums >= allele_threshold)

    vaf_mat = vaf_mat[, indices]
    variant_meta = variant_meta %>% dplyr::filter(col_index %in% indices)

    vaf_ped = tibble::as_tibble(as.matrix(vaf_mat)) %>%
                setNames(variant_meta$variant_label)

    sample_meta %>%
        dplyr::bind_cols(vaf_ped) %>%
        dplyr::select(-row_index)
}

#' @title Extract the VAFs of a specific set of variants
#' @param sample_meta_path a file path
#' @param variant_meta_path a file path
#' @param variant_id a character vector containing one or more requested variants in the format CHROM-POS-REF-ALT
#' @param return_genotypes a logical value. If TRUE, returns 0 or 1 values. If FALSE, returns VAFs.
#' @return a tibble where the first column is the sample identifier, and all remaining columns are the variant mutations or VAFs
#' @export
subset_vaf_matrix_by_id = function(sample_meta_path, variant_meta_path, variant_id, return_genotypes = FALSE, normalize_columns = FALSE, return_sparse_matrix = FALSE) {

    stopifnot(!missing(variant_id) & !missing(sample_meta_path) & !missing(variant_meta_path))
    stopifnot(is.character(variant_id))

    vaf_mat = convert_to_dgCMatrix() 
    if(return_genotypes) {
        vaf_mat = (vaf_mat > 0) * 1L # convert to matrix of 0 and 1 when value > 0 (i.e. sample i with variant j has the mutation
    }

    variant_meta = vroom::vroom(variant_meta_path) %>%
        tibble::rowid_to_column(var = "col_index") %>%
        dplyr::select(CHROM, POS, REF, ALT, col_index) %>%
        dplyr::mutate(variant_label = glue::glue("{CHROM}-{POS}-{REF}-{ALT}"))

    assert_not_empty(variant_meta)

    sample_meta = vroom::vroom(sample_meta_path) %>%
        tibble::rowid_to_column(var = "row_index") %>%
        dplyr::select(Sample, row_index)

    assert_not_empty(sample_meta)

    flog.info(glue::glue("identified {length(variant_id)} requested variants"))

    variant_meta = dplyr::filter(variant_meta, variant_label %in% variant_id)

    flog.info(glue::glue("after subsetting to requested variants, we have {nrow(variant_meta)} variants"))

    assert_not_empty(variant_meta)

    vaf_mat = vaf_mat[, variant_meta$col_index]


    if(normalize_columns) {
        col_sums = Matrix::colSums(vaf_mat)
        vaf_mat = vaf_mat %*% Matrix::diag(1 / col_sums)
        stopifnot(all(abs(Matrix::colSums(vaf_mat) - 1) < 1e-4))
    }

    if(return_sparse_matrix) {
        return(vaf_mat)
    } else {

        vaf_ped = tibble::as_tibble(as.matrix(vaf_mat)) %>%
                    setNames(variant_meta$variant_label)

        ped = sample_meta %>%
            dplyr::bind_cols(vaf_ped) %>%
            dplyr::select(-row_index)

        return(ped)
    }

}

#' @title Compute mutation burden per sample
#' @param sample_meta_path a file path
#' @param variant_meta_path a file path
#' @param variant_id a character vector containing one or more requested variants in the format CHROM-POS-REF-ALT
#' @param return_genotypes a logical value. If TRUE, returns 0 or 1 values. If FALSE, returns VAFs.
#' @return a tibble where the first column is the sample identifier, and all remaining columns are the variant mutations or VAFs, except for the AC column, which is the sum of the individual variants
#' @export
compute_mutation_counts_per_sample = function(sample_meta_path, variant_meta_path, variant_id, include_variants = FALSE, normalize_columns = FALSE) {
    sample_genotypes = subset_vaf_matrix_by_id(
                sample_meta_path, 
                variant_meta_path, 
                variant_id, 
                return_genotypes = TRUE,
                normalize_columns = normalize_columns
    ) 

    assert_not_empty(sample_genotypes)

    sample_column = "Sample"

    sample_genotypes = sample_genotypes %>%
        dplyr::mutate(
            AC = rowSums(dplyr::across(-{{sample_column}}))
        )

    if(include_variants) {
        return(sample_genotypes) 
    } else {
        return(sample_genotypes %>% dplyr::select({{sample_column}}, AC))
    }
}

subset_vaf_matrix_by_sample_id = function(sample_meta_path, variant_meta_path, sample_id, return_genotypes = FALSE) {

    stopifnot(!missing(sample_id) & !missing(sample_meta_path) & !missing(variant_meta_path))
    stopifnot(is.character(sample_id))

    vaf_mat = convert_to_dgCMatrix() 
    if(return_genotypes) {
        vaf_mat = (vaf_mat > 0) * 1L # convert to matrix of 0 and 1 when value > 0 (i.e. sample i with variant j has the mutation
    }

    assert_not_empty(variant_meta)

    sample_id_col  = "Sample"
    sample_meta = vroom::vroom(sample_meta_path) %>%
        tibble::rowid_to_column(var = "row_index") %>%
        dplyr::select({{sample_id_col}}, row_index)

    assert_not_empty(sample_meta)

    flog.info(glue::glue("identified {length(sample_id)} requested samples"))

    sample_meta = dplyr::filter(sample_meta, .data[[sample_id_col]] %in% sample_id)

    flog.info(glue::glue("after subsetting to requested samples, we have {nrow(sample_meta)} samples"))

    assert_not_empty(sample_meta)

    vaf_mat = vaf_mat[sample_meta$row_index, ]

    return(vaf_mat)
}
