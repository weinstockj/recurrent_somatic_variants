context("Vaf matrix subsetting")

test_that("can subset VAF matrix", {
    have_sparse = reticulate::py_module_available("sparse")
    have_numpy = reticulate::py_module_available("numpy")

    if(!(have_sparse & have_numpy)) {
        skip("sparse and numpy python modules not available for testing")
    }
    
    reticulate::py_run_string("import sparse as snp; import numpy as np; import gc")
    n_sims = 1000
    n_rows = 10
    n_cols = 100
    z_draws = rnorm(n = n_sims, mean = 0.3, sd = 0.01)
    z_mask  = rbinom(n = n_sims, size = 1, prob = .8)
    z_draws = z_draws * z_mask # simulate sparse VAF matrix
    r_vaf_mat = matrix(z_draws, nrow = n_rows, ncol = n_cols)
    write.table(r_vaf_mat, "intermediate__.csv", row.names = FALSE, col.names = FALSE)

    r_sample_meta = tibble::tibble(Sample = 1:n_rows, i = 1:n_rows) # need second column or else vroom throws and Error that delimiter isn't specified
    readr::write_tsv(r_sample_meta, "test_sample_meta__.tsv")
    r_variant_meta = tibble::tibble(
        variant_label = glue::glue("variant_{1:n_cols}"), i = 1:n_cols,
        CHROM = rep("chr1", n_cols),
        POS = 1:n_cols,
        REF = rep("G", n_cols),
        ALT = rep("C", n_cols) # dummy data needed by following function
    )
    readr::write_tsv(r_variant_meta, "test_variant_meta__.tsv")

    reticulate::py_run_string("np_mat = np.loadtxt('intermediate__.csv')")
    reticulate::py_run_string("sp_mat = snp.COO(np_mat)")
    reticulate::py_run_string("snp.save_npz('test.npz', sp_mat)")
    reticulate::py_run_string("vaf_table = snp.load_npz('test.npz')")
    
    test_index = sample(x = 1:n_cols, size = 1)
    test_variant_id = glue::glue("chr1-{test_index}-G-C")
    submatrix = subset_vaf_matrix_by_id("test_sample_meta__.tsv", "test_variant_meta__.tsv", test_variant_id)
    submatrix_geno = subset_vaf_matrix_by_id("test_sample_meta__.tsv", "test_variant_meta__.tsv", test_variant_id, TRUE)
    print(submatrix)

    # use as.numeric to strip off attributes
    expect_equal(r_vaf_mat[, test_index], as.numeric(submatrix[[2]]))
    expect_equal((r_vaf_mat[, test_index] > 0) * 1L, as.numeric(submatrix_geno[[2]]))

    # clean up 
    file.remove("test.npz")
    file.remove("intermediate__.csv")
    file.remove("test_sample_meta__.tsv")
    file.remove("test_variant_meta__.tsv")
    reticulate::py_run_string("del np_mat; del sp_mat; del vaf_table; gc.collect()")
})
