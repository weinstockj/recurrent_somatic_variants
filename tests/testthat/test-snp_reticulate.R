context("interface to sparse python package")

test_that("interface to reticulate sparse works", {
    have_sparse = reticulate::py_module_available("sparse")
    have_numpy = reticulate::py_module_available("numpy")

    if(!(have_sparse & have_numpy)) {
        skip("sparse and numpy python modules not available for testing")
    }
    
    reticulate::py_run_string("import sparse as snp; import numpy as np; import gc")
    reticulate::py_run_string("np_mat = np.random.rand(30, 20)")
    reticulate::py_run_string("sp_mat = snp.COO(np_mat)")
    reticulate::py_run_string("snp.save_npz('test.npz', sp_mat)")
    reticulate::py_run_string("vaf_table = snp.load_npz('test.npz')")

    # R indices [1, 2, 3, 4, 5] gets converted to [0, 1, 2, 3, 4]
    # R indices [2, 3] gets converted to [1, 2]
    submatrix = subset_vaf_table(1:5, 2:3)
    # print(submatrix)

    R_col_sums = colSums(submatrix)
    R_row_sums = rowSums(submatrix)
    # print(" R col sums")
    # print(dput(R_col_sums))

    python_col_sums = reticulate::py_eval("vaf_table[[0, 1, 2, 3, 4], :][:, [1, 2]].sum(axis = 0).todense()")
    python_row_sums = reticulate::py_eval("vaf_table[[0, 1, 2, 3, 4], :][:, [1, 2]].sum(axis = 1).todense()")

    # subset_vaf_table always returns a matrix even for scalar
    R_cell_value = as.numeric(subset_vaf_table(18, 8))
    python_cell_value = reticulate::py_eval("vaf_table[17, 7]")

    # print("python col sums")
    # print(dput(python_col_sums))
    # use as.numeric to strip off attributes
    expect_equal(R_col_sums, as.numeric(python_col_sums))

    expect_equal(R_cell_value, python_cell_value)
    # clean up 
    file.remove("test.npz")
    reticulate::py_run_string("del np_mat; del sp_mat; del vaf_table; gc.collect()")
})
