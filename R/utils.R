assert_not_empty = function(df) {
    stopifnot(is.data.frame(df))
    stopifnot(nrow(df) > 0)
}

inverse_normalize = function(x) {
    qnorm(rank(x, na.last = "keep") / (sum(!is.na(x)) + 1))
}

