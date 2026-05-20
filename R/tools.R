### common input checks across all tests
# check alternative specification
.check_alternative <- function(alternative, r) {

  if (length(alternative) != 1 || !alternative %in% c("two.sided", "less", "greater"))
    stop("The alternative must be either 'two.sided', 'less', or 'greater'")

}

# check df
.check_df <- function(df, information_message = "") {

  if(!is.numeric(df) || any(!is.finite(df)))
    stop(paste0("Degrees of freedom must contain finite numeric values. ", information_message))

  if(any(df <= 1))
    stop(paste0("Degrees of freedom must be greater than 1. ", information_message))
}
.check_n <- function(n, n_min = 1, information_message = ""){

  if(!is.numeric(n) || any(!is.finite(n)))
    stop(paste0("Sample size must contain finite numeric values. ", information_message))

  if(any(n <= n_min))
    stop(paste0("Sample size must be greater than ", n_min, " ", information_message))
}

.check_r <- function(r){

  if(!is.numeric(r) || length(r) != 1 || !is.finite(r) || r < 1)
    stop("r must be greater than or equal to 1")
}

.check_finite_numeric <- function(x, name){

  if(is.null(x) || length(x) == 0 || !is.numeric(x) || any(!is.finite(x)))
    stop(sprintf("`%s` must contain finite numeric values.", name))

  invisible(x)
}

.check_nonnegative_numeric <- function(x, name){

  .check_finite_numeric(x, name)
  if(any(x < 0))
    stop(sprintf("`%s` must be nonnegative.", name))

  invisible(x)
}

.check_positive_numeric <- function(x, name){

  .check_finite_numeric(x, name)
  if(any(x <= 0))
    stop(sprintf("`%s` must be positive.", name))

  invisible(x)
}

.recycle_stat_input <- function(x, name, target_length){

  if(length(x) == target_length)
    return(x)

  if(length(x) == 1)
    return(rep(x, target_length))

  stop(sprintf("The input length of `%s` must be 1 or match the statistic length.", name))
}
