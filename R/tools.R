### common input checks across all tests
# check alternative specification
.check_alternative <- function(alternative, r) {

  if (!alternative %in% c("two.sided", "less", "greater"))
    stop("The alternative must be either 'two.sided', 'less', or 'greater'")

}

# check df
.check_df <- function(df, information_message = "") {

  if(any(df <= 1))
    stop(paste0("Degrees of freedom must be greater than 1. ", information_message))
}
.check_n <- function(n, n_min = 1, information_message = ""){

  if(any(n <= n_min))
    stop(paste0("Sample size must be greater than ", n_min, " ", information_message))
}
