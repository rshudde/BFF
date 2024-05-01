### common input checks across all tests
# check alternative specification
.check_alternative <- function(alternative, r) {

  if (!alternative %in% c("two.sided", "less", "greater"))
    stop("The alternative must be either 'two.sided', 'less', or 'greater'")

}
# check r specification
.check_and_set_r <- function(r, stat) {

  if (is.null(r) && length(stat) == 1)
    r  <- 1

  if (!is.null(r) && r < 1)
    stop("r must be greater than 1")

  return(r)
}

