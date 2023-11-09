
#' Check whether x is a BFF object
#'
#' @param x an object to test
#'
#' @return returns a boolean.
#'
#' @export
is.BFF <- function(x){
  return(inherits(x, "BFF"))
}
