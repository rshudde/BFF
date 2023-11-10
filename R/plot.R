#' @title Plot Bayes Factor Function
#'
#' @description Creates a Bayes factor function plot
#' of an BFF object. The BFF object needs to be specified
#' with \code{omega = NULL}.
#'
#' @param x a BFF object
#'
#'
#' @return either a ggplot2 object if \code{plot = TRUE} or a data.frame
#' with a Bayes factor function if \code{plot = FALSE}
#'
#' @seealso [z_test_BFF()], [t_test_BFF()], [chi2_test_BFF()], [f_test_BFF()]
#'
#' @export
plot.BFF = function(x, xlab = NULL, ylab = NULL, main = NULL, ...) {

  if (is.null(x$BFF)) stop("Bayes factor function can be plotted only if a specific tau2 is not user set")

  plot_BFF(
    effect_size = x$BFF$omega,
    BFF = x$BFF$log_bf,
    xlab = xlab,
    ylab = ylab,
    main = main,
    r = x$r
  )
}
