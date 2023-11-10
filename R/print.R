#' @title Summarize BFF object
#'
#' @param object a BFF object
#' @param ... additional arguments (unused)
#'
#' @return prints summary of a BFF object.
#'
#' @seealso [z_test_BFF()], [t_test_BFF()], [chi2_test_BFF()], [f_test_BFF()]
#' @export
summary.BFF <- function(object, ...){
  print.BFF(object, ...)
}

#' @title Summarize BFF object
#'
#' @param x a BFF object
#' @param ... additional arguments (unused)
#'
#' @return prints summary of a BFF object.
#'
#' @seealso [z_test_BFF()], [t_test_BFF()], [chi2_test_BFF()], [f_test_BFF()]
#' @export
print.BFF = function(x, ...) {
  cat(paste0("\t\t", .test_type_name(x$test_type, x$one_sample)))
  cat("\n\n")
  cat(gettextf("%1$slog Bayes factor = %2$.2f\n", if(!x$omega_set) "maximized " else "", x$log_bf))
  cat(gettextf("%1$slog tau2 = %2$.2f\n", if(!x$omega_set) "maximized " else "", x$omega))
  cat(paste0("alternative = ", x$alternative))
}

.test_type_name = function(test_type, one_sample) {
  starting_strng = gettextf("Bayesian non-local %1$s %2$s",
                            if(one_sample) "one-sample" else "two-sample",
                            switch(test_type,
                                   "t_test"    = "t test",
                                   "z_test"    = "z test",
                                   "chi2_test" = "chi2 test",
                                   "f_test"    = "f test"))

}
