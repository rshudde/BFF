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
print.BFF <- function(x, ...) {
  cat(paste0("\t", .test_type_name(x$test_type, x$input$one_sample)))
  cat("\n\n")
  cat(gettextf("%1$slog Bayes factor = %2$.2f\n", if(!x$omega_set) "maximized " else "", x$log_bf))
  if(x$generic_test){
    cat(gettextf("%1$s tau2 = %2$.2f\n", if(!x$omega_set) "maximized " else "", x$omega))
  }else{
    cat(gettextf("%1$somega = %2$.2f (%3$s)\n", if(!x$omega_set) "maximized " else "", x$omega, .test_effect_size_name(x$test_type)))
  }

  if(!is.null(x$input$alternative)) cat(paste0("alternative = ", x$input$alternative.original))
}

.test_type_name <- function(test_type, one_sample) {
  starting_strng = gettextf("Bayesian non-local %1$s%2$s",
                            if(!is.null(one_sample)) {if(one_sample) "one-sample " else "two-sample "} else "",
                            switch(test_type,
                                   "t_test"           = "t test",
                                   "z_test"           = "z test",
                                   "chi2_test"        = "chi2 test",
                                   "f_test"           = "f test",
                                   "regression_test"  = "regression test",
                                   "correlation_test" = "correlation_test"))

}
.test_effect_size_name <- function(test_type){
  switch(test_type,
         "t_test"           = "Cohen's d",
         "z_test"           = "Cohen's d",
         "chi2_test"        = "Standardized effect size",
         "f_test"           = "Standardized effect size",
         "regression_test"  = "Partial correlation coefficient (eta squared)",
         "correlation_test" = "correlation coefficient")
}
