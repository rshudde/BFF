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

  # print in favor of alternative
  cat(gettextf("%1$slog Bayes factor = %2$.2f\n", if(!x$omega_set) "maximized (in favor of alternative) " else "", x$log_bf_h1))
  if(x$generic_test){
    cat(gettextf("%1$s tau2 = %2$.2f\n", if(!x$omega_set) "maximized (in favor of alternative) " else "", x$omega_h1))
  }else{
    cat(gettextf("%1$somega = %2$.2f (%3$s)\n", if(!x$omega_set) "maximized (in favor of alternative) " else "", x$omega_h1, .test_effect_size_name(x$test_type)))
  }

  # print in favor of null, only if no omega is set
  if (!x$omega_set) {
    cat(gettextf("%1$slog Bayes factor = %2$.2f\n", if(!x$omega_set) "minimized (in favor of null for medium/large effect sizes) " else "", x$log_bf_h0))
    if(x$generic_test){
      cat(gettextf("%1$s tau2 = %2$.2f\n", if(!x$omega_set) "minimized (in favor of null for medium/large effect sizes) " else "", x$omega_h0))
    }else{
      cat(gettextf("%1$somega = %2$.2f (%3$s)\n", if(!x$omega_set) "minimized (in favor of null for medium/large effect sizes) " else "", x$omega_h0, .test_effect_size_name(x$test_type)))
    }
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
         "f_test"           = "Cohen's f",
         "chi2_test"        = "Cohen's w",
         "regression_test"  = "Cohen's d",
         "correlation_test" = "correlation coefficient")
}
