## S3 class support functions
print.BFF = function(x, ...) {
  cat(paste0("\t\t", .test_type_name(x$test_type, x$one_sample)))
  cat("\n\n")
  cat(gettextf("%1$slog Bayes factor = %2$.2f\n", if(!x$omega_set) "maximized " else "", x$log_bf))
  cat(gettextf("%1$slog tau2 = %2$.2f\n", if(!x$omega_set) "maximized " else "", x$tau2))
  cat(paste0("alternative = ", x$alternative))
}

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

.test_type_name = function(test_type, one_sample) {
  starting_strng = gettextf("Bayesian non-local %1$s %2$s",
                            if(one_sample) "one-sample" else "two-sample",
                            switch(test_type,
                                   "t_test" = "t test",
                                   "z_test" = "z test",
                                   "chi2_test" = "chi2 test",
                                   "f_test" = "f test"))

}






## math support functions
double_factorial_even = function(n) {
  first = 2^(n/2)
  second = factorial(n/2)
  to_return = first * second
  return(to_return)
}

double_factorial_odd = function(n) {
  n1 = n+1
  numerator = factorial(n1)
  first = 2^(n1/2)
  second = factorial(n1/2)
  denomonator = first*second

  to_return = numerator / denomonator
  return(to_return)
}

# double factorial expression
double_factorial = function(n) {
  if (n %% 2 == 1) {
    to_return = double_factorial_odd(n)
  } else {
    to_return = double_factorial_even(n)
  }

  return(to_return)
}

# approximation of gamma for large n
sterling_gamma = function(n)
{
  # if (n<7)
  # {
  #   const = gamma(n)
  # } else {
  #   const = n*log(n) - n
  # }

  return(exp(lgamma(n)))

  # return(const)
}
