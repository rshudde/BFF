################# F functions if r is an integer and equal to 1
f_val_r1 = function(tau2, f_stat, df1, df2)
{
  v = df2 * (tau2 + 1)
  term_one = (tau2 + 1) ^ (-df1 / 2 - 1)
  term_two = (1 + df1 * f_stat / df2) / (1 + df1 * f_stat / v)
  term_three = 1 + (df1 + df2) * tau2 * f_stat / (v * (1 + df1 * f_stat /
                                                         v))

  BFF = term_one * term_two * term_three

  to_return = log(BFF)
  return(to_return)
}


# ################# F functions if r is an integer and greater than 1
# prod_val_f = function(k, r, n)
# {
#   term_one = k / 2 + r - 1 - n
#
#   return(term_one)
# }
#
# prod_val_function_f = function(i, k, r)
# {
#   val = 1
#   if (i > 0)
#   {
#     for (nn in 0:(i - 1))
#     {
#       val = val * prod_val_f(k = k, r = r, n = nn)
#     }
#   }
#   return(val)
# }
#
# sum_val_f = function(i, k, f, m, r, tau)
# {
#   w = k * f / m
#   c = 1 + tau
#   a = tau / c
#
#   one = choose(r, i)
#   two = prod_val_function_f(i = i, k = k, r = r)
#   three = (a * w / (1 + w / c)) ^ (r - i)
#   four = gamma((k + m) / 2 + r - i)
#
#   to_return = one * two * three * four
#   return(to_return)
# }
#
#
# sum_val_function_f = function(k, f, m, r, tau)
# {
#   val = 0
#   for (ii in 0:r)
#   {
#     val = val + sum_val_f(
#       i = ii,
#       k = k,
#       f = f,
#       m = m,
#       r = r,
#       tau = tau
#     )
#   }
#
#   return(val)
# }
#
# log_F = function(tau, f, k, m, r)
# {
#   num1 = sterling_gamma(k / 2)
#   den1 = sterling_gamma(k / 2 + r)
#   den2 = sterling_gamma(k / 2 + m / 2)
#   first_term = num1 / (den1 * den2)
#
#   second_term = (1 + tau) ^ (-k / 2 - r)
#
#   w = k * f / m
#   num3 = 1 + w
#   den3 = 1 + w / (1 + tau)
#   third_term = (num3 / den3) ^ ((k + m) / 2)
#
#   fourth_term = sum_val_function_f(
#     k = k,
#     f = f,
#     m = m,
#     r = r,
#     tau = tau
#   )
#
#   to_return = first_term * second_term * third_term * fourth_term
#
#   # log version
#   to_return = log(first_term) + log(second_term) + log(third_term) + log(fourth_term)
#   return(to_return)
# }

################# F functions if r is a fraction
log_F_frac = function(tau2, f_stat, k, m, r)
{
  one = 1 / ((1 + tau2) ^ (k / 2 + r))

  a = k / 2 + r
  b = (k + m) / 2
  c = k / 2
  d = k * f_stat * tau2 / ((1 + tau2) * (m + k * f_stat))

  two = Gauss2F1(a, b, c, d)

  to_return = log(one) + log(two)


  return(to_return)
}


####################### backend implementation
backend_f = function(r,
                     f_stat,
                     df1,
                     df2,
                     n = NULL,
                     omega = NULL,
                     tau2 = NULL)

{
  # same effect sizes for all tests
  if (!is.null(omega))
  {
    effect_size = omega
  } else {
    effect_size = seq(0.01, 1, by = 0.01)
  }

  # user_supplied_omega = TRUE
  # if (is.null(omega))
  #   user_supplied_omega = FALSE

  log_vals = rep(0, length(effect_size))

  tau2 = get_linear_tau2(n = n, w = effect_size, k = df1, r = r)

  log_vals = unlist(lapply(
    tau2,
    log_F_frac,
    f_stat = f_stat,
    r = r,
    k = df1,
    m = df2
  ))

  # stuff to return
  BFF = log_vals

  # check the results are finite
  if (!all(is.finite(BFF)))
  {
    warning(
      "Values entered produced non-finite numbers for some effect sizes.
      The most likely scenario is the evidence was so strongly in favor of the alternative that there was numeric overflow.
      Only effect sizes with non-NaN values are kept in the plots.
      Please contact the maintainer for more information."
    )
  }

  return(BFF)
}

maximize_f = function(r,
                      f_stat,
                      df1,
                      df2,
                      n = NULL,
                      omega = NULL) {

  logbf = stats::dcauchy(r)/(1-stats::pcauchy(1))
  for (t in range(1, length(f_stat))) {
    logbf = logbf + backend_f(r = r,
                              f_stat = f_stat[t],
                              n = n[t],
                              df1 = df1[t],
                              df2 = df2[t],
                              omega = omega, # technically not used
                              tau2 = omega^2*n[t])
  }

  return(logbf)
}

################# T function user interaction

#' f_test_BFF
#'
#' f_test_BFF constructs BFFs based on the t test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#'
#' @param f_stat T statistic
#' @param n sample size (if one sample test)
#' @param df1 sample size of group one for two sample test.
#' @param df2 sample size of group two for two sample test
#' @param r r value
#' @param omega omega values (can be a single entry or a vector of values)
#'
#' @return Returns an S3 object with Bayes Factor function results.
#'  \tabular{ll}{
#'    \code{BFF} \tab the object containing the log_bf (log bayes factor values) and corresponding omega values \cr
#'    \tab \cr
#'    \code{input} \tab the object containing the input values \cr
#'    \tab \cr
#'    \code{log_bf} \tab maximized bayes factor\cr
#'    \tab \cr
#'    \code{omega} \tab corresponding omega value for maximized bayes factor\cr
#'    \tab \cr
#'    \code{omega_set} \tab was an omega value provided?\cr
#'    \tab \cr
#'    \code{r} \tab r value (default is 1 if not provided by user) \cr
#'    \tab \cr
#'    \code{test_type} \tab type of BFF test\cr
#'    \tab \cr
#'    \code{generic_test} \tab FALSE \cr
#'    \tab \cr
#' }
#' @export
#'
#' @examples
#' fBFF = f_test_BFF(f_stat = 2.5, n = 50, df1 = 25, df2 = 48)
#' fBFF
#' plot(fBFF)
#'
f_test_BFF = function(f_stat,
                      n,
                      df1,
                      df2,
                      r = NULL,
                      omega = NULL)

{

  if (is.null(r) && length(f_stat) == 1) r = 1
  if (!is.null(r) && r < 1) {
    stop("r must be greater than 1")
  }

  # check that the correct lengths for everything is populated
  if (length(f_stat > 1)) {
    len_t = length(f_stat)

    if (is.null(n)) {
      if (length(df1) != len_t || length(df1) != len_t) {
        stop("If providing a vector of t statistics, degrees of freedom must also be supplied as vectors of equal length")
      }
    } else {
      if (length(n) != len_t) {
        stop("If providing a vector of t statistics, sample size must also be supplied as a vector of equal length")
      }
    }
  }

  # did user set
  omega_set = !is.null(omega)

  # should we maximize? If the t statistic is a vector and r is not provided, yes
  maximize = length(f_stat) > 1 && is.null(r)

  #####  same effect sizes for all tests
  effect_size = seq(0.01, 1, by = 0.01)

  ##### optimization logic
  if (maximize)
  {
    # set the "omega max" we are searching over. We are calling this omega
    # max because it is important to keep original value of omega for later
    if (is.null(omega)) {

      omega_max = effect_size
    } else {
      omega_max = omega
    }
    optimal_r = vector(length = length(omega_max))
    count = 1
    for (i in omega_max)
    {
      optimal_r[count] = stats::optimize(
        maximize_f,
        c(1, 20),
        tol = 0.001,
        f_stat = f_stat,
        df1 = df1,
        df2 = df2,
        n = n,
        omega = i,
        maximum = TRUE
      )$maximum
      count = count + 1
    }
    maximized_values = as.data.frame(cbind(omega_max, optimal_r))

    r = optimal_r
    results = vector()
    for (i in 1:length(optimal_r)) {
      results[i] = maximize_f(
        r = optimal_r[i],
        f_stat = f_stat,
        df1 = df1,
        df2 = df2,
        n = n,
        omega = omega_max[i]
      )
    }

  } else {
    results = backend_f(
      f_stat = f_stat,
      n = n,
      df1 = df1,
      df2 = df2,
      r = r,
      omega = omega
    )
  }


  ###### return logic
  BFF = results
  effect_size = effect_size
  idx_max = which.max(BFF)
  BFF_max_RMSE = BFF[idx_max]
  max_RMSE = effect_size[idx_max]

  output = list(
    log_bf = BFF_max_RMSE,
    omega = max_RMSE,
    omega_set = omega_set,
    test_type = "f_test",
    generic_test = FALSE,
    r = r, # r that is maximized or set by user
    input = list(
      f_stat = f_stat,
      df1     = df1,
      df2 = df2
    )
  )
  if (!omega_set) {
    output$BFF = list(log_bf = results, omega = effect_size)
  }

  class(output) = "BFF"
  return(output)
}




