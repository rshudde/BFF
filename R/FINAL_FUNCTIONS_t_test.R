# source("~/Desktop/Research/BFF/R/FINAL_SUPPORT_hypergeometric.R")
# source("~/Desktop/Research/BFF/R/FINAL_FUNCTIONS_tau2.R")
# source("~/Desktop/Research/BFF/R/FINAL_FUNCTIONS_plotting.R")
# source("~/Desktop/Research/BFF/R/FINAL_support_functions.R")
# library(gsl)

################# T functions if r is an integer and equal to 1
t_val_r1 = function(tau2, t_stat, df)
{
  r = 1 + t_stat ^ 2 / df
  s = 1 + t_stat ^ 2 / (df * (1 + tau2))
  q = tau2 * (df + 1) / (df * (1 + tau2))

  BF = (tau2 + 1) ^ (-3 / 2) * (r / s) ^ ((df + 1) / 2) * (1 + q * t_stat ^
                                                             2 / s)

  to_return = log(BF)
  return(to_return)
}


################# T functions if r is an integer and greater than 1
#
# get_w = function(tau, v, t = t)
# {
#   num = 1 + t ^ 2 / (v * (1 + tau))
#   den = 1 + t ^ 2 / v
#
#   to_return = num / den
#   return(to_return)
# }
#
# sum_val_t = function(r, m, tau, t, v, w)
# {
#   one = choose(2 * r, 2 * m)
#   two = (2 * tau * t ^ 2 / ((t ^ 2 + v) * (tau + 1) * w)) ^ m
#   three = sterling_gamma((v + 2 * m + 1) / 2) * double_factorial(2 * r - 2 *
#                                                                    m - 1)
#
#   to_return = one * two * three
#   return(to_return)
# }
#
# sum_val_function_t = function(r, tau, t, v, w)
# {
#   val = 0
#   for (mm in 0:r)
#   {
#     val = val + sum_val_t(
#       r = r,
#       m = mm,
#       tau = tau,
#       t = t,
#       v = v,
#       w = w
#     )
#   }
#   return(val)
# }
#
# log_T = function(t, r, tau, v)
# {
#   w = get_w(tau = tau, v = v, t = t)
#   num1 = 1
#   den1 = double_factorial(2 * r - 1) * (1 + tau) ^ (r + 1 / 2) * sterling_gamma((v +
#                                                                                    1) / 2) * w ^ ((v + 1) / 2)
#   first_term = num1 / den1
#
#   second_term = sum_val_function_t(
#     r = r,
#     tau = tau,
#     t = t,
#     v = v,
#     w = w
#   )
#
#   to_return = first_term * second_term
#
#   # log_version
#   to_return = log(first_term) + log(second_term)
#
#   return(to_return)
# }

################# T functions if r is a fraction
log_T_frac = function(tau2, t, v, r)
{
  tp1 = 1 + tau2 # one plus tau^2
  tau = sqrt(tau2)
  c = 1 / (tp1 ^ (r + 1 / 2)) # c
  y = tau * t / sqrt((t ^ 2 + v) * tp1)

  a1 = (v + 1) / 2
  b1 = r + 1 / 2
  c1 = 1 / 2
  first_hypergeo = Gauss2F1(a1, b1, c1, y ^ 2)

  four = sterling_gamma(v / 2 + 1) * sterling_gamma(r + 1)
  five = sterling_gamma((v + 1) / 2) * sterling_gamma(r + 1 / 2)
  gamma_term = four / five

  aa = v / 2 + 1
  bb = r + 1
  cc = 3 / 2
  second_hypergeo = Gauss2F1(aa, bb, cc, y ^ 2)

  to_return = c * (first_hypergeo + y * gamma_term * second_hypergeo)
  to_return = log(to_return)

  return(to_return)
}

log_T_frac_onesided = function(tau2, t, v, r)
{
  tp1 = 1 + tau2
  tau = sqrt(tau2)
  c = 1 / (tp1 ^ (r + 1 / 2)) # c
  y = tau * t / sqrt((t ^ 2 + v) * tp1)

  a1 = (v + 1) / 2
  b1 = r + 1 / 2
  c1 = 1 / 2
  first_hypergeo = Gauss2F1(a1, b1, c1, y ^ 2)

  four = sterling_gamma(v / 2 + 1) * sterling_gamma(r + 1)
  five = sterling_gamma((v + 1) / 2) * sterling_gamma(r + 1 / 2)
  gamma_term = four / five

  aa = v / 2 + 1
  bb = r + 1
  cc = 3 / 2
  second_hypergeo = Gauss2F1(aa, bb, cc, y ^ 2)

  to_return = c * (first_hypergeo + 2 * y * gamma_term * second_hypergeo)
  to_return = log(to_return)

  return(to_return)
}

####################### backend implementation
backend_t = function(r,
                     t_stat,
                     n = NULL,
                     df = NULL,
                     one_sample = TRUE,
                     n1 = NULL,
                     n2 = NULL,
                     savename = NULL,
                     r1 = FALSE,
                     tau2 = NULL)

{
  # same effect sizes for all tests
  effect_size = seq(0.01, 1, by = 0.01)

  user_supplied_tau2 = TRUE
  if (is.null(tau2))
    user_supplied_tau2 = FALSE

  r1 = r1
  frac_r = !r1

  log_vals = rep(0, length(effect_size))

  if (r1) {
    if (one_sample)
    {
      if (!user_supplied_tau2)
        tau2 = get_one_sample_tau2(n = n, w = effect_size)
    } else {
      if (!user_supplied_tau2)
        tau2 = get_two_sample_tau2(n1 = n1, n2 = n2, w = effect_size)
    }
    log_vals = unlist(lapply(tau2, t_val_r1, t_stat = t_stat, df = df))
  }

  if (frac_r) {
    if (one_sample)
    {
      if (!user_supplied_tau2)
        tau2 = get_one_sample_tau2(n = n, w = effect_size, r = r)
      log_vals = unlist(lapply(
        tau2,
        log_T_frac,
        r = r,
        v = df,
        t = t_stat
      ))
    } else {
      if (!user_supplied_tau2)
        tau2 = get_two_sample_tau2(
          n1 = n1,
          n2 = n2,
          w = effect_size,
          r = r
        )

      log_vals = unlist(lapply(
        tau2,
        log_T_frac_onesided,
        r = r,
        v = df,
        t = t_stat
      ))
    }
  }

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

################# T function user interaction

#' t_test_BFF
#'
#' t_test_BFF constructs BFFs based on the t test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#' Plot saved to working directory unless a full path is specified in the 'savename' variable of the function.
#'
#' @param t_stat T statistic
#' @param df degrees of freedom
#' @param n sample size (if one sample test)
#' @param one_sample is test one sided? Default is TRUE
#' @param n1 sample size of group one for two sample test
#' @param n2 sample size of group two for two sample test
#' @param savename optional, filename for saving the pdf of the final plot
#' @param maximize should the function be maximzied over all possible r values? Default is FALSE. Only set to TRUE if analyzing multiple studies
#' @param r r value
#' @param tau2 tau2 values (can be a single entry or a vector of values)
#' @param save should a copy of the plot be saved?
#' @param xlab optional, x label for plot
#' @param ylab optional, y label for plot
#' @param main optional, main label for plot
#'
#' @return Returns Bayes factor function results
#'  \tabular{ll}{
#'    \code{BFF} \tab The log of the Bayes Factor Function values \cr
#'    \tab \cr
#'    \code{effect_size} \tab Effect sizes tested (seq(0, 1, by = 0.01)) \cr
#'    \tab \cr
#'    \code{BFF_max_RMSE} \tab Maximum BFF value \cr
#'    \tab \cr
#'    \code{max_RMSE} \tab Effect size that maximizes BFF\cr
#'    \tab \cr
#'    \code{tau2} \tab tau^2 values tested\cr
#' }
#' @export
#'
#' @examples
#' tBFF = t_test_BFF(t_stat = 2.5, n = 50, df = 49, save = FALSE)
#' t_test_BFF(t_stat = 2.5, n = 50, df = 49, save = FALSE, tau2 = 0.5)
#' t_test_BFF(t_stat = 2.5, n = 50, df = 49, save = FALSE, tau2 = c(0.5, 0.2))
#' t_test_BFF(t_stat = 2.5, n1 = 50, n2 = 40, df = 88, save = FALSE, one_sample = FALSE)
#' t_test_BFF(t_stat = 2.5, n = 50, r = 2, df = 49, save = FALSE)
#' t_test_BFF(t_stat = 2.5, r = 2, n1 = 50, n2 = 30, df = 78, one_sample = FALSE, save = FALSE)
#' t_test_BFF(t_stat = 2.5, n = 50, r = 2.5, df = 49, save = FALSE)
#' t_test_BFF(t_stat=2.5, r = 2.5, n1 = 50, n2 = 30, df = 78, one_sample = FALSE, save=FALSE)
#' t_test_BFF(t_stat = 2.5, n = 50, df = 49, save = FALSE, maximize = TRUE)
#' t_test_BFF(t_stat = 2.5, n = 50, df = 49, save = FALSE, maximize = TRUE, tau2 = 0.5)
#' t_test_BFF(t_stat = 2.5, n = 50, df = 49, save = FALSE, maximize = TRUE, tau2 = c(0.5, 0.8))
#' tBFF$BFF_max_RMSE   # maximum BFF value
#' tBFF$max_RMSE       # effect size which maximizes the BFF value
#'
t_test_BFF = function(t_stat,
                      n = NULL,
                      df = NULL,
                      one_sample = TRUE,
                      n1 = NULL,
                      n2 = NULL,
                      savename = NULL,
                      maximize = FALSE,
                      r = 1,
                      tau2 = NULL,
                      save = TRUE,
                      xlab = NULL,
                      ylab = NULL,
                      main = NULL)

{
  if (is.null(n) &
      (is.null(n1) &
       is.null(n2)))
    stop("Either n or n1 and n2 is required")

  if (is.null(df))
    stop("df is required")

  #####  same effect sizes for all tests
  effect_size = seq(0.01, 1, by = 0.01)


  ##### is tau2 supplied as an argument?
  user_supplied_tau2 = TRUE
  if (is.null(tau2))
  {
    user_supplied_tau2 = FALSE
  }

  #####  call results
  r1 = FALSE
  if (r == 1)
    r1 = TRUE
  results = backend_t(
    t_stat = t_stat,
    n = n,
    df = df,
    r = r,
    n1 = n1,
    n2 = n2,
    tau2 = tau2,
    r1 = r1,
    one_sample = one_sample
  )

  #####  plotting if tau2 is not specified
  if (!user_supplied_tau2 && !maximize) {
    bff_plot = c()
    bff_plot[[1]] = results

    plot_BFF(
      effect_size = effect_size,
      BFF = bff_plot,
      save = save,
      savename = savename,
      xlab = xlab,
      ylab = ylab,
      main = main,
      r = r
    )
  }

  ##### optimzation logic
  if (maximize)
  {
    if (is.null(tau2))
      tau2 = seq(0, 1, 0.1)
    optimal_r = vector(length = length(tau2))
    count = 1
    for (i in tau2)
    {
      optimal_r[count] = optimize(
        backend_t,
        c(1, 20),
        tol = 0.001,
        t_stat = t_stat,
        n = n,
        n1 = n1,
        n2 = n2,
        df = df,
        one_sample = one_sample,
        r1 = FALSE,
        tau2 = i,
        maximum = TRUE
      )$maximum
      count = count + 1
    }
    maximized_values = as.data.frame(cbind(tau2, optimal_r))
  }


  ###### return logic
  BFF = results
  effect_size = effect_size
  idx_max = which.max(BFF)
  BFF_max_RMSE = BFF[idx_max]
  max_RMSE = effect_size[idx_max]

  if (maximize) {
    print(
      "The maximum r value for each specified tau2 is given. Re-run the test with the desired r to generate plots and get the BFF value."
    )
    to_return = maximized_values
  } else if (user_supplied_tau2) {
    to_return = list(BFF = BFF,
                     tau2 = tau2)
  } else {
    to_return = list(
      log_BFF = BFF,
      effect_size = effect_size,
      log_BFF_max_RMSE = BFF_max_RMSE,
      max_RMSE = max_RMSE
    )

  }
  return(to_return)

}
