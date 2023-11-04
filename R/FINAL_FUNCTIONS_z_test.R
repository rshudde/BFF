# rm(list = ls())
# source("~/Desktop/Research/BFF/R/FINAL_SUPPORT_hypergeometric.R")
# source("~/Desktop/Research/BFF/R/FINAL_FUNCTIONS_tau2.R")
# source("~/Desktop/Research/BFF/R/FINAL_FUNCTIONS_plotting.R")
# source("~/Desktop/Research/BFF/R/FINAL_support_functions.R")
################# Z functions if r is an integer and equal to 1
z_val_r1 = function(tau2, z_stat)
{
  term_one = (tau2 + 1) ^ (-3 / 2)
  term_two = 1 + tau2 * z_stat ^ 2 / (tau2 + 1)
  term_three = exp(tau2 * z_stat ^ 2 / (2 * (tau2 + 1)))

  to_return = term_one * term_two * term_three
  to_return = log(to_return)
  return(to_return)
}

################# Z functions if r is an integer and greater than 1
# sum_val_z = function(r, k, z, tau)
# {
#   # constructing third term of z statistic
#   one = choose(2 * r, 2 * k)
#   two = (tau * z / (1 + tau)) ^ (2 * k)
#   three = (tau / (1 + tau)) ^ (r - k)
#   four = double_factorial(2 * r - 2 * k - 1)
#
#   final_val = one * two * three * four
#   return(final_val)
# }
#
# sum_val_function_z = function(r, z, tau)
# {
#   val = 0
#   for (kk in 0:r)
#   {
#     val = val + sum_val_z(r = r,
#                           k = kk,
#                           z = z,
#                           tau = tau)
#   }
#   return(val)
# }
#
# log_Z = function(z, r, tau)
# {
#   num1 = 1
#   den1 = sqrt(1 + tau) * tau ^ (r) * double_factorial(2 * r - 1)
#   first_term = num1 / den1
#
#   num2 = tau * z ^ 2
#   den2 = 2 * (1 + tau)
#   second_term = exp(num2 / den2)
#
#   third_term = sum_val_function_z(r = r, z = z, tau = tau)
#
#   to_return = first_term * second_term * third_term
#
#   # log version
#   to_return = log(first_term) + log(second_term) + log(third_term)
#   return(to_return)
# }
################# Z functions if r is a f raction
log_Z_frac = function(tau2, z, r)
{
  c = 1 / ((1 + tau2) ^ (r + 1 / 2))

  third_term = tau2 * z ^ 2 / (2 * (1 + tau2))
  hyper_term = hypergeom1F1(r + 1 / 2, 1 / 2, third_term)$f # is this log or not already?

  to_return = c * hyper_term
  to_return = log(to_return)

  return(to_return)
}

log_Z_frac_onesided = function(tau2, z, r)
{
  tau = sqrt(tau2)
  c = 1 / ((1 + tau2) ^ (r + 1 / 2))
  y = tau * z / sqrt(2 * (1 + tau2))

  first_hyper = hypergeom1F1(r + 1 / 2, 1 / 2, y ^ 2)$f

  gamma_term = sterling_gamma(r + 1) / sterling_gamma(r + 1 / 2)

  second_hper = hypergeom1F1(r + 1, 3 / 2, y ^ 2)$f

  to_return = c * (first_hyper + 2 * y * gamma_term * second_hper)
  to_return = log(to_return)

  return(to_return)
}


####################### backend implementation
backend_z = function(r,
                     z_stat,
                     n,
                     one_sample = TRUE,
                     n1 = NULL,
                     n2 = NULL,
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
    log_vals = z_val_r1(tau2 = tau2, z_stat = z_stat)
  }

  if (frac_r) {
    if (one_sample)
    {
      if (!user_supplied_tau2)
        tau2 = get_one_sample_tau2(n = n, w = effect_size, r = r)
      log_vals = log_Z_frac_onesided(z = z_stat,
                                     r = r,
                                     tau2 = tau2)
    } else {
      if (!user_supplied_tau2)
        tau2 = get_two_sample_tau2(
          n1 = n1,
          n2 = n2,
          w = effect_size,
          r = r
        )
      log_vals = log_Z_frac(z = z_stat,
                            r = r,
                            tau2 = tau2)
    }
  }

  # stuff to return
  BFF = log_vals

  # check the results are finite
  if (!all(is.finite(BFF)))
  {
    stop(
      "Values entered produced non-finite numbers.
      The most likely scenario is the evidence was so strongly in favor of the
      alternative that there was numeric overflow. Please contact the maintainer for more information."
    )
  }

  return(BFF)

}


################# Z function user interaction

#' z_test_BFF
#'
#' z_test_BFF constructs BFFs based on the z test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#' Plot saved to working directory unless a full path is specified in the 'savename' variable of the function.
#'
#' @param z_stat z statistic
#' @param n sample size (if one sample test)
#' @param one_sample is test one sided? Default is TRUE
#' @param n1 sample size of group one for two sample test
#' @param n2 sample size of group two for two sample test
#' @param savename optional, filename for saving the pdf of the final plot
#' @param maximize Should the value of r be maximized? Default is FALSE. Only set to TRUE if analyzing multiple studies
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
#' zBFF = z_test_BFF(z_stat = 2.5, n = 50, save = FALSE)
#' z_test_BFF(z_stat = 2.5, n = 50, save = FALSE, tau2 = 0.5)
#' z_test_BFF(z_stat = 2.5, n = 50, save = FALSE, tau2 = c(0.5, 0.8))
#' z_test_BFF(z_stat = 2.5, n1 = 50, n2 = 35, one_sample = FALSE, save = FALSE) ##
#' z_test_BFF(z_stat = 2.5, n = 50, r = 2, save = FALSE)
#' z_test_BFF(z_stat = 2.5, r = 2, n1 = 50, n2 = 30, one_sample = FALSE, save = FALSE) ##
#' z_test_BFF(z_stat = 2.5, n = 50, r = 2.5, save = FALSE)
#' z_test_BFF(z_stat = 2.5, r = 2.5, n1 = 50, n2 = 30, one_sample = FALSE, save = FALSE) ##
#' z_test_BFF(z_stat=2.5, n = 50, maximize = TRUE)
#' z_test_BFF(z_stat=2.5, n = 50,  maximize = TRUE, tau2 = 0.5)
#' z_test_BFF(z_stat=2.5, n = 50,  maximize = TRUE, tau2 = c(0.5, 0.8))
#' zBFF$BFF_max_RMSE   # maximum BFF value
#' zBFF$max_RMSE       # effect size which maximizes the BFF value
#'
z_test_BFF = function(z_stat,
                      n = NULL,
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
  results = backend_z(
    z_stat = z_stat,
    n = n,
    one_sample = one_sample,
    r = r,
    tau2 = tau2,
    r1 = r1,
    n1 = n1,
    n2 = n2
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
        backend_z,
        c(1, 20),
        tol = 0.001,
        z_stat = z_stat,
        n = n,
        n1 = n1,
        n2 = n2,
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
