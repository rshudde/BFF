# source("~/Desktop/Research/BFF/R/FINAL_SUPPORT_hypergeometric.R")
# source("~/Desktop/Research/BFF/R/FINAL_FUNCTIONS_tau2.R")
# source("~/Desktop/Research/BFF/R/FINAL_FUNCTIONS_plotting.R")
# source("~/Desktop/Research/BFF/R/FINAL_support_functions.R")
# library(gsl)

################# chih2 functions if r is an integer and equal to 1
G_val_r1 = function(tau2, chi2_stat, df)
{
  BFF = (tau2 + 1) ^ (-df / 2 - 1) * (1 + tau2 * chi2_stat / (df * (tau2 + 1))) * exp(tau2 *
                                                                                        chi2_stat / (2 * (tau2 + 1)))

  to_return = log(BFF)
  return(to_return)
}


################# chi2 functions if r is an integer and greater than 1

# ################# Gamma functions / chi^2
# prod_val_g = function(k, r, n)
# {
#   one = k/2 + r - 1 - n
#   return(one)
# }
#
# prod_val_function_g = function(k, m, r)
# {
#   val = 1
#   if (m > 0)
#   {
#     for (nn in 0:(m-1))
#     {
#       val = val * prod_val_g(k=k, r=r, n=nn)
#     }
#   }
#   return(val)
# }
#
# sum_val_g = function(k, r, m, tau, x)
# {
#   one = choose(r, m)
#   two = prod_val_function_g(k = k, m=m, r=r)
#   three = (tau*x/(2*(1+tau)))^(r-m)
#   four = 1/((1+tau)^(k/2+r))
#
#   to_return = one*two*three*four
#   return(to_return)
# }
#
# sum_val_function_g = function(k, r, tau,x)
# {
#   val = 0
#   for (mm in 0:r)
#   {
#     val = val + sum_val_g(k=k, r = r, m =mm, tau=tau, x=x)
#   }
#
#   return(val)
# }
#
# log_G = function(tau, h, k, r)
# {
#   num1 = sterling_gamma(k/2)
#   den1 = sterling_gamma(k/2 + r)
#   first_term = num1/den1
#
#   second_term = exp(tau * h / (2*(1+tau)))
#
#   third_term = sum_val_function_g(k=k, r=r, tau=tau, x=h)
#
#   to_return = first_term * second_term * third_term
#
#   # log version
#   to_return = log(first_term) + log(second_term) + log(third_term)
#   return(to_return)
# }

################# T functions if r is a fraction
log_G_frac = function(tau2, h, k, r)
{
  tp1 = 1 + tau2
  one = 1 / (tp1 ^ (k / 2 + r))

  a = k / 2 + r
  b = k / 2
  c = tau2 * h / (2 * (1 + tau2))

  two = hypergeom1F1(a, b, c)$f

  to_return = log(one) + log(two)

  return(to_return)
}


####################### backend implementation
backend_chi2 = function(r,
                        chi2_stat,
                        n = NULL,
                        df = NULL,
                        pearsons = TRUE,
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
    if (pearsons) {
      if (!user_supplied_tau2)
        tau2 = get_count_tau2(n = n, k = df, w = effect_size)
    } else {
      if (!user_supplied_tau2)
        tau2 = get_LRT_tau2(n = n, k = df, w = effect_size)
    }
    log_vals = unlist(lapply(tau2, G_val_r1, chi2_stat = chi2_stat, df = df))
  }


  if (frac_r) {
    if (pearsons)
    {
      if (!user_supplied_tau2)
        tau2 = get_count_tau2(n = n,
                              w = effect_size,
                              k = df,
                              r = r)
    } else {
      if (!user_supplied_tau2)
        tau2 = get_LRT_tau2(n = n,
                            w = effect_size,
                            k = df,
                            r = r)
    }
    log_vals = unlist(lapply(
      tau2,
      log_G_frac,
      h = chi2_stat,
      r = r,
      k = df
    ))
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


################# T function user interaction

#' chi2_test_BFF
#'
#' chi2_test_BFF constructs BFFs based on the chi-squared test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#' Plot saved to working directory unless a full path is specified in the 'savename' variable of the function.
#'
#' @param chi2_stat chi^2 statistic
#' @param df degrees of freedom
#' @param n sample size
#' @param pearsons Is this a test of Pearson’s chi^2 test for goodness-of-fit? Default is TRUE. FALSE assumes a likelihood ratio test
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
#' chi2BFF = chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, save = FALSE)
#' chi2BFF = chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, save = FALSE, tau2 = 0.5)
#' chi2BFF = chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, save = FALSE, tau2 = c(0.5, 0.8))
#' chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, pearsons = FALSE, save = FALSE)
#' chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, r = 2, save = FALSE)
#' chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, r = 2, pearsons = FALSE, save = FALSE)
#' chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, r = 2.5, save = FALSE)
#' chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, r = 2.5, pearsons = FALSE, save = FALSE)
#' chi2_test_BFF(chi2_stat=2.5, n = 50, df = 49, maximize = TRUE)
#' chi2_test_BFF(chi2_stat=2.5, n = 50,  df = 49, maximize = TRUE, tau2 = 0.5)
#' chi2_test_BFF(chi2_stat=2.5, n = 50,  df = 49, maximize = TRUE, tau2 = c(0.5, 0.8))
#' chi2BFF$BFF_max_RMSE  # maximum BFF value
#' chi2BFF$max_RMSE      # effect size which maximizes the BFF
chi2_test_BFF = function(chi2_stat,
                         n = NULL,
                         df = NULL,
                         pearsons = TRUE,
                         savename = NULL,
                         maximize = FALSE,
                         r = 1,
                         tau2 = NULL,
                         save = TRUE,
                         xlab = NULL,
                         ylab = NULL,
                         main = NULL)

{
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
  results = backend_chi2(
    chi2_stat = chi2_stat,
    n = n,
    pearsons = pearsons,
    r = r,
    tau2 = tau2,
    r1 = r1,
    df = df
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
        backend_chi2,
        c(1, 20),
        tol = 0.001,
        chi2_stat = chi2_stat,
        n = n,
        df = df,
        pearsons = pearsons,
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
