# source("~/Desktop/Research/BFF/R/FINAL_SUPPORT_hypergeometric.R")
# source("~/Desktop/Research/BFF/R/FINAL_FUNCTIONS_tau2.R")
# source("~/Desktop/Research/BFF/R/FINAL_FUNCTIONS_plotting.R")
# source("~/Desktop/Research/BFF/R/FINAL_support_functions.R")
# library(gsl)

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
log_F_frac = function(tau2, f, k, m, r)
{
  one = 1 / ((1 + tau2) ^ (k / 2 + r))

  a = k / 2 + r
  b = (k + m) / 2
  c = k / 2
  d = k * f * tau2 / ((1 + tau2) * (m + k * f))

  two = Gauss2F1(a, b, c, d)

  to_return = log(one) + log(two)


  return(to_return)
}


####################### backend implementation
backend_f = function(r,
                     f_stat,
                     n,
                     df1,
                     df2,
                     r1,
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
    if (!user_supplied_tau2)
      tau2 = get_linear_tau2(n = n, k = df1, w = effect_size)
    log_vals = unlist(lapply(
      tau2,
      f_val_r1,
      f_stat = f_stat,
      df1 = df1,
      df2 = df2
    ))
  }

  if (frac_r) {
    if (!user_supplied_tau2)
      tau2 = get_linear_tau2(n = n,
                             k = df1,
                             w = effect_size,
                             r = r)

    log_vals = unlist(lapply(
      tau2,
      log_F_frac,
      f = f_stat,
      k = df1,
      m = df2,
      r = r
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


################# F function user interaction

#' f_test_BFF
#'
#' f_test_BFF constructs BFFs based on the F test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#' Plot saved to working directory unless a full path is specified in the 'savename' variable of the function.
#'
#' @param f_stat F statistic
#' @param df1 first degree of freedom
#' @param df2 first degree of freedom
#' @param n sample size
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
#' fBFF = f_test_BFF(f_stat = 2.5, n = 50, df1 = 20, df2 = 48, save = FALSE)
#' f_test_BFF(f_stat = 2.5, n = 50, df1 = 20, df2 = 48, save = FALSE, tau2 = 0.5)
#' f_test_BFF(f_stat = 2.5, n = 50, df1 = 20, df2 = 48, save = FALSE, tau2 = c(0.5, 0.8))
#' f_test_BFF(f_stat = 2.5, n = 50, df1 = 20, df2 = 48, r = 2, save = FALSE)
#' f_test_BFF(f_stat = 2.5, n = 50, df1 = 20, df2 = 48, r = 2.5, save = FALSE)
#' f_test_BFF(f_stat=2.5, n = 50, df1 = 20, df2 = 48, maximize = TRUE)
#' f_test_BFF(f_stat=2.5, n = 50, df1 = 20, df2 = 48, maximize = TRUE, tau2 = 0.5)
#' f_test_BFF(f_stat=2.5, n = 50, df1 = 20, df2 = 48, maximize = TRUE, tau2 = c(0.5, 0.8))
#' fBFF$BFF_max_RMSE  # maximum BFF value
#' fBFF$max_RMSE      # effect size which maximizes the BFF value
#'
f_test_BFF = function(f_stat,
                      n,
                      df1,
                      df2,
                      savename = NULL,
                      maximize = FALSE,
                      r = 1,
                      tau2 = NULL,
                      save = TRUE,
                      xlab = NULL,
                      ylab = NULL,
                      main = NULL)

{
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
  results = backend_f(
    f_stat = f_stat,
    n = n,
    df1 = df1,
    df2 = df2,
    r = r,
    tau2 = tau2,
    r1 = r1
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
        backend_f,
        c(1, 20),
        tol = 0.001,
        f_stat = f_stat,
        n = n,
        df1 = df1,
        df2 = df2,
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
