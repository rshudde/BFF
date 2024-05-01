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
                     df,
                     n = NULL,
                     LRT = FALSE,
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

  if (is.null(tau2))
  {
    if (LRT)
    {
      tau2 = get_count_tau2(n = n, k = df, w = effect_size)
    } else {
      tau2 = get_LRT_tau2(n = n, k = df, w = effect_size)
    }
  }

  log_vals = unlist(lapply(
    tau2,
    log_G_frac,
    h = chi2_stat,
    r = r,
    k = df
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

maximize_chi2 = function(r,
                      chi2_stat,
                      df,
                      n = NULL,
                      LRT = FALSE,
                      omega = NULL) {

  logbf = stats::dcauchy(r)/(1-stats::pcauchy(1))
  for (t in range(1, length(chi2_stat))) {
    logbf = logbf + backend_chi2(r = r,
                              chi2_stat = chi2_stat[t],
                              df = df[t],
                              n = n[t],
                              LRT = LRT,
                              omega = omega, # technically not used
                              tau2 = omega^2*n[t])
  }

  return(logbf)
}

################# T function user interaction

#' chi2_test_BFF
#'
#' chi2_test_BFF constructs BFFs based on the t test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#'
#' @param chi2_stat chi-square statistic
#' @param n sample size (if one sample test)
#' @param df degrees of freedom
#' @param LRT should LRT be performed? Default is FALSE
#' @param r r value
#' @param omega standardized effect size. (can be a single entry or a vector of values)
#'
#' @return Returns an S3 object of class `BFF` (see `BFF.object` for details).
#' @export
#'
#' @examples
#' chi2BFF = chi2_test_BFF(chi2_stat = 7.5, n = 25, df = 23)
#' chi2BFF
#' plot(chi2BFF)
#'
chi2_test_BFF = function(chi2_stat,
                      n,
                      df,
                      LRT = FALSE,
                      r = NULL,
                      omega = NULL)

{

  if (is.null(r) && length(chi2_stat) == 1) r = 1
  if (!is.null(r) && r < 1) {
    stop("r must be greater than 1")
  }

  # check that the correct lengths for everything is populated
  if (length(chi2_stat > 1)) {
    if (length(chi2_stat) != length(n)) {
      stop("If providing a vector of t statistics, sample size must also be supplied as vectors of equal length")
    }
    if (length(chi2_stat) != length(df)) {
      stop("If providing a vector of t statistics, degrees of freedom must also be supplied as vectors of equal length")
    }
  }

  for (k in df){
    if (k <= 1) {
      stop("Degrees of freedom must be greater than 1. If using a two sample test, n must be greater
         than 3, if using a one sample test, n must be greater than 2")
    }
  }

  # did user set
  omega_set = !is.null(omega)

  # should we maximize? If the t statistic is a vector and r is not provided, yes
  maximize = length(chi2_stat) > 1 && is.null(r)

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
        maximize_chi2,
        c(1, 20),
        tol = 0.001,
        chi2_stat = chi2_stat,
        df = df,
        n = n,
        LRT = LRT,
        omega = i,
        maximum = TRUE
      )$maximum
      count = count + 1
    }
    maximized_values = as.data.frame(cbind(omega_max, optimal_r))

    r = optimal_r
    results = vector()
    for (i in 1:length(optimal_r)) {
      results[i] = maximize_chi2(
        r = optimal_r[i],
        chi2_stat = chi2_stat,
        df = df,
        n = n,
        LRT = LRT,
        omega = omega_max[i]
      )
    }

  } else {
    results = backend_chi2(
      chi2_stat = chi2_stat,
      n = n,
      df = df,
      r = r,
      LRT = LRT,
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
    LRT = LRT,
    test_type = "chi2_test",
    generic_test = FALSE,
    r = r, # r that is maximized or set by user
    input = list(
      chi2_stat = chi2_stat,
      df     = df
    )
  )
  if (!omega_set) {
    output$BFF = list(log_bf = results, omega = effect_size)
  }

  class(output) = "BFF"
  return(output)
}




