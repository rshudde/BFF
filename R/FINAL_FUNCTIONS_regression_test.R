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
backend_reg = function(r,
                     t_stat,
                     df,
                     n,
                     k,
                     one_sided = TRUE,
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
    tau2 = get_regression_tau2(n = n, k = k, w = effect_size, r = r)
  }


  if (one_sided) {
    log_vals = unlist(lapply(
      tau2,
      log_T_frac_onesided,
      r = r,
      v = df,
      t = t_stat
    ))
  } else {
    log_vals = unlist(lapply(
      tau2,
      log_T_frac,
      r = r,
      v = df,
      t = t_stat
    ))
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

maximize_reg = function(r,
                      t_stat,
                      df,
                      n,
                      k = TRUE,
                      one_sided = TRUE,
                      omega = NULL) {

  logbf = dcauchy(r)/(1-pcauchy(1))
  for (t in range(1, length(t_stat))) {
    logbf = logbf + backend_reg(r = r,
                              t_stat = t_stat[t],
                              df = df[t],
                              n = n[t],
                              k = k[t],
                              one_sided = one_sided,
                              omega = omega, # technically not used
                              tau2 = omega^2*n[t])
  }

  return(logbf)
}

################# T function user interaction

#' regression_test_BFF
#'
#' regression_test_BFF constructs BFFs based on the t test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#'
#' @param t_stat T statistic
#' @param n sample size (if one sample test)
#' @param k number of predictors
#' @param r r value
#' @param tau2 tau2 values (can be a single entry or a vector of values)
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
#'    \code{omega} \tab omega values tested, can be a single number or vector\cr
#' }
#' @export
#'
#' @examples
#' regBFF = regression_test_BFF(t_stat = 2.5, n = 50, k = 3)
#' regBFF$BFF_max_RMSE   # maximum BFF omega
#' regBFF$max_RMSE       # effect size which maximizes the BFF value
#'
regression_test_BFF = function(t_stat,
                      n,
                      k,
                      alternative = "two.sided",
                      r = NULL,
                      omega = NULL)

{

  # check alternative
  if (!alternative %in% c("two.sided", "less", "greater")) {
    stop("The alternative must be either 'two.sided', 'less', or 'greater'")
  }

  if (is.null(r) && length(t_stat) == 1) r = 1
  if (!is.null(r) && r < 1) {
    stop("r must be greater than 1")
  }

  # check that the correct lengths for everything is populated
  if (length(t_stat > 1)) {
    len_t = length(t_stat)
    if (length(n) != len_t) {
      stop("If providing a vector of t statistics, sample size must also be supplied as a vector of equal length")
    }
    if (length(k) != len_t) {
      stop("If providing a vector of t statistics, k (number of predictors) must also be supplied as a vector of equal length")
    }
  }

  df = n - k - 1

  for (i in df){
    if (i <= 1) {
      stop("Degrees of freedom must be greater than 1. If using a two sample test, n must be greater
         than 3, if using a one sample test, n must be greater than 2")
    }
  }

  used_alternative = alternative
  if (alternative == "less")
  {
    t_stat = -t_stat
    used_alternative = "greater"
  }

  # did user set
  omega_set = !is.null(omega)

  # should we maximize? If the t statistic is a vector and r is not provided, yes
  maximize = length(t_stat) > 1 && is.null(r)

  #####  same effect sizes for all tests
  omega_sequence = seq(0.01, 1, by = 0.01)

  ##### optimization logic
  if (maximize)
  {
    # set the "omega max" we are searching over. We are calling this omega
    # max because it is important to keep original value of omega for later
    if (is.null(omega)) {

      omega_max = omega_sequence
    } else {
      omega_max = omega
    }
    optimal_r = vector(length = length(omega_max))
    count = 1
    for (i in omega_max)
    {
      optimal_r[count] = optimize(
        maximize_reg,
        c(1, 20),
        tol = 0.001,
        t_stat = t_stat,
        df = df,
        n = n,
        k = k,
        one_sided = used_alternative == "greater",
        omega = i,
        maximum = TRUE
      )$maximum
      count = count + 1
    }
    maximized_values = as.data.frame(cbind(omega_max, optimal_r))

    r = optimal_r
    results = vector()
    for (i in 1:length(optimal_r)) {
      results[i] = maximize_reg(
        r = optimal_r[i],
        t_stat = t_stat,
        df = df,
        n = n,
        k = k,
        one_sidedkused_alternative == "greater",
        omega = omega_max[i]
      )
    }

  } else {
    results = backend_reg(
      t_stat = t_stat,
      n = n,
      df = df,
      r = r,
      omega = omega,
      k = k,
      one_sided = used_alternative == "greater"
    )
  }


  ###### return logic
  if (!omega_set) {
    log_bf         <- c(0, results)
    omega_sequence <- c(0, omega_sequence)
    idx_max        <- which.max(log_bf)
    this_log_bf    <- log_bf[idx_max]
    this_omega     <- omega_sequence[idx_max]
  }else{
    this_log_bf    <- results
    this_omega     <- omega
  }

  output = list(
    log_bf       = this_log_bf,
    omega        = this_omega,
    omega_set    = omega_set,
    alternative  = alternative,
    test_type    = "regression_test",
    generic_test = FALSE,
    r            = r, # r that is maximized or set by user
    input = list(
      t_stat = t_stat,
      df     = df,
      k = k
    )
  )
  if (!omega_set) {
    output$BFF = list(log_bf = log_bf, omega = omega_sequence)
  }

  class(output) = "BFF"
  return(output)
}




