

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
                     n1 = NULL,
                     n2 = NULL,
                     n = NULL,
                     one_sample = TRUE,
                     one_sided = TRUE,
                     omega = NULL,
                     tau2 = NULL)

{
  # same effect sizes for all tests
  if (!is.null(omega))
  {
    omega_sequence = omega
  } else {
    omega_sequence = seq(0.01, 1, by = 0.01)
  }

  # user_supplied_omega = TRUE
  # if (is.null(omega))
  #   user_supplied_omega = FALSE

  log_vals = rep(0, length(omega_sequence))

  if (is.null(tau2))
  {
    if (one_sample)
    {
      tau2 = get_one_sample_tau2(n = n, w = omega_sequence, r = r)
    } else if (!one_sample)
      tau2 = get_two_sample_tau2(n1 = n1,
                                 n2 = n2,
                                 w = omega_sequence,
                                 r = r)
  }


  if (one_sided) {
    log_vals = unlist(lapply(
      tau2,
      log_Z_frac_onesided,
      r = r,
      z = z_stat
    ))
  } else {
    log_vals = unlist(lapply(
      tau2,
      log_Z_frac,
      r = r,
      z = z_stat
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

maximize_z = function(r,
                      z_stat,
                      n = NULL,
                      one_sample = TRUE,
                      one_sided = TRUE,
                      n1 = NULL,
                      n2 = NULL,
                      omega = NULL) {

  logbf = stats::dcauchy(r)/(1-stats::pcauchy(1))
  for (t in range(1, length(z_stat))) {
    logbf = logbf + backend_z(r = r,
                              z_stat = z_stat[t],
                              n = n[t],
                              one_sample = one_sample,
                              one_sided = one_sided,
                              n1 = n1[t],
                              n2 = n2[t],
                              omega = omega, # technically not used
                              tau2 = omega^2*n[t])
  }

  return(logbf)
}

################# T function user interaction

#' z_test_BFF
#'
#' z_test_BFF constructs BFFs based on the z test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#'
#' @param z_stat Z statistic
#' @param n sample size (if one sample test)
#' @param one_sample is test one sided? Default is FALSE
#' @param alternative the alternative. options are "two.sided" or "less" or "greater"
#' @param n1 sample size of group one for two sample test. Must be provided if one_sample = FALSE
#' @param n2 sample size of group two for two sample test. Must be provided if one_sample = FALSE
#' @param r r value
#' @param omega standardized effect size. For the z-test, this is often called Cohen's d (can be a single entry or a vector of values)
#'
#' @return Returns an S3 object of class `BFF` (see `BFF.object` for details).
#' @export
#'
#' @examples
#' zBFF = z_test_BFF(z_stat = 2.5, n = 50, one_sample = TRUE)
#' zBFF
#' plot(zBFF)

z_test_BFF = function(z_stat,
                      n = NULL,
                      one_sample = FALSE,
                      alternative = "two.sided",
                      n1 = NULL,
                      n2 = NULL,
                      r = NULL,
                      omega = NULL)

{

  ### input checks
  .check_alternative(alternative)
  r <- .check_and_set_r(r, z_stat)

  # check that the correct lengths for everything is populated
  if (length(z_stat > 1)) {
    len_t = length(z_stat)

    if (is.null(n)) {
      if (length(n1) != len_t || length(n2) != len_t) {
        stop("If providing a vector of t statistics, sample sizes must also be supplied as vectors of equal length")
      }
    } else {
      if (length(n) != len_t) {
        stop("If providing a vector of t statistics, sample size must also be supplied as a vector of equal length")
      }
    }
  }

  # check if one_sample is FALSE that n1 and n2 are provided
  if ((!one_sample) && (is.null(n1) || is.null(n1))) {
    stop("if one_sample is FALSE, both n1 and n2 must be provided")
  }

  used_alternative = alternative
  if (alternative == "less")
  {
    z_stat = -z_stat
    used_alternative = "greater"
  }

  # did user set
  omega_set = !is.null(omega)

  # should we maximize? If the t statistic is a vector and r is not provided, yes
  maximize = length(z_stat) > 1 && is.null(r)

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
      optimal_r[count] = stats::optimize(
        maximize_z,
        c(1, 20),
        tol = 0.001,
        z_stat = z_stat,
        n = n,
        one_sample = one_sample,
        one_sided = used_alternative == "greater",
        n1 = n1,
        n2 = n1,
        omega = i,
        maximum = TRUE
      )$maximum
      count = count + 1
    }
    maximized_values = as.data.frame(cbind(omega_max, optimal_r))

    r = optimal_r
    results = vector()
    for (i in 1:length(optimal_r)) {
      results[i] = maximize_z(
        r = optimal_r[i],
        z_stat = z_stat,
        n = n,
        one_sample = one_sample,
        one_sided =  used_alternative == "greater",
        n1 = n1,
        n2 = n2,
        omega = omega_max[i]
      )
    }

  } else {
    results = backend_z(
      z_stat = z_stat,
      n = n,
      r = r,
      n1 = n1,
      n2 = n2,
      omega = omega,
      one_sample = one_sample,
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
    test_type    = "z_test",
    generic_test = FALSE,
    r            = r, # r that is maximized or set by user
    input = list(
      z_stat      = z_stat,
      df          = NULL, # no df in a z test
      n1          = n1,
      n2          = n2,
      one_sample  = one_sample
    )
  )
  if (!omega_set) {
    output$BFF = list(log_bf = log_bf, omega = omega_sequence)
  }


  class(output) = "BFF"
  return(output)
}




