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

  to_return = c * (first_hypergeo + 2*y * gamma_term * second_hypergeo)
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
backend_t <- function(
    r,
    input,
    omega = NULL){

  # function can deal with a vector of t-statistics and vector of omegas
  # however, only one r is supported at a time
  if(length(r) != 1)
    stop("internal error: 'r' must be of length 1.")

  # compute tau2 from omega
  # if multiple omegas and t-stats are supplied, each element of tau2
  # corresponds a vector of tau2 for the corresponding t-statistics
  # i.e., tau2[omega][t-stat]
  tau2 <- lapply(omega, function(x){
    if(input$one_sample){
      tau2 <- get_one_sample_tau2(n = input$n, w = x, r = r)
    }else{
      tau2 <- get_two_sample_tau2(n1 = input$n1, n2 = input$n2, w = x, r = r)
    }
  })

  # compute log_BF
  log_BF <- sapply(tau2, function(x){
    sum(sapply(seq_along(input$t_stat), function(i){
      if(input$alternative == "greater"){
        log_T_frac_onesided(
          tau2 = x[i],
          r    = r,
          v    = input$df[i],
          t    = input$t_stat[i]
        )
      }else{
        log_T_frac(
          tau2 = x[i],
          r    = r,
          v    = input$df[i],
          t    = input$t_stat[i]
        )
      }
    }))
  })

  # check the results are finite
  if (!all(is.finite(log_BF)))
    warning(
      "Values entered produced non-finite numbers for some effect sizes.
      The most likely scenario is the evidence was so strongly in favor of the alternative that there was numeric overflow.
      Only effect sizes with non-NaN values are kept in the plots.
      Please contact the maintainer for more information."
    )

  return(log_BF)
}

maximize_t <- function(
    r,
    input,
    omega = NULL){

  logbf <- stats::dcauchy(r)/(1-stats::pcauchy(1))
  logbf <- logbf + sum(backend_t(
    r         = r,
    input     = input,
    omega     = omega))

  return(logbf)
}

################# T function user interaction

#' t_test_BFF
#'
#' t_test_BFF constructs BFFs based on the t test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#'
#' @param t_stat T statistic
#' @param n sample size (if one sample test)
#' @param one_sample is test one sided? Default is FALSE
#' @param alternative the alternative. options are "two.sided" or "less" or "greater"
#' @param n1 sample size of group one for two sample test. Must be provided if one_sample = FALSE
#' @param n2 sample size of group two for two sample test. Must be provided if one_sample = FALSE
#' @param r r value
#' @param omega standardized effect size. For the t-test, this is often called Cohen's d (can be a single entry or a vector of values)
#'
#' @return Returns an S3 object of class `BFF` (see `BFF.object` for details).
#' @export
#'
#' @examples
#' tBFF = t_test_BFF(t_stat = 2.5, n = 50, one_sample = TRUE)
#' tBFF
#' plot(tBFF)
t_test_BFF <- function(
    t_stat,
    n = NULL,
    n1 = NULL,
    n2 = NULL,
    one_sample = FALSE,
    alternative = "two.sided",
    r = NULL,
    omega = NULL,
    omega_sequence = if(is.null(omega)) seq(0.01, 1, by = 0.01)){


  ### input checks and processing
  r     <- .check_and_set_r(r, t_stat)
  input <- .process_input.t.test(t_stat, n, n1, n2, one_sample, alternative)

  ### computation
  # if only one t-statistic is present or if r is set, we can directly compute the BFF
  # otherwise, we need to find the r that maximizes the BFF
  if(length(input$t_stat) == 1 || !is.null(r)){

    results   <- backend_t(
      r         = r,
      input     = input,
      omega     = if(!is.null(omega)) omega else omega_sequence
    )

  }else{

    # compute optimal r for each omega
    optimal_r <- sapply(if(!is.null(omega)) omega else omega_sequence, function(x){
      stats::optimize(
        maximize_t,
        interval  = c(1, 20),
        tol       = 0.001,
        input     = input,
        omega     = x,
        maximum   = TRUE
      )$maximum
    })

    # use each r to compute BFF
    results <- sapply(seq_along(optimal_r), function(i){
      backend_t(
        r         = optimal_r[i],
        input     = input,
        omega     = if(!is.null(omega)) omega else omega_sequence[i]
      )
    })

  }

  ###### return logic
  if(is.null(omega)){
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
    omega_set    = !is.null(omega),
    test_type    = "t_test",
    generic_test = FALSE,
    r            = r, # r that is maximized or set by user
    input        = input
  )
  if(is.null(omega)){
    output$BFF = list(log_bf = log_bf, omega = omega_sequence)
  }

  class(output) = "BFF"
  return(output)
}


.process_input.t.test <- function(t_stat, n, n1, n2, one_sample, alternative){

  .check_alternative(alternative)

  # one vs. two-sample test processing
  if(one_sample){

    if(is.null(t_stat) || is.null(n))
      stop("Both t_stat and and n must be provided for one-sample (`one_sample = TRUE`) test.")
    if(length(t_stat) != length(n))
      stop("The input length of t_stat and n must be the same.")

    df <- n - 1
    .check_df(df, "(Total sample size must be greater than 2.)")
  }else{

    if(is.null(t_stat) || is.null(n1) || is.null(n2))
      stop("Both t_stat, n1, and n2 must be provided for two-sample (`one_sample = FALSE`) test.")
    if(length(t_stat) != length(n1) || length(t_stat) != length(n2))
      stop("The input length of t_stat, n1, and n2 must be the same.")

    df <- n1 + n2 - 2
    .check_df(df, "(Total sample size must be greater than 3.)")
  }

  # computation is implemented only for alternative = "two-sided" or "greater"
  # if lower, reverse the sign of t_stat, set alternative to "greater",
  # and remember that the original alternative was "less"
  if (alternative == "less"){
    t_stat      <- -t_stat
    alternative <- "greater"
    alternative.original <- "less"
  }else{
    alternative.original <- alternative
  }

  return(list(
    t_stat     = t_stat,
    n          = n,
    n1         = n1,
    n2         = n2,
    df         = df,
    one_sample = one_sample,
    alternative          = alternative,
    alternative.original = alternative.original
  ))
}
