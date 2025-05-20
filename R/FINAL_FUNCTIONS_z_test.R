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

####################### backend implementation (SPL, default)
BFF_z_test = function(tau2, z_stat, r, two_sided)
{

  a = get_a(tau2=tau2, r=r)
  y = get_y_onesided_z_test(tau2=tau2, z=z_stat)


  if (two_sided) {
    final_BF = a * hypergeom1F1(r + 1/2, 1/2, tau2 * z_stat^2 / (2 * (1 + tau2)) )$f
  } else {
    first_hypergeo = hypergeom1F1(r +1/2, 1/2, y^2)$f
    second_hypergeo = hypergeom1F1(r + 1, 3/2, y^2)$f
    const = 2 * y * gamma_approx(r+1) / gamma_approx(r + 1/2)
    final_BF = a*(first_hypergeo + const*second_hypergeo)
  }
  to_return = log(final_BF)
  return(to_return)
}


####################### backend implementation
backend_z <- function(
    input,
    r,
    omega = NULL){

  # compute tau2 from omega
  # if multiple omegas and t-stats are supplied, each element of tau2
  # corresponds a vector of tau2 for the corresponding t-statistics
  # i.e., tau2[omega][t-stat]
  tau2 <- lapply(omega, function(x){
    if(input$one_sample){
      tau2 <- get_one_sample_tau2(n = input$n, w = x, r=r)
    }else{
      tau2 <- get_two_sample_tau2(n1 = input$n1, n2 = input$n2, w = x, r=r)
    }
  })

  # compute log_BF
  log_BF <- sapply(tau2, function(x){
    sum(sapply(seq_along(input$z_stat), function(i){
      BFF_z_test(
        tau2 = x[i],
        z_stat    = input$z_stat[i],
        r = r,
        two_sided = input$alternative == "two.sided"
      )
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


################# Z function user interaction

#' z_test_BFF
#'
#' z_test_BFF constructs BFFs based on the z test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#'
#' @param z_stat Z statistic
#' @param n sample size (if one sample test)
#' @param n1 sample size of group one for two sample test. Must be provided if one_sample = FALSE
#' @param n2 sample size of group two for two sample test. Must be provided if one_sample = FALSE
#' @param one_sample is test one sided? Default is FALSE
#' @param alternative the alternative. options are "two.sided" or "less" or "greater"
#' @param omega standardized effect size. For the z-test, this is often called Cohen's d (can be a single entry or a vector of values)
#' @param omega_sequence sequence of standardized effect sizes. If no omega is provided, omega_sequence is set to be seq(0.01, 1, by = 0.01)
#' @param r variable controlling dispersion of non-local priors. Default is 1. r must be >= 1
#'
#' @return Returns an S3 object of class `BFF` (see `BFF.object` for details).
#' @export
#'
#' @examples
#' zBFF = z_test_BFF(z_stat = 2.5, n = 50, one_sample = TRUE)
#' zBFF
#' plot(zBFF)

z_test_BFF <- function(
    z_stat,
    n = NULL,
    n1 = NULL,
    n2 = NULL,
    one_sample = FALSE,
    alternative = "two.sided",
    omega = NULL,
    omega_sequence = if(is.null(omega)) seq(0.01, 1, by = 0.01),
    r=1)

{
  ### input checks and processing
  input <- .process_input.z.test(z_stat, n, n1, n2, one_sample, alternative, r)

  ### computation
  # calculate BF
  results   <- backend_z(
    input     = input,
    r         = r,
    omega     = if(!is.null(omega)) omega else omega_sequence
  )

  ## compute minimum BFF for anything larger than small effect sizes
  if (is.null(omega)) {
    minimums = get_min_omega_bff(omega = omega_sequence, bff = results, cutoff = 0.2)
  }  else
  {
    minimums = c(NULL, NULL)
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
    log_bf_h1       = this_log_bf,
    omega_h1        = this_omega,
    log_bf_h0     = minimums[1],
    omega_h0      = minimums[2],
    omega_set    = !is.null(omega),
    test_type    = "z_test",
    generic_test = FALSE,
    r            = r,
    input        = input
  )
  if(is.null(omega)){
    output$BFF = list(log_bf = log_bf, omega = omega_sequence)
  }

  class(output) = "BFF"
  return(output)
}


.process_input.z.test <- function(z_stat, n, n1, n2, one_sample, alternative, r){


  if (r < 1)
    stop("r must be greater than or equal to 1")

  .check_alternative(alternative)

  # one vs. two-sample test processing
  if(one_sample){

    if(is.null(z_stat) || is.null(n))
      stop("Both z_stat and and n must be provided for one-sample (`one_sample = TRUE`) test.")
    if(length(z_stat) != length(n))
      stop("The input length of z_stat and n must be the same.")

    df <- n - 1
    .check_df(df, "(Total sample size must be greater than 2.)")
  }else{

    if(is.null(z_stat) || is.null(n1) || is.null(n2))
      stop("Both z_stat, n1, and n2 must be provided for two-sample (`one_sample = FALSE`) test.")
    if(length(z_stat) != length(n1) || length(z_stat) != length(n2))
      stop("The input length of z_stat, n1, and n2 must be the same.")

    df <- n1 + n2 - 2
    .check_df(df, "(Total sample size must be greater than 3.)")
  }

  # computation is implemented only for alternative = "two-sided" or "greater"
  # if lower, reverse the sign of z_stat, set alternative to "greater",
  # and remember that the original alternative was "less"
  if (alternative == "less"){
    z_stat      <- -z_stat
    alternative <- "greater"
    alternative.original <- "less"
  }else{
    alternative.original <- alternative
  }

  return(list(
    z_stat     = z_stat,
    n          = n,
    n1         = n1,
    n2         = n2,
    df         = df,
    one_sample = one_sample,
    alternative          = alternative,
    alternative.original = alternative.original
  ))
}

