################# T functions if r is an integer and equal to 1
reg_t_val_r1 = function(tau2, t_stat, df)
{
  r = 1 + t_stat ^ 2 / df
  s = 1 + t_stat ^ 2 / (df * (1 + tau2))
  q = tau2 * (df + 1) / (df * (1 + tau2))

  BF = (tau2 + 1) ^ (-3 / 2) * (r / s) ^ ((df + 1) / 2) * (1 + q * t_stat ^
                                                             2 / s)

  to_return = log(BF)
  return(to_return)
}

####################### backend implementation (SPL, default)
BFF_reg_test = function(tau2, t_stat, df, r, two_sided)
{

  a = get_a(tau2=tau2, r=r)
  c = get_c(tau2=tau2, df=df, r=r)
  y = get_y_t_test(tau2=tau2, t=t_stat, df=df)

  first_hypergeo = Gauss2F1((df+1)/2, r + 1/2, 1/2, y^2)
  second_hypergeo = Gauss2F1(df/2 + 1, r + 1, 3/2, y^2)

  const = ifelse(two_sided, 1, 2)

  final_BF = a*(first_hypergeo + const*c*y*second_hypergeo)
  to_return = log(final_BF)
  return(to_return)
}

####################### backend implementation
backend_reg <- function(
    input,
    r,
    omega = NULL){

  # compute tau2 from omega
  # if multiple omegas and t-stats are supplied, each element of tau2
  # corresponds a vector of tau2 for the corresponding t-statistics
  # i.e., tau2[omega][t-stat]
  tau2 <- lapply(omega, function(x){
      tau2 <- get_regression_tau2(n = input$n, k = input$k, w = x, r = r)
  })

  # compute log_BF
  log_BF <- sapply(tau2, function(x){
    sum(sapply(seq_along(input$t_stat), function(i){
      BFF_reg_test(
        tau2 = x[i],
        t_stat    = input$t_stat[i],
        df = input$df[i],
        r = r,
        two_sided = input$alternative.original == "two.sided"
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


################# T function user interaction

#' regression_test_BFF
#'
#' regression_test_BFF constructs BFFs based on the t test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#'
#' @param t_stat T statistic
#' @param alternative is the alternative a one.sided or two.sided test? default is two.sided
#' @param n sample size (if one sample test)
#' @param k number of predictors
#' @param omega standadized effect size. For the regression test, this is also known as Cohen's f^@ (can be a single entry or a vector of values)
#' @param omega_sequence sequence of standardized effect sizes. If no omega is provided, omega_sequence is set to be seq(0.01, 1, by = 0.01)
#' @param r variable controlling dispersion of non-local priors. Default is 1. r must be >= 1
#'
#' @return Returns an S3 object of class `BFF` (see `BFF.object` for details).
#' @export
#'
#' @examples
#' regBFF = regression_test_BFF(t_stat = 1.5, n = 50, k = 3)
#' regBFF
#' plot(regBFF)
#'
regression_test_BFF <- function(
    t_stat,
    n = NULL,
    k = NULL,
    alternative = "two.sided",
    omega = NULL,
    omega_sequence = if(is.null(omega)) seq(0.01, 1, by = 0.01),
    r = 1){


  ### input checks and processing
  input <- .process_input.reg.test(t_stat, n, k, alternative, r)

  ### computation
  # calculate BF
  results   <- backend_reg(
    input     = input,
    r         = r,
    omega     = if(!is.null(omega)) omega else omega_sequence
  )

  ## compute minimum BFF for anything larger than small effect sizes
  if (is.null(omega)) {
    minimums = get_min_omega_bff(omega = omega_sequence, bff = results, cutoff = 0.02)
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
    test_type    = "regression_test",
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


.process_input.reg.test <- function(t_stat, n, k, alternative, r){

  if (r < 1)
    stop("r must be greater than or equal to 1")

  .check_alternative(alternative)

  df <- n -k - 1

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
    df         = df,
    k = k,
    alternative          = alternative,
    alternative.original = alternative.original
  ))
}
