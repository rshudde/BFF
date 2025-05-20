################# F functions if r is an integer and equal to 1 (PNAS)
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

####################### backend implementation (SPL, default)
BFF_f_test = function(tau2, f_stat, k, m, r)
{
  b = get_b(tau2 = tau2, r = r, k = k)

  numerator = k * tau2 * f_stat
  denomonator = (1 + tau2) * (m + k*f_stat)
  term_four = numerator / denomonator
  hypergeo = Gauss2F1(k/2 + r, (k+m)/2, k/2, term_four)

  final_BF =b * hypergeo
  to_return = log(final_BF)
  return(to_return)
}

####################### backend implementation
backend_f <- function(
    input,
    r,
    omega = NULL){

  # compute tau2 from omega
  # if multiple omegas and t-stats are supplied, each element of tau2
  # corresponds a vector of tau2 for the corresponding t-statistics
  # i.e., tau2[omega][t-stat]
  tau2 <- lapply(omega, function(x){
      tau2 = get_linear_tau2(n = input$n, w = x, k = input$df1, r = r)
  })

  # compute log_BF
  log_BF <- sapply(tau2, function(x){
    sum(sapply(seq_along(input$f_stat), function(i){
      BFF_f_test(
        tau2 = x[i],
        f_stat    = input$f_stat[i],
        k = input$df1,
        m = input$df2,
        r = r
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

#' f_test_BFF
#'
#' f_test_BFF constructs BFFs based on the t test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#'
#' @param f_stat T statistic
#' @param n sample size (if one sample test)
#' @param df1 sample size of group one for two sample test.
#' @param df2 sample size of group two for two sample test
#' @param omega standardized effect size. For the f-test, this is often called Cohen's f (can be a single entry or a vector of values)
#' @param omega_sequence sequence of standardized effect sizes. If no omega is provided, omega_sequence is set to be seq(0.01, 1, by = 0.01)
#' @param r variable controlling dispersion of non-local priors. Default is 1. r must be >= 1
#'
#' @return Returns an S3 object of class `BFF` (see `BFF.object` for details).
#' @export
#'
#' @examples
#' fBFF = f_test_BFF(f_stat = 1.5, n = 50, df1 = 25, df2 = 48)
#' fBFF
#' plot(fBFF)
#'
f_test_BFF = function(f_stat,
                      n,
                      df1,
                      df2,
                      omega = NULL,
                      omega_sequence = if(is.null(omega)) seq(0.01, 1, by = 0.01),
                      r = 1)


{

  ### input checks and processing
  input <- .process_input.f.test(f_stat, n, df1, df2, r)

  ### computation
  # calculate BF
  results   <- backend_f(
    input     = input,
    r         = r,
    omega     = if(!is.null(omega)) omega else omega_sequence
  )

  ## compute minimum BFF for anything larger than small effect sizes
  if (is.null(omega)) {
    minimums = get_min_omega_bff(omega = omega_sequence, bff = results, cutoff = 0.1)
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
    test_type    = "f_test",
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



.process_input.f.test <- function(f_stat, n, df1, df2, r) {

  if (r < 1)
    stop("r must be greater than or equal to 1")

  return(list(
    f_stat     = f_stat,
    n          = n,
    df1         = df1,
    df2         = df2
  ))
}


