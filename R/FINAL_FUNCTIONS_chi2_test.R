################# chih2 functions if r is an integer and equal to 1
G_val_r1 = function(tau2, chi2_stat, df)
{
  BFF = (tau2 + 1) ^ (-df / 2 - 1) * (1 + tau2 * chi2_stat / (df * (tau2 + 1))) * exp(tau2 *
                                                                                        chi2_stat / (2 * (tau2 + 1)))

  to_return = log(BFF)
  return(to_return)
}


####################### backend implementation (SPL, default)
BFF_chi2_test = function(tau2, chi2_stat, k, r)
{

  b = get_b(tau2=tau2, r = r, k = k)

  term_three = tau2 * chi2_stat / (2*(1 + tau2^2))
  hypergeo = hypergeom1F1(k/2 + r, k/2, term_three)$f

  final_BF = b*hypergeo
  to_return = log(final_BF)
  return(to_return)
}


####################### backend implementation
backend_chi2 <- function(
    input,
    r,
    omega = NULL){

  # compute tau2 from omega
  # if multiple omegas and t-stats are supplied, each element of tau2
  # corresponds a vector of tau2 for the corresponding t-statistics
  # i.e., tau2[omega][t-stat]
  tau2 <- lapply(omega, function(x){
    if(input$LRT){
      tau2 <- get_LRT_tau2(n = input$n, k = input$df, w = x, r = r)
    }else{
      tau2 <- get_count_tau2(n = input$n, k = input$df, w = x, r = r)
    }
  })

  # compute log_BF
  log_BF <- sapply(tau2, function(x){
    sum(sapply(seq_along(input$chi2_stat), function(i){
      BFF_chi2_test(
          tau2 = x[i],
          chi2_stat    = input$chi2_stat[i],
          k = input$df,
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
#' @param omega standardized effect size. For the chi^2-test, this is often called Cohen's w (can be a single entry or a vector of values)
#' @param omega_sequence sequence of standardized effect sizes. If no omega is provided, omega_sequence is set to be seq(0.01, 1, by = 0.01)
#' @param r variable controlling dispersion of non-local priors. Default is 1. r must be >= 1
#'
#' @return Returns an S3 object of class `BFF` (see `BFF.object` for details).
#' @export
#'
#' @examples
#' chi2BFF = chi2_test_BFF(chi2_stat = 6.5, n = 10, df = 9)
#' chi2BFF
#' plot(chi2BFF)
#'
chi2_test_BFF = function(chi2_stat,
                      n,
                      df,
                      LRT = FALSE,
                      omega = NULL,
                      omega_sequence = if(is.null(omega)) seq(0.01, 1, by = 0.01),
                      r = 1)

{
  ### input checks and processing
  input <- .process_input.chi2.test(chi2_stat, n, LRT, df, r)

  ### computation
  # calculate BF
  results   <- backend_chi2(
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
    test_type    = "chi2_test",
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




.process_input.chi2.test <- function(chi2_stat, n, LRT, df, r){

  if (r < 1)
    stop("r must be greater than or equal to 1")

  return(list(
    chi2_stat     = chi2_stat,
    n          = n,
    df         = df,
    LRT = LRT
  ))
}



