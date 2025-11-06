# TODO ask val about tau2 for one sample test
# TODO ask val about warning "In dt(t, df = df, ncp = lambda, log = TRUE) : full precision may not have been achieved in pnt{final}"

# x <- function(i){
#   if (i < 10) warning("A warning")
#   return(i)
# }
#
# a = tryCatch(x(5),warning=function(w) return(list(x(5),w)))
# b = tryCatch(x(15),warning=function(w) return(list(x(5),w)))


integrand_invm = function(lambda,t,tau2,nu,df, default_max){
  # set initial values of nu an lambda^2 based on input and control for overflow
  nu_half = 0.5*nu
  lambda2 = lambda^2
  lambda2[lambda2<(10^(-25))] = 10^(-25) # control for overflow

  # do the actual integration
  # is_warning = FALSE
  # dt_term = tryCatch(dt(t,df=df,ncp=lambda,log=TRUE),warning=function(w) return(list(x(5),w)))
  # if (length(dt_term) == 2) {
  #   dt_term = dt_term[[1]]
  #   is_warning = TRUE
  # }

  arg = -tau2/lambda^2-(0.5*(nu+1))*log(lambda2)+nu_half*log(tau2)- lgamma(nu_half) + stats::dt(t,df=df,ncp=lambda,log=TRUE)

  arg[arg<(-default_max)]=(-default_max) # cut off anything that is too small
  x = exp(arg)
  return(x)
}

# backend_t_invm = function(t,n1,n2,nu,omega, default_max = 700){ # Two-sided t test with IM(nu,tau(omega))#  prior}

BFF_t_test_invm = function(tau2, t_stat, nu, df, default_max) {
  BFF = tryCatch(stats::integrate(integrand_invm,
                                lower=-Inf,
                                upper=Inf,
                                t=t_stat,
                                tau2=tau2,
                                nu=nu,
                                df=df,
                                default_max=default_max,
                                rel.tol=.Machine$double.eps^.125),
                 warning = function(w)
                   # return(list(5, w)))
                 return(list(suppressWarnings(stats::integrate(integrand_invm,
                                       lower=-Inf,
                                       upper=Inf,
                                       t=t_stat,
                                       tau2=tau2,
                                       nu=nu,
                                       df=df,
                                       default_max=default_max,
                                       rel.tol=.Machine$double.eps^.125))$value,w)))
  is_warning = FALSE
  if (length(BFF) == 2) {
    BFF = BFF[[1]]
    is_warning = TRUE
    # print("The non-central t may have precision input issues stemming from the dt() function. The estimate is still returned")
  }
  # calculate logs and control for overflow
  # first, catch if any errors in dt
  dt_term = stats::dt(t_stat,df,ncp=0,log=TRUE)

  log_BF = log(BFF)-dt_term
  to_return = min(log_BF,default_max)

  return(list("to_return" = to_return, "is_warning" = is_warning))
}

backend_t_invm <- function(
    input,
    omega = NULL){

    # compute tau2 from omega
    # if multiple omegas and t-stats are supplied, each element of tau2
    # corresponds a vector of tau2 for the corresponding t-statistics
    # i.e., tau2[omega][t-stat]
    # TODO assuming same tau2 for one and two sample. ask Val
    tau2 <- lapply(omega, function(x){
      if(input$one_sample){
        tau2 <- get_one_sample_invm_t_tau2(n = input$n, w = x, nu = input$nu)
      }else{
        tau2 <- get_two_sample_invm_t_tau2(n1 = input$n1, n2 = input$n2, w = x, nu = input$nu)
      }
    })

    # # compute log_BF
    # log_BF <- sapply(tau2, function(x){
    #   sum(sapply(seq_along(input$t_stat), function(i){
    #     BFF_t_test_invm(
    #       input = input,
    #       tau2=x[i]
    #     )
    #   }))
    # })

    # compute log_BF
    log_BF <- sapply(tau2, function(x){
      sum(sapply(seq_along(input$t_stat), function(i){
        BFF_t_test_invm(
          tau2 = x[i],
          t_stat    = input$t_stat[i],
          nu = input$nu,
          df = input$df,
          default_max = input$default_max
        )$to_return
      }))
    })


  # recompute to check for any warnings
  # TODO have a more elegant solution for this, this is not great yet
  count = 1
  warnings_list = vector()
  for (i in 1:length(tau2)) {
    for (j in 1:length(input$t_stat)) {
      warnings_list[count] = BFF_t_test_invm(
        tau2 = tau2[[i]],
        t_stat    = input$t_stat[j],
        nu = input$nu,
        df = input$df,
        default_max = input$default_max
      )$is_warning

      count = count + 1
    }
  }

  if (any(warnings_list)) {
    print("The non-central t has potential precision issues stemming from dt(). The estimate is still returned")
  }

  return(log_BF)
}


################# t function user interaction

#' t_test_BFF_invm
#'
#' t_test_BFF constructs BFFs based on the t test using the inverse moment prior. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
#' By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale.
#'
#' @param t_stat t statistic
#' @param n sample size (if one sample test)
#' @param nu hyperparemeter for the inverse moment prior
#' @param n1 sample size of group one for two sample test. Must be provided if one_sample = FALSE
#' @param n2 sample size of group two for two sample test. Must be provided if one_sample = FALSE
#' @param one_sample is test one sided? Default is FALSE
#' @param alternative the alternative. options are "two.sided" or "less" or "greater"
#' @param omega standardized effect size. For the t-test, this is often called Cohen's d (can be a single entry or a vector of values)
#' @param omega_sequence sequence of standardized effect sizes. If no omega is provided, omega_sequence is set to be seq(0.01, 1, by = 0.01)
#' @param default_max set to 700, required for integration to correctly calculate. Authors do not suggest changing this value, as it is for computation only
#'
#' @return Returns an S3 object of class `BFF` (see `BFF.object` for details).
#' @export
#'
#' @examples
#' tBFF = t_test_BFF_invm(t_stat = 0.5, n = 50, nu = 1, one_sample = TRUE)
#' tBFF
#' plot(tBFF)

t_test_BFF_invm <- function(
    t_stat,
    n = NULL,
    nu = 1,
    n1 = NULL,
    n2 = NULL,
    one_sample = FALSE,
    alternative = "two.sided",
    omega = NULL,
    omega_sequence = if(is.null(omega)) seq(0.01, 1, by = 0.01),
    default_max = 700
    )

{
  ### input checks and processing
  input <- .process_input.t.test.invm(t_stat, n, n1, n2, nu, one_sample, alternative, default_max)

  ### computation
  # calculate BF
  results   <- backend_t_invm(
    input     = input,
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
    test_type    = "t_test",
    generic_test = FALSE,
    r            = NULL,
    input        = input
  )
  if(is.null(omega)){
    output$BFF = list(log_bf = log_bf, omega = omega_sequence)
  }

  class(output) = "BFF"
  return(output)
}


.process_input.t.test.invm <- function(t_stat, n, n1, n2, nu, one_sample, alternative, default_max){


  # if (r < 1)
  #   stop("r must be greater than or equal to 1")
  #
  # if (nu < 0)
  #   stop("nu must be greater than or equal to 0")

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
    nu         = nu,
    df         = df,
    one_sample = one_sample,
    alternative          = alternative,
    alternative.original = alternative.original,
    default_max = default_max
  ))
}


# Code from Val
# InvmomentBF = function(t,n1,n2,nu,omega){ # Two-sided t test with IM(nu,tau(omega))
#
#   #  prior}
#
#   M = length(t)
#
#   bf = c(rep(0,M))
#
#   tau = n1*n2*omega^2*(nu+1)/(n1+n2) # M vector
#
#   df = n1+n2-2
#
#
#
#   for(m in 1:M ){
#
#     x = integrate(integrand,lower=-Inf,upper=Inf,t=t[m],
#
#                   tau=tau[m],nu=nu,df=df[m],rel.tol=.Machine$double.eps^.125)$value
#
#     arg = log(x)-dt(t[m],df[m],ncp=0,log=TRUE)
#
#     arg = min(arg,700)
#
#     bf[m] = exp(arg)
#
#   }
#
#   return(bf)
#
# }
#
#
#
# integrand = function(lambda,t,tau,nu,df){
#
#   nuhalf = 0.5*nu
#
#   lambda2 = lambda^2
#
#   lambda2[lambda2<(10^(-25))] = 10^(-25)
#
#   arg = -tau/lambda^2-(0.5*(nu+1))*log(lambda2)+nuhalf*log(tau)-
#
#     lgamma(nuhalf)+dt(t,df=df,ncp=lambda,log=TRUE)
#
#   arg[arg<(-700)]=(-700)
#
#   x = exp(arg)
#
#   return(x)
#
# }

