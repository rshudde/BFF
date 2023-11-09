### generic function non-local normal moment distribution ----
# It would be nice if we export this to the user as they can
# use it for many other different things

#' @title Non-local Normal Moment Distribution
#'
#' @param x vector of quantiles.
#' @param tau2 the tau2 parameter
#' @param r the r parameter
#' @param log logical; if \code{TRUE}, probabilities
#' \code{p} are given as \code{log(p)}.
#'
#' @return \code{dnlnm} gives the density of non-local
#' normal moment distribution
#'
#' @name nlnm
#' @export dnlnm
#' @examples
NULL

# this allows us to add other functions for the nlnm distribution
# in accordance with the common R naming scheme:
# r(nlnm) = for generating random numbers
# q(nlnm) = for quantiles
# p(nlnm) = for cdf

#' @rdname nlnm
dnlnm <- function(x, tau2, r, log = FALSE){

  # TODO: add input checks

  lik <- (abs(x))^(2*r) * exp(-x^2/(2*tau2))/((2*tau2)^(r+0.5)*gamma(r+0.5))

  if(log)
    return(log(lik))
  else
    return(lik)
}


### test-specific prior and posterior distributions ----
# do not need to be exported as we use them for plotting

# let's use the following naming scheme
# (test type).(prior/posterior)
#
# test types:
# - z_test
# - t_test
# - chi2_test
# - f_test
#
# one/two sample as well as one/two sided version dispatched internally via an argument



### t_test
.t_test.prior     <- function(x, tau2, r, one_sample = FALSE, one_sided = FALSE){

  if(one_sample){
    stop("TODO")
  }else{
    lik <- dnlnm(x, tau2, r)
  }

  # apply side restriction
  if(one_sided){
    lik <- ifelse(x >= 0, 2*lik, 0)
  }

  return(lik)
}
.t_test.posterior <- function(x, tau2, r, t_stat, df, one_sample = FALSE, one_sided = FALSE){

  lik_prior <- .t_test.prior(x = x, tau2 = tau2, r = r, one_sample = one_sample, one_sided = one_sided)

  if(one_sample){
    if(one_sided){
      stop("TODO")
    }else{
      stop("TODO")
    }
  }else{

    lambda <- 2 # TODO???

    I1 <- .I1.t_test2(tau2 = tau2, r = r, df = df, t = x)
    I2 <- .I2.t_test2(tau2 = tau2, r = r, df = df, t = x)
    m0 <- .m0.t_test2(t = x, df = df, lambda = lambda)

    if(one_sided){
      m1  <- 2*m0*(I1 + I2)
      lik <- (ifelse(x > 0, 2*stats::dt(x = x, df = df, ncp = lambda), 0) * lik_prior)/m1
    }else{
      m1  <- m0*(I1 + I2)
      lik <- (stats::dt(x = x, df = df, ncp = lambda) * lik_prior)/m1
    }
  }

  return(lik)
}

### helper functions for two-sample t-test (i.e., t_test2)
# marginal under null
.m0.t_test2 <- function(t, df, lambda){
  m0_1s_t <- stats::dt(t, df = df, ncp = lambda)
  return(m0_1s_t)
}
# I1 in closed form expression of marginal under alternative (in supplemental material)
.I1.t_test2 <- function(tau2, r, df, t){
  c     <- 1/(2*((1+tau2)^(r+0.5)))
  gauss <- Gauss2F1(a = (df+1)/2, b = (r + 0.5), c = 0.5, x = (tau2*(t^2))/((t^2 + df)*(tau2 + 1)))
  I1    <-  c*gauss
  return(I1)
}
# I2 in closed form expression of marginal under alternative (in supplemental material)
.I2.t_test2 <- function(tau2, r, df, t){
  c     <- ((t*sqrt(tau2))*(gamma(df/2 + 1))*gamma(r+1))/((sqrt(t^2 + df))*((tau2+1)^(r+1))*gamma((df+1)/2)*gamma(r + 0.5))
  gauss <- Gauss2F1(a = ((df/2)+1), b = (r+1), c = 1.5, x = (((t^2)*tau2)/((t^2 + df)*(1+tau2))))
  I2    <- c*gauss
  return(I2)
}
