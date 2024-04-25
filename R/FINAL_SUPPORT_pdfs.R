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
NULL

# this allows us to add other functions for the nlnm distribution
# in accordance with the common R naming scheme:
# r(nlnm) = for generating random numbers
# q(nlnm) = for quantiles
# p(nlnm) = for cdf

#' @rdname nlnm
dnlnm <- function(x, tau2, r, log = FALSE){

  # TODO: add input checks

  lik <- ((abs(x))^(2*r) * exp((-x^2)/(2*tau2)))/(((2*tau2)^(r+0.5))*gamma(r+0.5))

  if(log){
    return(log(lik))
  }else{
    return(lik)
  }
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
.t_test.prior     <- function(tau2, r, effect_size, n = NULL, n1 = NULL, n2 = NULL, one_sample = FALSE, one_sided = FALSE){
  if(one_sided){
    lik = .dnlnm.t_test1(n = n, effect_size = effect_size, tau2 = tau2, r = r)*sqrt(n)
  }else{
    if(one_sample){
      stop("TODO")
    }else{
      lik = .dnlnm_t_test2(n1 = n1, n2 = n2, effect_size = effect_size, tau2 = tau2, r = r) * ((sqrt(2*n1*n2))/(sqrt(n1+n2)))
    }
  }
  return(lik)
}

.t_test.posterior <- function(t_stat, tau2, r, effect_size, n = NULL, n1 = NULL, n2 = NULL, one_sample = FALSE, one_sided = FALSE){

  if(one_sided){
    # TODO: please check whether the code for one-sample and two sample identical, if not, dispatch
    df <- n - 1
    lik_prior = .t_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n = n, n1 = n1, n2 = n2, one_sample = one_sample, one_sided = one_sided)
    m0 = .m0.t_test1(t = t_stat, df = df)
    I1 = .I1.t_test1(tau2 = tau2, r = r, df = df, t = t_stat)
    I2 = .I2.t_test1(tau2 = tau2, r = r, df = df, t = t_stat)
    m1 = 2*m0*(I1 + I2)
    post_lik = ((.dt.t_test1(t = t_stat, df = df, n = n, effect_size = effect_size) * lik_prior)/m1)
  }else{
    if(one_sample){
      # one-sample two-sided t-test
      # TODO: needs to be implemented
      stop("TODO")
    }else{
      # two-sample two-sided t-test
      df <- n1 + n2 - 2
      lik_prior <- .t_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n1 = n1, n2 = n2, one_sample = one_sample, one_sided = one_sided)
      m0 <- .m0.t_test2(t = t_stat, df = df)
      I1 <- .I1.t_test2(tau2 = tau2, r = r, df = df, t = t_stat)
      I2 <- .I2.t_test2(tau2 = tau2, r = r, df = df, t = t_stat)
      m1 <- m0*(I1 + I2)
      post_lik <- ((.dt.t_test2(t = t_stat, df = df, n1 = n1, n2 = n2, effect_size = effect_size) * lik_prior)/m1)
    }
}
return(post_lik)
}

### helper functions for one-sided t-test (i.e., t_test1)
# likelihood
.dt.t_test1 <- function(t, df, n, effect_size){
  ifelse(t>0, 2*stats::dt(x = t, df = df, ncp = sqrt(n)*effect_size), 0)
}

# for prior
.dnlnm.t_test1 <- function(n, effect_size, tau2, r) {
  density = ifelse(sqrt(n)*effect_size >= 0, dnlnm(x = sqrt(n)*effect_size, tau2 = tau2, r = r)*2, 0)
  return(density)
}

# marginal under null
.m0.t_test1 <- function(t, df){
  m0_1s_t = stats::dt(x = t, df = df, ncp = 0)
  return(m0_1s_t)
}

# I1 in closed form expression of marginal under alternative (in supplemental material)
.I1.t_test1 <- function(tau2, r, df, t){
  c = 1/(2*((1+tau2)^(r+0.5)))
  gauss = Gauss2F1(a = (df+1)/2, b = (r + 0.5), c = 0.5, x = (tau2*(t^2))/((t^2 + df)*(tau2 + 1)))
  I1 = c*gauss
  return(I1)
}

# I2 in closed form expression of marginal under alternative (in supplemental material)
.I2.t_test1 <-  function(tau2, r, df, t){
  c = ((t*sqrt(tau2))*(gamma((df/2) + 1))*gamma(r+1))/((sqrt(t^2 + df))*((tau2+1)^(r+1))*gamma((df+1)/2)*gamma(r + 0.5))
  gauss = Gauss2F1(a = ((df/2)+1), b = (r+1), c = 1.5, x = (((t^2)*tau2)/((t^2 + df)*(1+tau2))))
  I2 = c*gauss
  return(I2)
}

### helper functions for two-sample t-test (i.e., t_test2)
# likelihood
.dt.t_test2 <- function(t, df, n1, n2, effect_size){
  stats::dt(x = t, df = df, ncp =  ((sqrt(2*n1*n2))/(sqrt(n1+n2))) * effect_size)
}

# marginal under null
.m0.t_test2 <- function(t, df, lambda){
  m0_2s_t = stats::dt(x = t, df = df, ncp = 0)
  return(m0_2s_t)
}

# for prior
.dnlnm_t_test2 <- function(n1, n2, effect_size, tau2, r){
  density = dnlnm((((sqrt(2*n1*n2))/(sqrt(n1+n2))) * effect_size), tau2, r)
  return(density)
}

# I1 in closed form expression of marginal under alternative (in supplemental material)
.I1.t_test2 <- function(tau2, r, df, t){
  c = 1/((1+tau2)^(r + 0.5))
  gauss = Gauss2F1(a = (df +1)/2, b = (r = 0.5), c = 0.5, x = (tau2*(t^2))/(((t^2)+df)*(1 + tau2)))
  I1 = c*gauss
  return(I1)
}
# I2 in closed form expression of marginal under alternative (in supplemental material)
.I2.t_test2 <- function(tau2, r, df, t){
  c = (((sqrt(tau2))*t)*(gamma((df/2) + 1))*(gamma(r + 1)))/((sqrt((t^2) + df))*((1 + tau2)^(r + 1))*(gamma((df + 1)/2))*(gamma(r + 0.5)))
  gauss = Gauss2F1(a = ((df/2)+1), b = (r + 1), c = 1.5, x = (tau2*(t^2))/((1 + tau2)*((t^2) + df)))
  I2 = c*gauss
  return(I2)
}

