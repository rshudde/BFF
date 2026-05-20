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
  scale <- .t_test_ncp_scale(n = n, n1 = n1, n2 = n2, one_sample = one_sample)
  lambda <- scale * effect_size
  density <- dnlnm(x = lambda, tau2 = tau2, r = r) * scale

  if(one_sided){
    density <- ifelse(effect_size >= 0, 2 * density, 0)
  }

  return(density)
}

.t_test.posterior <- function(t_stat, tau2, r, effect_size, n = NULL, n1 = NULL, n2 = NULL, one_sample = FALSE, one_sided = FALSE){
  df <- .t_test_df(n = n, n1 = n1, n2 = n2, one_sample = one_sample)
  scale <- .t_test_ncp_scale(n = n, n1 = n1, n2 = n2, one_sample = one_sample)

  lik_prior <- .t_test.prior(
    tau2       = tau2,
    r          = r,
    effect_size = effect_size,
    n          = n,
    n1         = n1,
    n2         = n2,
    one_sample = one_sample,
    one_sided  = one_sided
  )
  lik_t <- suppressWarnings(stats::dt(x = t_stat, df = df, ncp = scale * effect_size))
  m1 <- .m1.t_test(t = t_stat, tau2 = tau2, r = r, df = df, two_sided = !one_sided)

  post_lik <- (lik_t * lik_prior) / m1
  return(post_lik)
}

.t_test_df <- function(n = NULL, n1 = NULL, n2 = NULL, one_sample = FALSE){
  if(one_sample){
    return(n - 1)
  }
  return(n1 + n2 - 2)
}

.t_test_ncp_scale <- function(n = NULL, n1 = NULL, n2 = NULL, one_sample = FALSE){
  if(one_sample){
    return(sqrt(n))
  }
  return(sqrt(n1 * n2 / (n1 + n2)))
}

.m1.t_test <- function(t, tau2, r, df, two_sided){
  stats::dt(x = t, df = df, ncp = 0) *
    exp(BFF_t_test(tau2 = tau2, t_stat = t, r = r, two_sided = two_sided, df = df))
}

### helper functions for one-sided t-test (i.e., t_test1)
# likelihood
.dt.t_test1 <- function(t, df, n, effect_size){
  suppressWarnings(stats::dt(x = t, df = df, ncp = .t_test_ncp_scale(n = n, one_sample = TRUE) * effect_size))
}

# for prior
.dnlnm.t_test1 <- function(n, effect_size, tau2, r) {
  density = ifelse(effect_size >= 0, dnlnm(x = .t_test_ncp_scale(n = n, one_sample = TRUE) * effect_size, tau2 = tau2, r = r)*2, 0)
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
  suppressWarnings(stats::dt(x = t, df = df, ncp = .t_test_ncp_scale(n1 = n1, n2 = n2, one_sample = FALSE) * effect_size))
}

# marginal under null
.m0.t_test2 <- function(t, df, lambda){
  m0_2s_t = stats::dt(x = t, df = df, ncp = 0)
  return(m0_2s_t)
}

# for prior
.dnlnm_t_test2 <- function(n1, n2, effect_size, tau2, r){
  density = dnlnm(.t_test_ncp_scale(n1 = n1, n2 = n2, one_sample = FALSE) * effect_size, tau2, r)
  return(density)
}

# I1 in closed form expression of marginal under alternative (in supplemental material)
.I1.t_test2 <- function(tau2, r, df, t){
  c = 1/((1+tau2)^(r + 0.5))
  gauss = Gauss2F1(a = (df +1)/2, b = (r + 0.5), c = 0.5, x = (tau2*(t^2))/(((t^2)+df)*(1 + tau2)))
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


### z_test
.z_test.prior <- function(tau2, r, effect_size, n = NULL, n1 = NULL, n2 = NULL, one_sample = FALSE, one_sided = FALSE){
  scale <- .t_test_ncp_scale(n = n, n1 = n1, n2 = n2, one_sample = one_sample)
  lambda <- scale * effect_size
  density <- dnlnm(x = lambda, tau2 = tau2, r = r) * scale

  if(one_sided){
    density <- ifelse(effect_size >= 0, 2 * density, 0)
  }

  return(density)
}

.z_test.posterior <- function(z_stat, tau2, r, effect_size, n = NULL, n1 = NULL, n2 = NULL, one_sample = FALSE, one_sided = FALSE){
  scale <- .t_test_ncp_scale(n = n, n1 = n1, n2 = n2, one_sample = one_sample)

  lik_prior <- .z_test.prior(
    tau2       = tau2,
    r          = r,
    effect_size = effect_size,
    n          = n,
    n1         = n1,
    n2         = n2,
    one_sample = one_sample,
    one_sided  = one_sided
  )
  lik_z <- stats::dnorm(x = z_stat, mean = scale * effect_size, sd = 1)
  m1 <- .m1.z_test(z = z_stat, tau2 = tau2, r = r, two_sided = !one_sided)

  post_lik <- (lik_z * lik_prior) / m1
  return(post_lik)
}

.m1.z_test <- function(z, tau2, r, two_sided){
  stats::dnorm(x = z, mean = 0, sd = 1) *
    exp(BFF_z_test(tau2 = tau2, z_stat = z, r = r, two_sided = two_sided))
}


### regression_test
.regression_test.prior <- function(tau2, r, effect_size, n, k, one_sided = FALSE){
  scale <- .regression_test_ncp_scale(n = n, k = k)
  lambda <- scale * effect_size
  density <- dnlnm(x = lambda, tau2 = tau2, r = r) * scale

  if(one_sided){
    density <- ifelse(effect_size >= 0, 2 * density, 0)
  }

  return(density)
}

.regression_test.posterior <- function(t_stat, tau2, r, effect_size, n, k, one_sided = FALSE){
  df <- .regression_test_df(n = n, k = k)
  scale <- .regression_test_ncp_scale(n = n, k = k)

  lik_prior <- .regression_test.prior(
    tau2       = tau2,
    r          = r,
    effect_size = effect_size,
    n          = n,
    k          = k,
    one_sided  = one_sided
  )
  lik_t <- suppressWarnings(stats::dt(x = t_stat, df = df, ncp = scale * effect_size))
  m1 <- .m1.regression_test(t = t_stat, tau2 = tau2, r = r, df = df, two_sided = !one_sided)

  post_lik <- (lik_t * lik_prior) / m1
  return(post_lik)
}

.regression_test_df <- function(n, k){
  n - k - 1
}

.regression_test_ncp_scale <- function(n, k){
  sqrt(.regression_test_df(n = n, k = k))
}

.m1.regression_test <- function(t, tau2, r, df, two_sided){
  stats::dt(x = t, df = df, ncp = 0) *
    exp(BFF_reg_test(tau2 = tau2, t_stat = t, df = df, r = r, two_sided = two_sided))
}


### chi2_test
.chi2_test.prior <- function(tau2, r, effect_size, n, df){
  lambda <- .chi2_test_ncp(effect_size = effect_size, n = n, df = df)
  jacobian <- .chi2_test_ncp_jacobian(effect_size = effect_size, n = n, df = df)

  density <- stats::dgamma(
    x     = lambda,
    shape = df/2 + r,
    rate  = 1/(2*tau2)
  ) * jacobian
  density <- ifelse(effect_size >= 0, density, 0)

  return(density)
}

.chi2_test.posterior <- function(chi2_stat, tau2, r, effect_size, n, df){
  lambda <- .chi2_test_ncp(effect_size = effect_size, n = n, df = df)
  lik_prior <- .chi2_test.prior(
    tau2       = tau2,
    r          = r,
    effect_size = effect_size,
    n          = n,
    df         = df
  )
  lik_chi2 <- suppressWarnings(stats::dchisq(x = chi2_stat, df = df, ncp = lambda))
  m1 <- .m1.chi2_test(chi2_stat = chi2_stat, tau2 = tau2, r = r, df = df)

  post_lik <- (lik_chi2 * lik_prior) / m1
  return(post_lik)
}

.chi2_test_ncp <- function(effect_size, n, df){
  n * df * effect_size^2
}

.chi2_test_ncp_jacobian <- function(effect_size, n, df){
  2 * n * df * effect_size
}

.m1.chi2_test <- function(chi2_stat, tau2, r, df){
  stats::dchisq(x = chi2_stat, df = df, ncp = 0) *
    exp(BFF_chi2_test(tau2 = tau2, chi2_stat = chi2_stat, k = df, r = r))
}


### f_test
.f_test.prior <- function(tau2, r, effect_size, n, df1){
  lambda <- .f_test_ncp(effect_size = effect_size, n = n, df1 = df1)
  jacobian <- .f_test_ncp_jacobian(effect_size = effect_size, n = n, df1 = df1)

  density <- stats::dgamma(
    x     = lambda,
    shape = df1/2 + r,
    rate  = 1/(2*tau2)
  ) * jacobian
  density <- ifelse(effect_size >= 0, density, 0)

  return(density)
}

.f_test.posterior <- function(f_stat, tau2, r, effect_size, n, df1, df2){
  lambda <- .f_test_ncp(effect_size = effect_size, n = n, df1 = df1)
  lik_prior <- .f_test.prior(
    tau2       = tau2,
    r          = r,
    effect_size = effect_size,
    n          = n,
    df1        = df1
  )
  lik_f <- suppressWarnings(stats::df(x = f_stat, df1 = df1, df2 = df2, ncp = lambda))
  m1 <- .m1.f_test(f_stat = f_stat, tau2 = tau2, r = r, df1 = df1, df2 = df2)

  post_lik <- (lik_f * lik_prior) / m1
  return(post_lik)
}

.f_test_ncp <- function(effect_size, n, df1){
  n * df1 * effect_size^2 / 2
}

.f_test_ncp_jacobian <- function(effect_size, n, df1){
  n * df1 * effect_size
}

.m1.f_test <- function(f_stat, tau2, r, df1, df2){
  stats::df(x = f_stat, df1 = df1, df2 = df2, ncp = 0) *
    exp(BFF_f_test(tau2 = tau2, f_stat = f_stat, k = df1, m = df2, r = r))
}

