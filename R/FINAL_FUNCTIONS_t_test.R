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

get_w = function(tau, v, t = t)
{
  num = 1 + t ^ 2 / (v * (1 + tau))
  den = 1 + t ^ 2 / v

  to_return = num / den
  return(to_return)
}

sum_val_t = function(r, m, tau, t, v, w)
{
  one = choose(2 * r, 2 * m)
  two = (2 * tau * t ^ 2 / ((t ^ 2 + v) * (tau + 1) * w)) ^ m
  three = sterling_gamma((v + 2 * m + 1) / 2) * double_factorial(2 * r - 2 *
                                                                   m - 1)

  to_return = one * two * three
  return(to_return)
}

sum_val_function_t = function(r, tau, t, v, w)
{
  val = 0
  for (mm in 0:r)
  {
    val = val + sum_val_t(
      r = r,
      m = mm,
      tau = tau,
      t = t,
      v = v,
      w = w
    )
  }
  return(val)
}

log_T = function(t, r, tau, v)
{
  w = get_w(tau = tau, v = v, t = t)
  num1 = 1
  den1 = double_factorial(2 * r - 1) * (1 + tau) ^ (r + 1 / 2) * sterling_gamma((v +
                                                                                   1) / 2) * w ^ ((v + 1) / 2)
  first_term = num1 / den1

  second_term = sum_val_function_t(
    r = r,
    tau = tau,
    t = t,
    v = v,
    w = w
  )

  to_return = first_term * second_term

  # log_version
  to_return = log(first_term) + log(second_term)

  return(to_return)
}

################# T functions if r is a fraction
log_T_frac = function(t, v, r, tau)
{
  tp1 = 1 + tau
  one = 1 / (tp1 ^ (r + 1 / 2))

  a1 = (v + 1) / 2
  b1 = r + 1 / 2
  c1 = 1 / 2
  d1 = tau * t ^ 2 / ((t ^ 2 + v) * tp1)
  two = Gauss2F1(a1, b1, c1, d1)

  three = t * sqrt(tau) / (sqrt(t ^ 2 + v) * tp1 ^ (r + 1))
  four = sterling_gamma(v / 2 + 1) / sterling_gamma((v + 1) / 2)
  five = sterling_gamma(r + 1) / sterling_gamma(r + 1 / 2)

  aa = v / 2 + 1
  bb = r + 1
  cc = 3 / 2
  dd = d1
  six = Gauss2F1(aa, bb, cc, dd)

  to_return = one * two + three * four * five * six
  to_return = log(to_return)

  return(to_return)
}

log_T_frac_onesided = function(t, v, tau, r)
{
  one = 1 / ((1 + tau) ^ (r + 1 / 2))

  a1 = (v + 1) / 2
  b1 = r + 1 / 2
  c1 = 1 / 2
  d1 = tau * t ^ 2 / ((t ^ 2 + v) * (tau + 1))
  two = Gauss2F1(a1, b1, c1, d1)

  three = 2 * t * sqrt(tau) / (sqrt(t ^ 2 + v) * (tau + 1) ^ (r + 1))
  four = sterling_gamma(v / 2 + 1) / sterling_gamma((v + 1) / 2)
  five = sterling_gamma(r + 1) / sterling_gamma(r + 1 / 2)

  aa = v / 2 + 1
  bb = r + 1
  cc = 3 / 2
  dd = d1
  six = Gauss2F1(aa, bb, cc, dd)

  to_return = one * two + three * four * five * six
  to_return = log(to_return)

  return(to_return)
}



################# T function user interaction

#' t.test.BFF
#'
#' T test using BFF methods. Setting r > 1 uses the higher order moments. Fractional moments are set with r > 1 and r not an integer.
#' All results are on the log scale. Plot saved to working directory unless a full path is specified in the 'savename' variable.
#'
#' @param t_stat T statistic
#' @param df degrees of freedom
#' @param n sample size (if one sample test)
#' @param one_sample is test one sided? Default is TRUE
#' @param n1 sample size of group one for two sample test
#' @param n2 sample size of group two for two sample test
#' @param savename optional, filename for saving the pdf of the final plot
#' @param r r value
#' @param save should a copy of the plot be saved?
#' @param xlab optional, x label for plot
#' @param ylab optional, y label for plot
#' @param main optional, main label for plot
#'
#' @return Returns Bayes factor function results
#'  \tabular{ll}{
#'    \code{BFF} \tab Bayes Factor Function values \cr
#'    \tab \cr
#'    \code{effect_size} \tab Effect sizes tested (seq(0, 1, by = 0.01)) \cr
#'    \tab \cr
#'    \code{BFF_max_RMSE} \tab Maximum BFF value \cr
#'    \tab \cr
#'    \code{max_RMSE} \tab Effect size that maximizes BFF\cr
#'    \tab \cr
#'    \code{tau2} \tab tau^2 values tested\cr
#' }
#' @export
#'
#' @examples
#' t.test.BFF(t_stat=2.5, n = 50, df = 20, save=FALSE)
#' t.test.BFF(t_stat=2.5, n1 = 50, n2 = 40, df = 20, save=FALSE, one_sample = FALSE)
t.test.BFF = function(t_stat,
                      n = NULL,
                      df = NULL,
                      one_sample = TRUE,
                      n1 = NULL,
                      n2 = NULL,
                      savename = NULL,
                      r = 1,
                      save = TRUE,
                      xlab = NULL,
                      ylab = NULL,
                      main = NULL)

{
  if (is.null(n) &
      (is.null(n1) &
       is.null(n2)))
    stop("Either n or n1 and n2 is required")

  if (is.null(df))
    stop("df is required")

  # same effect sizes for all tests
  effect_size = seq(0.01, 1, by = 0.01)

  r1 = FALSE
  frac_r = FALSE
  integer_r = FALSE

  if (r == 1) {
    r1 = TRUE
  } else if (r > 1)
  {
    remainder = floor(r) - r
    if (abs(remainder) < 1e-5) {
      integer_r = TRUE
    } else {
      frac_r = TRUE
    }
  }

  log_vals = rep(0, length(effect_size))
  if (r1) {
    if (one_sample)
    {
      tau2 = get_one_sample_tau2(n = n, w = effect_size)
    } else {
      tau2 = get_two_sample_tau2(n1 = n1, n2 = n2, w = effect_size)
    }
    log_vals = unlist(lapply(tau2, t_val_r1, t_stat=t_stat, df=df))
  }


  if (integer_r) {
    if (one_sample)
    {
      tau2 = get_tau_z_t_one_sample_frac(n = n, w = effect_size, r = r)
    } else {
      tau2 = get_tau_z_t_two_sample_frac(
        n1 = n,
        n2 = n2,
        w = effect_size,
        r = r
      )
    }
    log_vals = log_T(t = t_stat,
                     r = r,
                     tau = tau2,
                     v = df)
  }

  if (frac_r) {
    if (one_sample)
    {
      tau2 = get_tau_z_t_one_sample_frac(n = n, w = effect_size, r = r)
      log_vals = log_T_frac_onesided(t = t_stat, r = r, tau = tau2, v = df)
    } else {
      tau2 = get_tau_z_t_two_sample_frac(
        n1 = n,
        n2 = n2,
        w = effect_size,
        r = r
      )
      log_vals = log_T_frac_onesided(
        t = t_stat,
        r = r,
        tau = tau2,
        v = df
      )
    }
  }

  # stuff to return
  BFF = log_vals
  effect_size = effect_size
  idx_max = which.max(BFF)
  BFF_max_RMSE = BFF[idx_max]
  max_RMSE = effect_size[idx_max]

  # check the results are finite
  if (!all(is.finite(BFF)))
  {
    stop(
      "Values entered produced non-finite numbers. The most likely scenario is the evidence was so strongly in favor of the alternative that there was numeric overflow. Please contact the maintainer for more information."
    )
  }

  # plotting
  bff_plot = c()
  bff_plot[[1]] = BFF

  plot_BFF(
    effect_size = effect_size,
    BFF = bff_plot,
    save = save,
    savename = savename,
    xlab = xlab,
    ylab = ylab,
    main = main,
    r = r
  )

  return(
    list(
      BFF = BFF,
      effect_size = effect_size,
      BFF_max_RMSE = BFF_max_RMSE,
      max_RMSE = max_RMSE,
      tau2 = tau2
    )
  )

}