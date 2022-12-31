# functions to set tau2 - user does not interact with these
get_one_sample_tau2 = function(n, w)
{
  to_return = n * w^2 / 2
  return(to_return)
}

get_two_sample_tau2 = function(n1, n2, w)
{
  to_return = n1*n2*w^2 / (n1 + n2)
  return(to_return)
}

get_count_tau2 = function(n, w)
{
  to_return = n*w^2
  return(to_return)
}

get_LRT_tau2 = function(n, w)
{
  to_return = n*w^2
  return(to_return)
}

get_linear_tau2 = function(n, w)
{
  to_return = n*w^2 / 2
  return(to_return)
}

plot_BFF = function(effect_size, BFF, savename, show)
{
  if (max(BFF) > 2*log(5000)) warning("Very strog evidence for the alternative, plots may not render correctly.")
  if(show) pdf(savename)
  plot(effect_size, BFF,
       type='l',
       xlab= expression(paste("RMSE ", tilde(omega))),
       ylab = "Bayes Factor Against Null Hypothesis",
       main = "BFF",
       lwd=3,cex=2,yaxt="n")
  if (min(BFF) < 0) abline(h = log(1))
  rect(-1, -log(10000), 0.1, log(1e10), col=adjustcolor("red", 0.1))
  rect(0.1, -log(10000), .35, log(1e10), col=adjustcolor("orange", 0.1))
  rect(.35, -log(10000), .65, log(1e10), col=adjustcolor("blue", 0.1))
  rect(.65, -log(10000), 2, log(1e10), col=adjustcolor("green", 0.1))
  axis(2, at=c(log(5000),log(2000),log(1000),log(500),log(200),log(100),log(50),log(20),log(10),log(5),log(2),log(1),
               -log(2),-log(5),-log(10),-log(20),-log(50),-log(100),-log(200),-log(500),-log(1000),-log(2000),-log(5000)),
       labels = c("5000:1","2000:1","1000:1","500:1","200:1","100:1","50:1","20:1","10:1","5:1","2:1","1:1",
                 "1:2","1:5","1:10","1:20","1:50","1:100","1:200","1:500","1:1000","1:2000","1:5000"))
  if(show) dev.off()
}
#' BFF_z_test
#'
#' Bayes Factor function test for the t statistic. Computes the Bayes factor in favor of the
#' alternative given a z statistic and sample size. The plot shown
#' when running the function is saved to "BFF_plot.pdf."
#'
#' @param z_stat t statistic
#' @param n Sample size
#' @param one_sample Is this a one or two sample z-test? Default is FALSE
#' @param n1 Sample size of group 1 if one_sample is FALSE
#' @param n2 Sample size of group 2 if one_sample is FALSE
#' @param savename Name of pdf file to save. Requires .pdf extension. Required if saving plot
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
#' }
#' @export
#'
#' @examples
#' BFF_z_test(1.4, 500)
BFF_z_test = function(z_stat, n = NULL, one_sample = TRUE, n1 = NULL, n2 = NULL, savename = NULL)
{
  if (is.null(n) & (is.null(n1) & is.null(n2))) stop("Either n or n1 and n2 is required")
  # get correct tau^2 value
  effect_size = seq(0, 1, by = 0.01)
  # if (is.null(tau2))
  # {
  if (one_sample)
  {
    tau2 = get_one_sample_tau2(n=n, w=effect_size)
  } else {
    if (is.null(n)) n = n1
    tau2 = get_two_sample_tau2(n1=n, n2=n2, w=effect_size)
  }
  # }

  term_one = (tau2 + 1)^(-3/2)
  term_two = 1 + tau2 * z_stat^2 / (tau2 + 1)
  term_three = exp(tau2 * z_stat^2/ (2 * (tau2 + 1)))

  to_return = term_one * term_two * term_three

  # stuff to return
  BFF = to_return
  effect_size = effect_size
  idx_max = which.max(BFF)
  BFF_max_RMSE = BFF[idx_max]
  max_RMSE = effect_size[idx_max]

  # plotting
  if(!is.null(savename)) plot_BFF(effect_size = effect_size, BFF = BFF, savename = savename, show = FALSE)
  plot_BFF(effect_size = effect_size, BFF = BFF, savename = savename, show = TRUE)

  return(list(BFF=BFF,
              effect_size=effect_size,
              BFF_max_RMSE=BFF_max_RMSE,
              max_RMSE=max_RMSE,
              tau2 = tau2))
}

#' BFF_t_test
#'
#' Bayes Factor function test for the t statistic. Computes the Bayes factor in favor of the
#' alternative given a chi^2 statistic, the degrees of freedom, and sample size. The plot shown
#' when running the function is saved to "BFF_plot.pdf."
#'
#' @param t_stat t statistic
#' @param df Degrees of freedom
#' @param n Sample size
#' @param one_sample Is this a one or two sample z-test? Default is FALSE
#' @param n1 Sample size of group 1 if one_sample is FALSE
#' @param n2 Sample size of group 2 if one_sample is FALSE
#' @param savename Name of pdf file to save. Requires .pdf extension. Required if saving plot
#'
#' @return Returns Bayes factor function results
#'  \tabular{ll}{
#'    \code{BFF} \tab Bayes Factor Function values \cr
#'    \tab \cr
#'    \code{effect_size} \tab Effect sizes tested (seq(0, 1, by = 0.01)) \cr
#'    \tab \cr
#'    \code{max_BFF} \tab Maximum BFF value \cr
#'    \tab \cr
#'    \code{max_RMSE} \tab Effect size that maximizes BFF\cr
#' }
#' @export
#'
#' @examples
#' BFF_t_test(1.4, 10, n = 100)
BFF_t_test = function(t_stat, df, n = NULL, one_sample = TRUE, n1 = NULL, n2 = NULL, savename = NULL)
{
  if (is.null(n) & (is.null(n1) & is.null(n2))) stop("Either n or n1 and n2 is required")
  effect_size = seq(0, 1, by = 0.01)

  # get correct tau^2 value
  # if (is.null(tau2))
  # {
  if (one_sample)
  {
    tau2 = get_one_sample_tau2(n=n, w=effect_size)
  } else {
    if (is.null(n)) n = n1
    tau2 = get_two_sample_tau2(n1=n1, n2=n2, w=effect_size)
  }
  # }

  r = 1 + t_stat^2 / df
  s = 1 + t_stat^2 / (df*(1+ tau2))
  q = tau2 * (df + 1) / (df*(1+ tau2))

  BF = (tau2 + 1)^(-3/2) * (r/s)^((df+1)/2) * (1 + q*t_stat^2 / s)

  to_return = BF

  # stuff to return
  BFF = to_return
  effect_size = effect_size
  idx_max = which.max(BFF)
  BFF_max_RMSE = BFF[idx_max]
  max_RMSE = effect_size[idx_max]

  # plotting
  if(!is.null(savename))  plot_BFF(effect_size = effect_size, BFF = BFF, savename = savename, show = FALSE)
  plot_BFF(effect_size = effect_size, BFF = BFF, savename = savename, show = TRUE)

  return(list(BFF=BFF,
              effect_size=effect_size,
              max_BFF=BFF_max_RMSE,
              max_RMSE=max_RMSE,
              tau2 = tau2))
}

#' BFF_chi2_test
#'
#' Bayes Factor function test for the chi^2 statistic. Computes the Bayes factor in favor of the
#' alternative given a chi^2 statistic, the degrees of freedom, and sample size. The plot shown
#' when running the function is saved to "BFF_plot.pdf."
#'
#' @param chi_stat chi^2 statistic
#' @param df Degrees of freedom
#' @param n Sample size
#' @param count Is this a test of Pearson’s chi^2 test for goodness-of-fit? Default is TRUE. FALSE assumes a likelihiood ratio test
#' @param savename Name of pdf file to save. Requires .pdf extension. Required if saving plot
#'
#' @return Returns Bayes factor function results
#'  \tabular{ll}{
#'    \code{BFF} \tab Bayes Factor Function values \cr
#'    \tab \cr
#'    \code{effect_size} \tab Effect sizes tested (seq(0, 1, by = 0.01)) \cr
#'    \tab \cr
#'    \code{max_BFF} \tab Maximum BFF value \cr
#'    \tab \cr
#'    \code{max_RMSE} \tab Effect size that maximizes BFF\cr
#' }
#' @export
#'
#' @examples
#' BFF_chi2_test(chi_stat = 3.7, df = 10, n = 100)
BFF_chi2_test = function(chi_stat, df, n, count = TRUE, savename = NULL)
{
  effect_size = seq(0, 1, by = 0.01)

  # if (is.null(tau2))
  # {
  if (count)
  {
    tau2 = get_count_tau2(n=n, w=effect_size)
  } else {
    tau2 = get_LRT_tau2(n=n, w=effect_size)
  }
  # }

  BF = (tau2 + 1)^(-df/2 - 1) * (1 + tau2 * chi_stat / (df * (tau2 + 1))) * exp(tau2*chi_stat / (2*(tau2+1)))

  BFF = BF
  effect_size = effect_size
  idx_max = which.max(BFF)
  BFF_max_RMSE = BFF[idx_max]
  max_RMSE = effect_size[idx_max]

  # plotting
  if(!is.null(savename))  plot_BFF(effect_size = effect_size, BFF = BFF, savename = savename, show = FALSE)
  plot_BFF(effect_size = effect_size, BFF = BFF, savename = savename, show = TRUE)

  return(list(BFF=BFF,
              effect_size=effect_size,
              max_BFF=BFF_max_RMSE,
              max_RMSE=max_RMSE,
              tau2 = tau2))
}

#' BFF_F_test
#'
#' Bayes Factor function test for the F statistic. Computes the Bayes factor in favor of the
#' alternative given an F statistic, the degrees of freedom, and sample size. The plot shown
#' when running the function is saved to "BFF_plot.pdf."
#'
#' @param f_stat F statistic
#' @param df1 Degrees of freedom
#' @param df2 Degrees of freedom
#' @param n Sample size
#' @param savename Name of pdf file to save. Requires .pdf extension. Required if saving plot
#'
#' @return Returns Bayes factor function results
#'  \tabular{ll}{
#'    \code{BFF} \tab Bayes Factor Function values \cr
#'    \tab \cr
#'    \code{effect_size} \tab Effect sizes tested (seq(0, 1, by = 0.01)) \cr
#'    \tab \cr
#'    \code{max_BFF} \tab Maximum BFF value \cr
#'    \tab \cr
#'    \code{max_RMSE} \tab Effect size that maximizes BFF\cr
#' }
#' @export
#'
#' @examples
#' BFF_F_test(f_stat = 4.6, df1 = 4, df2 = 10, n = 100)
BFF_F_test = function(f_stat, df1, df2, n, savename = NULL)
{
  effect_size = seq(0, 1, by = 0.01)
  tau2 = get_linear_tau2(n, effect_size)

  v = df2 * (tau2 + 1)
  term_one = (tau2 + 1)^(-df1/2 - 1)
  term_two = (1 + df1 * f_stat / df2) / (1 + df1 * f_stat / v)
  term_three = 1 + (df1 + df2) * tau2 * f_stat / (v * (1+df1 * f_stat/v))

  to_return = term_one * term_two * term_three

  BFF = to_return
  effect_size = effect_size
  idx_max = which.max(BFF)
  BFF_max_RMSE = BFF[idx_max]
  max_RMSE = effect_size[idx_max]

  # plotting
  if(!is.null(savename))  plot_BFF(effect_size = effect_size, BFF = BFF, savename = savename, show = FALSE)
  plot_BFF(effect_size = effect_size, BFF = BFF, savename = savename, show = TRUE)

  return(list(BFF=BFF,
              effect_size=effect_size,
              max_BFF=BFF_max_RMSE,
              max_RMSE=max_RMSE,
              tau2 = tau2))
  return(to_return)
}
