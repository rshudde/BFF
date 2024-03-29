% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FINAL_FUNCTIONS_t_test.R
\name{t_test_BFF}
\alias{t_test_BFF}
\title{t_test_BFF}
\usage{
t_test_BFF(
  t_stat,
  n = NULL,
  one_sample = FALSE,
  alternative = "two.sided",
  n1 = NULL,
  n2 = NULL,
  r = NULL,
  omega = NULL
)
}
\arguments{
\item{t_stat}{T statistic}

\item{n}{sample size (if one sample test)}

\item{one_sample}{is test one sided? Default is FALSE}

\item{n1}{sample size of group one for two sample test.}

\item{n2}{sample size of group two for two sample test}

\item{r}{r value}

\item{tau2}{tau2 values (can be a single entry or a vector of values)}
}
\value{
Returns Bayes factor function results
 \tabular{ll}{
   \code{BFF} \tab The log of the Bayes Factor Function values \cr
   \tab \cr
   \code{effect_size} \tab Effect sizes tested (seq(0, 1, by = 0.01)) \cr
   \tab \cr
   \code{BFF_max_RMSE} \tab Maximum BFF value \cr
   \tab \cr
   \code{max_RMSE} \tab Effect size that maximizes BFF\cr
   \tab \cr
   \code{omega} \tab omega values tested, can be a single number or vector\cr
}
}
\description{
t_test_BFF constructs BFFs based on the t test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
All results are on the log scale.
}
\examples{
tBFF = t_test_BFF(t_stat = 2.5, n = 50)
t_test_BFF(t_stat = 2.5, n = 50, omega = 0.5)
t_test_BFF(t_stat = 2.5, n = 50, omega = c(0.5, 0.2))
t_test_BFF(t_stat = 2.5, n1 = 50, n2 = 40, one_sample = FALSE)
t_test_BFF(t_stat = 2.5, n = 50, r = 2)
t_test_BFF(t_stat = 2.5, r = 2, n1 = 50, n2 = 30, one_sample = FALSE)
t_test_BFF(t_stat = 2.5, n = 50, r = 2.5)
t_test_BFF(t_stat=2.5, r = 2.5, n1 = 50, n2 = 30,  one_sample = FALSE)
t_test_BFF(t_stat = 2.5, n = 50)
t_test_BFF(t_stat = 2.5, n = 50, omega = 0.5)
t_test_BFF(t_stat = 2.5, n = 50, tau2 = c(0.5, 0.8))
tBFF$BFF_max_RMSE   # maximum BFF omega
tBFF$max_RMSE       # effect size which maximizes the BFF value

}
