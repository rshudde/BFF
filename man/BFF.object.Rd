% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BFF-package.R
\name{BFF.object}
\alias{BFF.object}
\title{BFF Object}
\description{
The BFF object is a class that contains the results of a Bayes Factor Function.
The BFF object contains the following fields:
\tabular{ll}{
 \code{BFF} \tab the object containing the log_bf (log Bayes factor values) and corresponding omega values \cr
 \tab \cr
 \code{input} \tab the object containing the input values \cr
 \tab \cr
 \code{log_bf_h1} \tab maximized Bayes factor in favor of the alternative, or the corresponding Bayes factor if an omega is supplied\cr
 \tab \cr
 \code{omega_h1} \tab corresponding omega value for maximized Bayes factor in favor of the alternative, orthe corresponding omega for the Bayes factor if an omega is supplied\cr
 \tab \cr
 \code{log_bf_h0} \tab minimized Bayes factor in favor of the null (for medium or large effect sizes, only returned if a specific omega is not set) \cr
 \tab \cr
 \code{omega_h0} \tab corresponding omega value for minimized Bayes factor in favor of the null (for medium or large effect sizes, only returned if a specific omega is not set)\cr
 \tab \cr
 \code{alternative} \tab alternative hypothesis used in calculations \cr
 \tab \cr
 \code{omega_set} \tab was an omega value provided?\cr
 \tab \cr
 \code{r} \tab r value (default is 1 if not provided by user) \cr
 \tab \cr
 \code{test_type} \tab type of BFF test\cr
 \tab \cr
 \code{generic_test} \tab FALSE \cr
 \tab \cr
}
}
