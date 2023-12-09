# BFF
Bayes Factor Functions

<!-- badges: start -->
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/BFF)](https://cran.r-project.org/package=BFF)[![](https://cranlogs.r-pkg.org/badges/BFF)](https://CRAN.R-project.org/package=BFF)
<!-- badges: end -->

This package provides the Bayes Factor values for different effect sizes from 0 to 1. A small effect size is usually considered from  0.2 to 0.5, medium effect sizes from 0.5 to 0.8, and large effect sizes as greater than 0.8.

Using this package is very similar to using the familiar t, z, chi^2, and F tests in R. You will need the same information - the test statistic, degrees of freedom, and sample size. A graph is produced that shows the BFF curve over the different effect sizes. 

For evaluating evidence from multiple studies (see 'Bayes factor functions', 2023 (arxiv)), the parameter 'r' can also be set. The default value for r is 1, but 'r' can be suggested that maximizes the bayes factor at each tau by setting the 'maximization' argument in each test to "TRUE." 


Installation
------------

The R package 'BFF' is available from CRAN, use the commands below to install the most recent Github version.

```{r, eval = FALSE}
# Plain installation
devtools::install_github("rshudde/BFF") # BFF package

```

Example
-------

```{r}
library(BFF)

z_BFF_one = z_test_BFF(z_stat = 2.5, n = 50) #one sample z-test
z_BFF_two = z_test_BFF(z_stat = 2.5, one_sample = FALSE, n1 = 50, n2 = 50) #two sample z-test

plot(z_BFF_two) #to view the plot of BFF vs the maximized omega (here for the two sample z-test)

```
