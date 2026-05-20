integrate_density <- function(f, lower, upper = Inf){
  suppressWarnings(stats::integrate(
    f = f,
    lower = lower,
    upper = upper,
    rel.tol = 1e-7,
    subdivisions = 1000
  )$value)
}

expect_posterior_plot_data <- function(plot_data){
  testthat::expect_true(is.data.frame(plot_data))
  testthat::expect_equal(colnames(plot_data), c("x", "prior", "posterior"))
  testthat::expect_true(all(is.finite(plot_data$prior)))
  testthat::expect_true(all(is.finite(plot_data$posterior)))
  testthat::expect_true(all(plot_data$prior >= 0))
  testthat::expect_true(all(plot_data$posterior >= 0))
}
