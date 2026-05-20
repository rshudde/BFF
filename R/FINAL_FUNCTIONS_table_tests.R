####################### table-test user interfaces

#' prop_test_BFF
#'
#' Constructs a BFF for a two-sample test of proportions from success counts
#' and sample sizes. The test statistic is the uncorrected Pearson chi-square
#' statistic, or the likelihood-ratio G-squared statistic when \code{LRT = TRUE}.
#' Effect-size transformations for odds ratios, risk ratios, risk differences,
#' and arcsine h are first-order local/asymptotic transformations around the
#' pooled independence margins, not exact finite-sample parameterizations away
#' from the independence null. Because the underlying chi-square BFF is
#' sign-blind, signed transformed scales use the sign of the supplied
#' \code{omega} or \code{omega_sequence} to select a branch; odds/risk ratios
#' below 1 select the negative branch and values above 1 select the positive
#' branch. Prior and posterior density plots use the fitted signed branch when
#' available and otherwise split density symmetrically over both signs.
#'
#' @param x vector of two success counts.
#' @param n vector of two group sample sizes.
#' @param LRT should the likelihood-ratio chi-square statistic be used? Default is \code{FALSE}.
#' @param omega prior-mode effect size on the scale selected by \code{effect_size}.
#' @param omega_sequence sequence of prior-mode effect sizes. If no \code{omega}
#' is provided, the default sequence is chosen for the selected effect-size scale.
#' In that case, \code{log_bf_h1} and \code{omega_h1} report the largest log
#' Bayes factor and corresponding effect size on this evaluated grid.
#' @param r variable controlling dispersion of non-local priors. Default is 1. r must be >= 1.
#' @param effect_size scale used for \code{omega} and \code{omega_sequence}. Defaults to \code{"logOR"}. Alternatives include \code{"OR"}, \code{"logRR"}, \code{"risk_ratio"}, \code{"risk_difference"}, \code{"arcsine_h"}, \code{"cohens_w"}, and \code{"phi"}.
#'
#' @return Returns an S3 object of class `BFF` (see `BFF.object` for details).
#' @export
#'
#' @examples
#' propBFF <- prop_test_BFF(x = c(30, 20), n = c(100, 100))
#' propBFF
#' plot(propBFF)
#'
prop_test_BFF <- function(
    x,
    n,
    LRT = FALSE,
    omega = NULL,
    omega_sequence = NULL,
    r = 1,
    effect_size = NULL){

  omega_sequence_missing <- is.null(omega_sequence)
  effect_size_supplied   <- !is.null(effect_size)

  input <- .process_input.prop.test(x = x, n = n, LRT = LRT, r = r)
  effect_size <- .effect_size_normalize("prop_test", effect_size)

  if(is.null(omega) && omega_sequence_missing){
    omega_sequence <- .effect_size_default_sequence("prop_test", effect_size)
  }

  .table_BFF_fit(
    input                = input,
    test_type            = "prop_test",
    r                    = r,
    omega                = omega,
    omega_sequence       = omega_sequence,
    effect_size          = effect_size,
    effect_size_supplied = effect_size_supplied,
    store_effect_size    = TRUE
  )
}

#' contingency_table_BFF
#'
#' Constructs a BFF for an independence test from a contingency table. The test
#' statistic is the uncorrected Pearson chi-square statistic, or the
#' likelihood-ratio G-squared statistic when \code{LRT = TRUE}. For 2x2 tables,
#' the function stores the table margins required for first-order
#' local/asymptotic transformations to odds-ratio, risk-ratio, risk-difference,
#' and arcsine-h scales; these are not exact finite-sample parameterizations
#' away from the independence null. For signed transformed scales, the
#' chi-square BFF is sign-blind: the sign of the supplied \code{omega} or
#' \code{omega_sequence} selects the branch, odds/risk ratios below 1 select
#' the negative branch, and values above 1 select the positive branch. Prior
#' and posterior density plots use the fitted signed branch when available and
#' otherwise split density symmetrically over both signs. For larger tables,
#' use the usual chi-square association scales.
#'
#' @param table matrix-like object of nonnegative integer counts.
#' @param LRT should the likelihood-ratio chi-square statistic be used? Default is \code{FALSE}.
#' @param omega prior-mode effect size on the scale selected by \code{effect_size}.
#' @param omega_sequence sequence of prior-mode effect sizes. If no \code{omega}
#' is provided, the default sequence is chosen for the selected effect-size scale.
#' In that case, \code{log_bf_h1} and \code{omega_h1} report the largest log
#' Bayes factor and corresponding effect size on this evaluated grid.
#' @param r variable controlling dispersion of non-local priors. Default is 1. r must be >= 1.
#' @param effect_size scale used for \code{omega} and \code{omega_sequence}. Defaults to the package's internal \code{omega} RMSES scale. Alternatives include \code{"cohens_w"}, \code{"phi"}, \code{"cramers_v"}, \code{"tschuprow_t"}, \code{"contingency_coefficient"}, and, for 2x2 tables, \code{"logOR"}, \code{"OR"}, \code{"logRR"}, \code{"risk_ratio"}, \code{"risk_difference"}, and \code{"arcsine_h"}.
#'
#' @return Returns an S3 object of class `BFF` (see `BFF.object` for details).
#' @export
#'
#' @examples
#' tab <- matrix(c(30, 70, 20, 80), nrow = 2, byrow = TRUE)
#' tabBFF <- contingency_table_BFF(tab, effect_size = "logOR")
#' tabBFF
#' plot(tabBFF)
#'
contingency_table_BFF <- function(
    table,
    LRT = FALSE,
    omega = NULL,
    omega_sequence = NULL,
    r = 1,
    effect_size = NULL){

  omega_sequence_missing <- is.null(omega_sequence)
  effect_size_supplied   <- !is.null(effect_size)

  input <- .process_input.contingency.table(table = table, LRT = LRT, r = r)
  effect_size <- .effect_size_normalize("contingency_table", effect_size)

  if(is.null(omega) && omega_sequence_missing){
    omega_sequence <- .effect_size_default_sequence("contingency_table", effect_size)
  }

  .table_BFF_fit(
    input                = input,
    test_type            = "contingency_table",
    r                    = r,
    omega                = omega,
    omega_sequence       = omega_sequence,
    effect_size          = effect_size,
    effect_size_supplied = effect_size_supplied,
    store_effect_size    = effect_size_supplied
  )
}

.table_BFF_fit <- function(
    input,
    test_type,
    r,
    omega = NULL,
    omega_sequence,
    effect_size,
    effect_size_supplied,
    store_effect_size = FALSE){

  if(store_effect_size || effect_size_supplied){
    input$effect_size <- effect_size
  }

  omega_input <- unname(if(!is.null(omega)) omega else omega_sequence)
  omega_sign <- unname(.effect_size_branch_sign(
    value       = omega_input,
    test_type   = test_type,
    effect_size = effect_size,
    input       = input
  ))

  omega_internal <- unname(.effect_size_to_internal(
    value         = omega_input,
    test_type     = test_type,
    effect_size   = effect_size,
    input         = input,
    table_dim     = input$table_dim,
    table_margins = input$table_margins
  ))
  tau2 <- lapply(omega_input, function(x){
    unname(.effect_size_prior_mode_tau2(
      value         = x,
      test_type     = test_type,
      effect_size   = effect_size,
      input         = input,
      r             = r,
      table_dim     = input$table_dim,
      table_margins = input$table_margins
    ))
  })

  results <- backend_chi2(
    input = input,
    r     = r,
    omega = omega_internal,
    tau2  = tau2
  )

  if(is.null(omega)){
    cutoff <- if(effect_size_supplied || store_effect_size) .effect_size_minimum_cutoff(
        test_type     = test_type,
        effect_size   = effect_size,
        input         = input,
        table_dim     = input$table_dim,
        table_margins = input$table_margins,
        default       = 0.1
      ) else 0.1
    idx_min <- get_min_omega_bff_index(
      omega  = omega_internal,
      bff    = results,
      cutoff = cutoff
    )
    if(is.na(idx_min)){
      minimums <- c(NA_real_, NA_real_)
      minimum_sign <- NA_real_
    }else{
      minimums <- c(results[idx_min], omega_internal[idx_min])
      minimum_sign <- omega_sign[idx_min]
    }
  }else{
    minimums <- c(NULL, NULL)
    minimum_sign <- NULL
  }

  if(is.null(omega)){
    log_bf         <- c(0, results)
    omega_internal <- c(0, omega_internal)
    omega_sign     <- c(1, omega_sign)
    tau2_output    <- c(list(rep(0, length(input$chi2_stat))), tau2)
    idx_max        <- which.max(log_bf)
    this_log_bf    <- log_bf[idx_max]
    this_omega     <- omega_internal[idx_max]
    this_sign      <- omega_sign[idx_max]
    this_tau2      <- tau2_output[[idx_max]]
  }else{
    this_log_bf    <- results
    this_omega     <- omega_internal
    this_sign      <- omega_sign
    this_tau2      <- tau2
  }

  output <- list(
    log_bf_h1   = this_log_bf,
    omega_h1    = this_omega,
    log_bf_h0   = minimums[1],
    omega_h0    = minimums[2],
    effect_size_sign_h1 = this_sign,
    effect_size_sign_h0 = minimum_sign,
    omega_set   = !is.null(omega),
    tau2_h1     = this_tau2,
    test_type   = test_type,
    generic_test = FALSE,
    r           = r,
    input       = input
  )
  if(is.null(omega)){
    output$BFF <- list(log_bf = log_bf, omega = omega_internal, effect_size_sign = omega_sign, tau2 = tau2_output)
  }

  class(output) <- "BFF"
  output
}

.process_input.prop.test <- function(x, n, LRT, r){
  .check_r(r)

  x <- as.numeric(x)
  n <- as.numeric(n)

  if(length(x) != 2 || length(n) != 2)
    stop("`x` and `n` must be numeric vectors of length two.")
  if(any(!is.finite(x)) || any(!is.finite(n)))
    stop("`x` and `n` must contain finite values.")
  if(any(x < 0) || any(n <= 0) || any(x > n))
    stop("`x` must contain success counts between 0 and `n`.")
  if(any(x != floor(x)) || any(n != floor(n)))
    stop("`x` and `n` must contain integer counts.")

  count_table <- cbind(success = x, failure = n - x)
  rownames(count_table) <- paste0("group", seq_along(x))

  input <- .process_input.contingency.table(table = count_table, LRT = LRT, r = r)
  input$x <- x
  input$n_by_group <- n
  input
}

.process_input.contingency.table <- function(table, LRT, r){
  .check_r(r)

  if(!is.logical(LRT) || length(LRT) != 1 || is.na(LRT))
    stop("`LRT` must be TRUE or FALSE.")

  count_table <- as.matrix(table)
  storage.mode(count_table) <- "numeric"

  if(length(dim(count_table)) != 2 || any(dim(count_table) < 2))
    stop("`table` must have at least two rows and two columns.")
  if(any(!is.finite(count_table)))
    stop("`table` must contain finite counts.")
  if(any(count_table < 0))
    stop("`table` counts must be nonnegative.")
  if(any(count_table != floor(count_table)))
    stop("`table` counts must be integers.")

  n <- sum(count_table)
  if(n <= 0)
    stop("`table` must contain at least one count.")

  row_totals <- rowSums(count_table)
  col_totals <- colSums(count_table)
  if(any(row_totals <= 0) || any(col_totals <= 0))
    stop("`table` must not contain empty rows or columns.")

  expected <- outer(row_totals, col_totals) / n
  pearson <- sum((count_table - expected)^2 / expected)
  g2 <- .likelihood_ratio_chisq(count_table, expected)
  df <- prod(dim(count_table) - 1)

  input <- list(
    chi2_stat = if(isTRUE(LRT)) g2 else pearson,
    n         = n,
    df        = df,
    LRT       = isTRUE(LRT),
    table     = count_table,
    table_dim = dim(count_table)
  )

  if(all(dim(count_table) == c(2, 2))){
    input$table_margins <- c(row1 = row_totals[1] / n, col1 = col_totals[1] / n)
    input$observed_effects <- .table_2x2_observed_effects(count_table)
  }

  input
}

.likelihood_ratio_chisq <- function(observed, expected){
  positive <- observed > 0
  2 * sum(observed[positive] * log(observed[positive] / expected[positive]))
}

.table_2x2_observed_effects <- function(count_table){
  a <- count_table[1, 1]
  b <- count_table[1, 2]
  c <- count_table[2, 1]
  d <- count_table[2, 2]
  n <- sum(count_table)

  row_totals <- rowSums(count_table)
  col_totals <- colSums(count_table)

  p1 <- a / row_totals[1]
  p2 <- c / row_totals[2]
  determinant <- a * d - b * c
  phi_denominator <- sqrt(prod(row_totals) * prod(col_totals))
  signed_phi <- determinant / phi_denominator
  ratio_table <- if(any(count_table == 0)){
    # Haldane-Anscombe correction keeps ratio summaries finite for zero cells.
    count_table + 0.5
  }else{
    count_table
  }
  ratio_row_totals <- rowSums(ratio_table)
  ratio_p1 <- ratio_table[1, 1] / ratio_row_totals[1]
  ratio_p2 <- ratio_table[2, 1] / ratio_row_totals[2]
  ratio_or <- (ratio_table[1, 1] * ratio_table[2, 2]) /
    (ratio_table[1, 2] * ratio_table[2, 1])
  ratio_rr <- ratio_p1 / ratio_p2

  list(
    sign = sign(determinant),
    signed_phi = signed_phi,
    cohens_w = abs(signed_phi),
    logOR = log(ratio_or),
    OR = ratio_or,
    risk_difference = p1 - p2,
    logRR = log(ratio_rr),
    risk_ratio = ratio_rr,
    arcsine_h = 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2)),
    row1 = row_totals[1] / n,
    col1 = col_totals[1] / n
  )
}
