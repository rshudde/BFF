.effect_size_key <- function(effect_size){
  key <- tolower(as.character(effect_size))
  key <- gsub("\\^", "", key)
  key <- gsub("[^a-z0-9]+", "_", key)
  key <- gsub("^_+|_+$", "", key)
  key
}

.effect_size_default <- function(test_type){
  switch(
    test_type,
    "t_test"            = "cohens_d",
    "z_test"            = "cohens_d",
    "chi2_test"         = "omega",
    "contingency_table" = "omega",
    "prop_test"         = "logOR",
    "f_test"            = "omega",
    "regression_test"   = "cohens_f",
    stop("Unknown BFF test type.")
  )
}

.effect_size_chi2_family <- function(test_type){
  test_type %in% c("chi2_test", "contingency_table", "prop_test")
}

.effect_size_normalize_chi2 <- function(key, allow_2x2 = FALSE, label = "chi-square"){
  out <- switch(
    key,
    "omega" = ,
    "rmses" = ,
    "rmsea" = ,
    "root_mean_square_standardized_effect_size" = "omega",
    "cohen_w" = ,
    "cohens_w" = ,
    "w" = "cohens_w",
    "phi" = "phi",
    "cramer_v" = ,
    "cramers_v" = ,
    "cramer" = "cramers_v",
    "tschuprow_t" = ,
    "tschuprows_t" = "tschuprow_t",
    "contingency_coefficient" = ,
    "pearson_contingency_coefficient" = ,
    "pearson_c" = ,
    "pearsons_c" = ,
    "c" = "contingency_coefficient",
    if(allow_2x2) switch(
      key,
      "log_or" = ,
      "log_odds_ratio" = ,
      "logor" = "logOR",
      "odds_ratio" = ,
      "or" = "OR",
      "log_rr" = ,
      "log_risk_ratio" = ,
      "log_relative_risk" = ,
      "logrr" = "logRR",
      "risk_ratio" = ,
      "relative_risk" = ,
      "rr" = "risk_ratio",
      "risk_difference" = ,
      "rd" = "risk_difference",
      "arcsine_h" = ,
      "cohens_h" = ,
      "h" = "arcsine_h",
      NULL
    ) else NULL
  )

  if(is.null(out)){
    stop(sprintf("Unsupported effect size for %s BFF objects.", label))
  }

  out
}

.effect_size_normalize <- function(test_type, effect_size = NULL){
  if(is.null(effect_size) || identical(effect_size, "default")){
    return(.effect_size_default(test_type))
  }

  key <- .effect_size_key(effect_size)

  out <- switch(
    test_type,
    "t_test" = switch(
      key,
      "omega" = ,
      "cohen_d" = ,
      "cohens_d" = ,
      "d" = "cohens_d",
      stop("Unsupported effect size for t-test BFF objects.")
    ),
    "z_test" = switch(
      key,
      "omega" = ,
      "cohen_d" = ,
      "cohens_d" = ,
      "d" = "cohens_d",
      stop("Unsupported effect size for z-test BFF objects.")
    ),
    "chi2_test" = .effect_size_normalize_chi2(key, allow_2x2 = FALSE, label = "chi-square"),
    "contingency_table" = .effect_size_normalize_chi2(key, allow_2x2 = TRUE, label = "contingency-table"),
    "prop_test" = .effect_size_normalize_chi2(key, allow_2x2 = TRUE, label = "proportions-test"),
    "f_test" = switch(
      key,
      "omega" = ,
      "rmses" = ,
      "root_mean_square_standardized_effect_size" = "omega",
      "cohen_f" = ,
      "cohens_f" = ,
      "f" = "cohens_f",
      "cohen_f2" = ,
      "cohens_f2" = ,
      "partial_f2" = ,
      "partial_f_squared" = ,
      "f2" = ,
      "f_squared" = "cohens_f2",
      "partial_eta2" = ,
      "partial_eta_squared" = ,
      "eta2" = ,
      "eta_squared" = "partial_eta2",
      "partial_r2" = ,
      "partial_r_squared" = ,
      "r2" = "partial_r2",
      stop("Unsupported effect size for F-test BFF objects.")
    ),
    "regression_test" = switch(
      key,
      "omega" = ,
      "delta" = ,
      "signed_cohen_f" = ,
      "signed_cohens_f" = ,
      "cohen_f" = ,
      "cohens_f" = ,
      "f" = "cohens_f",
      "partial_r" = ,
      "r" = "partial_r",
      "partial_r2" = ,
      "partial_r_squared" = ,
      "r2" = "partial_r2",
      "cohen_f2" = ,
      "cohens_f2" = ,
      "partial_f2" = ,
      "partial_f_squared" = ,
      "f2" = ,
      "f_squared" = "cohens_f2",
      stop("Unsupported effect size for regression-test BFF objects.")
    ),
    stop("Unknown BFF test type.")
  )

  out
}

.effect_size_for_object <- function(x, effect_size = NULL){
  if(is.null(effect_size) && !is.null(x$input$effect_size)){
    effect_size <- x$input$effect_size
  }
  .effect_size_normalize(x$test_type, effect_size)
}

.effect_size_label <- function(test_type, effect_size = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)

  switch(
    effect_size,
    "cohens_d" = "Cohen's d",
    "omega" = switch(
      test_type,
      "chi2_test" = "omega (RMSES)",
      "contingency_table" = "omega (RMSES)",
      "prop_test" = "omega (RMSES)",
      "f_test" = "omega (RMSES)",
      "omega"
    ),
    "cohens_w" = "Cohen's w",
    "phi" = "phi coefficient",
    "cramers_v" = "Cramer's V",
    "tschuprow_t" = "Tschuprow's T",
    "contingency_coefficient" = "contingency coefficient",
    "logOR" = "log odds ratio",
    "OR" = "odds ratio",
    "logRR" = "log risk ratio",
    "risk_ratio" = "risk ratio",
    "risk_difference" = "risk difference",
    "arcsine_h" = "arcsine h",
    "cohens_f" = switch(
      test_type,
      "regression_test" = "signed Cohen's f",
      "Cohen's f"
    ),
    "cohens_f2" = "Cohen's f^2",
    "partial_eta2" = "partial eta^2",
    "partial_r" = "partial r",
    "partial_r2" = "partial R^2",
    effect_size
  )
}

.effect_size_plot_label <- function(x, effect_size = NULL){
  if(is.null(effect_size) && is.null(x$input$effect_size)){
    return(.test_effect_size_name(x$test_type))
  }
  .effect_size_label(x$test_type, .effect_size_for_object(x, effect_size))
}

.effect_size_table_dim <- function(input = NULL, table_dim = NULL){
  if(is.null(table_dim) && !is.null(input$table_dim)){
    table_dim <- input$table_dim
  }

  if(is.null(table_dim)){
    stop("`table_dim` must be supplied for this chi-square effect-size transformation.")
  }

  if(!is.numeric(table_dim) || length(table_dim) != 2 || any(table_dim <= 1) || any(table_dim != floor(table_dim))){
    stop("`table_dim` must be an integer vector of length two, e.g., c(rows, columns).")
  }

  if(!is.null(input) && !is.null(input$df)){
    expected_df <- prod(table_dim - 1)
    input_df <- unique(as.numeric(input$df))
    if(length(input_df) != 1 || is.na(input_df) || input_df != expected_df){
      stop("`df` must equal (rows - 1) * (columns - 1) for this `table_dim`.")
    }
  }

  table_dim
}

.effect_size_table_margins <- function(input = NULL, table_margins = NULL){
  if(is.null(table_margins) && !is.null(input$table_margins)){
    table_margins <- input$table_margins
  }

  if(is.null(table_margins)){
    stop("`table_margins` must be supplied for this 2x2 effect-size transformation.")
  }

  if(!is.numeric(table_margins) || length(table_margins) != 2 ||
     any(!is.finite(table_margins)) || any(table_margins <= 0) || any(table_margins >= 1)){
    stop("`table_margins` must be two marginal probabilities in (0, 1), e.g., c(row1 = 0.5, col1 = 0.5).")
  }

  as.numeric(table_margins)
}

.effect_size_2x2_scale <- function(effect_size, input = NULL, table_dim = NULL, table_margins = NULL){
  table_dim <- .effect_size_table_dim(input, table_dim)
  if(!all(table_dim == c(2, 2))){
    stop(sprintf("effect_size = \"%s\" is available only for 2x2 tables.", effect_size))
  }

  table_margins <- .effect_size_table_margins(input, table_margins)
  row1 <- table_margins[1]
  col1 <- table_margins[2]
  row2 <- 1 - row1
  col2 <- 1 - col1

  # Local independence-null relationships: phi = scale * transformed_effect + o(effect).
  switch(
    effect_size,
    "logOR" = ,
    "OR" = sqrt(row1 * row2 * col1 * col2),
    "risk_difference" = sqrt(row1 * row2 / (col1 * col2)),
    "logRR" = ,
    "risk_ratio" = sqrt(row1 * row2 * col1 / col2),
    "arcsine_h" = sqrt(row1 * row2),
    stop(sprintf("Unsupported 2x2 effect-size scale for effect_size = \"%s\".", effect_size))
  )
}

.effect_size_check_interval <- function(value, effect_size, lower = 0, upper = Inf, upper_open = FALSE, allow_negative = FALSE){
  value <- as.numeric(value)
  if(any(!is.finite(value))){
    stop(sprintf("`omega` and `omega_sequence` must contain finite values for effect_size = \"%s\".", effect_size))
  }
  finite_value <- value[is.finite(value)]

  if(!allow_negative && any(finite_value < lower)){
    stop(sprintf("`omega` and `omega_sequence` must be >= %s for effect_size = \"%s\".", lower, effect_size))
  }

  if(allow_negative && any(abs(finite_value) >= upper)){
    stop(sprintf("Absolute values for `omega` and `omega_sequence` must be < %s for effect_size = \"%s\".", upper, effect_size))
  }

  if(!allow_negative && is.finite(upper)){
    bad_upper <- if(upper_open) finite_value >= upper else finite_value > upper
    if(any(bad_upper)){
      relation <- if(upper_open) "<" else "<="
      stop(sprintf("`omega` and `omega_sequence` must be %s %s for effect_size = \"%s\".", relation, upper, effect_size))
    }
  }

  invisible(value)
}

.effect_size_check_positive <- function(value, effect_size){
  value <- as.numeric(value)
  if(any(!is.finite(value))){
    stop(sprintf("`omega` and `omega_sequence` must contain finite values for effect_size = \"%s\".", effect_size))
  }
  finite_value <- value[is.finite(value)]
  if(any(finite_value <= 0)){
    stop(sprintf("`omega` and `omega_sequence` must be > 0 for effect_size = \"%s\".", effect_size))
  }
  invisible(value)
}

.effect_size_default_sequence <- function(test_type, effect_size = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)

  switch(
    effect_size,
    "partial_eta2" = ,
    "partial_r" = ,
    "partial_r2" = ,
    "contingency_coefficient" = seq(0.01, 0.99, by = 0.01),
    "logOR" = seq(0.01, 4, by = 0.01),
    "OR" = exp(seq(0.01, 1.99, by = 0.01)),
    "logRR" = seq(0.01, 3, by = 0.01),
    "risk_ratio" = exp(seq(0.01, 1.99, by = 0.01)),
    "risk_difference" = seq(0.01, 0.99, by = 0.01),
    "arcsine_h" = seq(0.01, 0.99 * pi, length.out = 400),
    seq(0.01, 1, by = 0.01)
  )
}

.effect_size_check_common_transform_scale <- function(test_type, effect_size = NULL, input = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)

  if(.effect_size_chi2_family(test_type) &&
     effect_size %in% c("cohens_w", "phi", "cramers_v", "tschuprow_t", "contingency_coefficient") &&
     length(unique(input$df)) > 1){
    stop(sprintf(
      "effect_size = \"%s\" requires a common `df` for vectorized chi-square inputs. Use effect_size = \"omega\" or fit the studies separately.",
      effect_size
    ))
  }

  if(test_type == "f_test" && effect_size != "omega" && length(unique(input$df1)) > 1){
    stop(sprintf(
      "effect_size = \"%s\" requires a common `df1` for vectorized F-test inputs. Use effect_size = \"omega\" or fit the studies separately.",
      effect_size
    ))
  }

  invisible(TRUE)
}

.effect_size_minimum_cutoff <- function(test_type, effect_size = NULL, input = NULL, table_dim = NULL, table_margins = NULL, default){
  cutpoints <- .get_effect_size_cutpoints(test_type, effect_size)

  if(length(cutpoints) == 0){
    return(default)
  }

  .effect_size_to_internal(
    value       = cutpoints[1],
    test_type   = test_type,
    effect_size = effect_size,
    input       = input,
    table_dim   = table_dim,
    table_margins = table_margins
  )
}

.effect_size_nonzero_sign <- function(value, default = 1){
  value <- as.numeric(value)
  out <- rep(default, length(value))
  finite <- is.finite(value)
  out[finite & value < 0] <- -1
  out[finite & value > 0] <- 1
  out[!is.finite(out) | out == 0] <- 1
  out
}

.effect_size_branch_sign <- function(value, test_type, effect_size = NULL, input = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)
  value <- as.numeric(value)
  out <- rep(1, length(value))

  if(.effect_size_chi2_family(test_type)){
    return(switch(
      effect_size,
      "logOR" = ,
      "logRR" = ,
      "risk_difference" = ,
      "arcsine_h" = .effect_size_nonzero_sign(value),
      "OR" = ,
      "risk_ratio" = {
        out[is.finite(value) & value < 1] <- -1
        out
      },
      out
    ))
  }

  if(test_type == "regression_test" && effect_size %in% c("cohens_f", "partial_r")){
    alternative <- input$alternative.original
    if(identical(alternative, "less")){
      return(rep(-1, length(value)))
    }
    if(identical(alternative, "greater")){
      return(rep(1, length(value)))
    }
    return(.effect_size_nonzero_sign(value))
  }

  out
}

.effect_size_recycle_branch_sign <- function(branch_sign, value){
  value <- as.numeric(value)
  if(is.null(branch_sign)){
    return(rep(1, length(value)))
  }

  branch_sign <- as.numeric(branch_sign)
  branch_sign[!is.finite(branch_sign) | branch_sign == 0] <- 1

  if(length(branch_sign) == 1){
    return(rep(branch_sign, length(value)))
  }
  if(length(branch_sign) == length(value)){
    return(branch_sign)
  }

  stop("Internal effect-size branch sign length mismatch.")
}

.effect_size_to_internal <- function(value, test_type, effect_size = NULL, input = NULL, table_dim = NULL, table_margins = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)
  value <- as.numeric(value)

  switch(
    test_type,
    "t_test" = ,
    "z_test" = {
      .effect_size_check_interval(value, effect_size, allow_negative = TRUE, upper = Inf)
      abs(value)
    },
    "chi2_test" = ,
    "contingency_table" = ,
    "prop_test" = {
      df <- unique(input$df)
      if(length(df) != 1) df <- input$df
      switch(
        effect_size,
        "omega" = {
          .effect_size_check_interval(value, effect_size)
          value
        },
        "cohens_w" = ,
        "phi" = {
          .effect_size_check_interval(value, effect_size)
          value / sqrt(df)
        },
        "cramers_v" = {
          table_dim <- .effect_size_table_dim(input, table_dim)
          .effect_size_check_interval(value, effect_size)
          value * sqrt(min(table_dim - 1)) / sqrt(df)
        },
        "tschuprow_t" = {
          table_dim <- .effect_size_table_dim(input, table_dim)
          .effect_size_check_interval(value, effect_size)
          value * ((table_dim[1] - 1) * (table_dim[2] - 1))^(1/4) / sqrt(df)
        },
        "contingency_coefficient" = {
          .effect_size_check_interval(value, effect_size, upper = 1, upper_open = TRUE)
          value / (sqrt(df) * sqrt(1 - value^2))
        },
        "logOR" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          .effect_size_check_interval(value, effect_size, allow_negative = TRUE, upper = Inf)
          abs(value) * scale
        },
        "OR" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          .effect_size_check_positive(value, effect_size)
          abs(log(value)) * scale
        },
        "logRR" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          .effect_size_check_interval(value, effect_size, allow_negative = TRUE, upper = Inf)
          abs(value) * scale
        },
        "risk_ratio" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          .effect_size_check_positive(value, effect_size)
          abs(log(value)) * scale
        },
        "risk_difference" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          .effect_size_check_interval(value, effect_size, allow_negative = TRUE, upper = 1)
          abs(value) * scale
        },
        "arcsine_h" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          .effect_size_check_interval(value, effect_size, allow_negative = TRUE, upper = pi)
          abs(value) * scale
        }
      )
    },
    "f_test" = {
      df1 <- unique(input$df1)
      if(length(df1) != 1) df1 <- input$df1
      switch(
        effect_size,
        "omega" = {
          .effect_size_check_interval(value, effect_size)
          value
        },
        "cohens_f" = {
          .effect_size_check_interval(value, effect_size)
          value * sqrt(2 / df1)
        },
        "cohens_f2" = {
          .effect_size_check_interval(value, effect_size)
          sqrt(2 * value / df1)
        },
        "partial_eta2" = ,
        "partial_r2" = {
          .effect_size_check_interval(value, effect_size, upper = 1, upper_open = TRUE)
          sqrt(2 * value / (df1 * (1 - value)))
        }
      )
    },
    "regression_test" = {
      switch(
        effect_size,
        "cohens_f" = {
          .effect_size_check_interval(value, effect_size, allow_negative = TRUE, upper = Inf)
          abs(value)
        },
        "partial_r" = {
          .effect_size_check_interval(value, effect_size, allow_negative = TRUE, upper = 1)
          abs(value) / sqrt(1 - value^2)
        },
        "partial_r2" = {
          .effect_size_check_interval(value, effect_size, upper = 1, upper_open = TRUE)
          sqrt(value / (1 - value))
        },
        "cohens_f2" = {
          .effect_size_check_interval(value, effect_size)
          sqrt(value)
        }
      )
    },
    stop("Unknown BFF test type.")
  )
}

.effect_size_prior_mode_tau2 <- function(value, test_type, effect_size = NULL, input = NULL, r, table_dim = NULL, table_margins = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)
  value <- as.numeric(value)

  switch(
    test_type,
    "t_test" = ,
    "z_test" = {
      if(input$one_sample){
        return(get_one_sample_tau2(n = input$n, w = abs(value), r = r))
      }
      return(get_two_sample_tau2(n1 = input$n1, n2 = input$n2, w = abs(value), r = r))
    },
    "chi2_test" = ,
    "contingency_table" = ,
    "prop_test" = {
      df <- unique(input$df)
      if(length(df) != 1) df <- input$df
      lambda_scale <- input$n * df
      alpha <- df / 2 + r

      switch(
        effect_size,
        "omega" = ,
        "cohens_w" = ,
        "phi" = ,
        "cramers_v" = ,
        "tschuprow_t" = ,
        "logOR" = ,
        "logRR" = ,
        "risk_difference" = ,
        "arcsine_h" = {
          internal <- .effect_size_to_internal(
            value = value, test_type = test_type, effect_size = effect_size,
            input = input, table_dim = table_dim, table_margins = table_margins
          )
          lambda_scale * internal^2 / (2 * alpha - 1)
        },
        "contingency_coefficient" = {
          .effect_size_check_interval(value, effect_size, upper = 1, upper_open = TRUE)
          ifelse(value == 0, 0, {
            b <- (2 * alpha - 1) * (1 - value^2)^2 / (2 * value^2) +
              (alpha + 1) * (1 - value^2)
            input$n / (2 * b)
          })
        },
        "OR" = ,
        "risk_ratio" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          .effect_size_check_positive(value, effect_size)
          log_value <- log(value)
          A <- ifelse(log_value == 0, Inf, ((2 * alpha - 1) / log_value - 1) / (2 * log_value))
          if(any(is.finite(A) & A <= 0)){
            stop(sprintf(
              "No finite non-local prior can have its mode at this value on the %s scale; use the corresponding log scale or a value closer to 1.",
              .effect_size_label(test_type, effect_size)
            ))
          }
          ifelse(is.infinite(A), 0, lambda_scale * scale^2 / (2 * A))
        }
      )
    },
    "f_test" = {
      df1 <- unique(input$df1)
      if(length(df1) != 1) df1 <- input$df1
      alpha <- df1 / 2 + r
      lambda_scale <- input$n * df1 / 2

      switch(
        effect_size,
        "omega" = ,
        "cohens_f" = {
          internal <- .effect_size_to_internal(value, "f_test", effect_size, input)
          lambda_scale * internal^2 / (2 * alpha - 1)
        },
        "cohens_f2" = {
          .effect_size_check_interval(value, effect_size)
          input$n * value / (2 * (alpha - 1))
        },
        "partial_eta2" = ,
        "partial_r2" = {
          .effect_size_check_interval(value, effect_size, upper = 1, upper_open = TRUE)
          ifelse(value == 0, 0, {
            b <- (alpha - 1) * (1 - value)^2 / value + (alpha + 1) * (1 - value)
            input$n / (2 * b)
          })
        }
      )
    },
    "regression_test" = {
      scale2 <- .regression_test_df(n = input$n, k = input$k)

      switch(
        effect_size,
        "cohens_f" = {
          .effect_size_check_interval(value, effect_size, allow_negative = TRUE, upper = Inf)
          scale2 * abs(value)^2 / (2 * r)
        },
        "partial_r" = {
          .effect_size_check_interval(value, effect_size, allow_negative = TRUE, upper = 1)
          abs_value <- abs(value)
          ifelse(abs_value == 0, 0, {
            a <- r * (1 - abs_value^2)^2 / abs_value^2 +
              (r + 3 / 2) * (1 - abs_value^2)
            scale2 / (2 * a)
          })
        },
        "partial_r2" = {
          .effect_size_check_interval(value, effect_size, upper = 1, upper_open = TRUE)
          ifelse(value == 0, 0, {
            a <- (r - 1 / 2) * (1 - value)^2 / value +
              (r + 3 / 2) * (1 - value)
            scale2 / (2 * a)
          })
        },
        "cohens_f2" = {
          .effect_size_check_interval(value, effect_size)
          scale2 * value / (2 * (r - 1 / 2))
        }
      )
    },
    stop("Unknown BFF test type.")
  )
}

.effect_size_from_internal <- function(value, test_type, effect_size = NULL, input = NULL, table_dim = NULL, table_margins = NULL, branch_sign = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)
  value <- as.numeric(value)
  branch_sign <- .effect_size_recycle_branch_sign(branch_sign, value)

  switch(
    test_type,
    "t_test" = ,
    "z_test" = value,
    "chi2_test" = ,
    "contingency_table" = ,
    "prop_test" = {
      df <- unique(input$df)
      if(length(df) != 1) df <- input$df
      switch(
        effect_size,
        "omega" = value,
        "cohens_w" = ,
        "phi" = sqrt(df) * value,
        "cramers_v" = {
          table_dim <- .effect_size_table_dim(input, table_dim)
          sqrt(df) * value / sqrt(min(table_dim - 1))
        },
        "tschuprow_t" = {
          table_dim <- .effect_size_table_dim(input, table_dim)
          sqrt(df) * value / ((table_dim[1] - 1) * (table_dim[2] - 1))^(1/4)
        },
        "contingency_coefficient" = {
          w <- sqrt(df) * value
          w / sqrt(1 + w^2)
        },
        "logOR" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          branch_sign * value / scale
        },
        "OR" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          exp(branch_sign * value / scale)
        },
        "logRR" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          branch_sign * value / scale
        },
        "risk_ratio" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          exp(branch_sign * value / scale)
        },
        "risk_difference" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          branch_sign * value / scale
        },
        "arcsine_h" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          branch_sign * value / scale
        }
      )
    },
    "f_test" = {
      df1 <- unique(input$df1)
      if(length(df1) != 1) df1 <- input$df1
      switch(
      effect_size,
      "omega" = value,
      "cohens_f" = sqrt(df1 / 2) * value,
      "cohens_f2" = df1 * value^2 / 2,
      "partial_eta2" = ,
      "partial_r2" = {
        f2 <- df1 * value^2 / 2
        f2 / (1 + f2)
      }
    )},
    "regression_test" = switch(
      effect_size,
      "cohens_f" = branch_sign * value,
      "partial_r" = branch_sign * value / sqrt(1 + value^2),
      "partial_r2" = value^2 / (1 + value^2),
      "cohens_f2" = value^2
    ),
    stop("Unknown BFF test type.")
  )
}

.effect_size_internal_branches <- function(x, test_type, effect_size = NULL, input = NULL, alternative = NULL, table_dim = NULL, table_margins = NULL, branch_sign = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)
  x <- as.numeric(x)

  make_branch <- function(effect_size, jacobian, valid){
    list(effect_size = effect_size, jacobian = jacobian, valid = valid)
  }

  if(.effect_size_chi2_family(test_type)){
    df <- unique(input$df)
    if(length(df) != 1) df <- input$df
    if(!is.null(branch_sign) && effect_size %in% c("logOR", "OR", "logRR", "risk_ratio", "risk_difference", "arcsine_h")){
      branch_sign <- .effect_size_recycle_branch_sign(branch_sign, x)
      branch_sign <- ifelse(branch_sign < 0, -1, 1)
      return(switch(
        effect_size,
        "logOR" = ,
        "logRR" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          valid <- ifelse(branch_sign > 0, x >= 0, x <= 0)
          make_branch(abs(x) * scale, rep(scale, length(x)), valid)
        },
        "OR" = ,
        "risk_ratio" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          log_x <- log(x)
          valid <- x > 0 & ifelse(branch_sign > 0, log_x >= 0, log_x <= 0)
          make_branch(abs(log_x) * scale, scale / x, valid)
        },
        "risk_difference" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          valid <- abs(x) < 1 & ifelse(branch_sign > 0, x >= 0, x <= 0)
          make_branch(abs(x) * scale, rep(scale, length(x)), valid)
        },
        "arcsine_h" = {
          scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
          valid <- abs(x) < pi & ifelse(branch_sign > 0, x >= 0, x <= 0)
          make_branch(abs(x) * scale, rep(scale, length(x)), valid)
        }
      ))
    }
    return(switch(
      effect_size,
      "omega" = make_branch(x, rep(1, length(x)), x >= 0),
      "cohens_w" = ,
      "phi" = make_branch(x / sqrt(df), rep(1 / sqrt(df), length(x)), x >= 0),
      "cramers_v" = {
        table_dim <- .effect_size_table_dim(input, table_dim)
        jacobian <- sqrt(min(table_dim - 1)) / sqrt(df)
        make_branch(x * jacobian, rep(jacobian, length(x)), x >= 0)
      },
      "tschuprow_t" = {
        table_dim <- .effect_size_table_dim(input, table_dim)
        jacobian <- ((table_dim[1] - 1) * (table_dim[2] - 1))^(1/4) / sqrt(df)
        make_branch(x * jacobian, rep(jacobian, length(x)), x >= 0)
      },
      "contingency_coefficient" = {
        valid <- x >= 0 & x < 1
        inverse <- x / (sqrt(df) * sqrt(1 - x^2))
        jacobian <- 1 / (sqrt(df) * (1 - x^2)^(3/2))
        make_branch(inverse, jacobian, valid)
      },
      "logOR" = {
        scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
        make_branch(abs(x) * scale, rep(scale / 2, length(x)), rep(TRUE, length(x)))
      },
      "OR" = {
        scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
        make_branch(abs(log(x)) * scale, scale / (2 * x), x > 0)
      },
      "logRR" = {
        scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
        make_branch(abs(x) * scale, rep(scale / 2, length(x)), rep(TRUE, length(x)))
      },
      "risk_ratio" = {
        scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
        make_branch(abs(log(x)) * scale, scale / (2 * x), x > 0)
      },
      "risk_difference" = {
        scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
        make_branch(abs(x) * scale, rep(scale / 2, length(x)), abs(x) < 1)
      },
      "arcsine_h" = {
        scale <- .effect_size_2x2_scale(effect_size, input, table_dim, table_margins)
        make_branch(abs(x) * scale, rep(scale / 2, length(x)), abs(x) < pi)
      }
    ))
  }

  if(test_type == "f_test"){
    df1 <- unique(input$df1)
    if(length(df1) != 1) df1 <- input$df1
    return(switch(
      effect_size,
      "omega" = make_branch(x, rep(1, length(x)), x >= 0),
      "cohens_f" = make_branch(x * sqrt(2 / df1), rep(sqrt(2 / df1), length(x)), x >= 0),
      "cohens_f2" = {
        valid <- x > 0
        make_branch(sqrt(2 * x / df1), 1 / (sqrt(2 * df1) * sqrt(x)), valid)
      },
      "partial_eta2" = ,
      "partial_r2" = {
        valid <- x > 0 & x < 1
        make_branch(
          sqrt(2 * x / (df1 * (1 - x))),
          1 / (sqrt(2 * df1) * sqrt(x) * (1 - x)^(3/2)),
          valid
        )
      }
    ))
  }

  if(test_type == "regression_test"){
    one_sided <- !is.null(alternative) && alternative != "two.sided"
    sign <- if(identical(alternative, "less")) -1 else 1

    return(switch(
      effect_size,
      "cohens_f" = make_branch(sign * x, rep(1, length(x)), rep(TRUE, length(x))),
      "partial_r" = {
        valid <- abs(x) < 1
        delta <- x / sqrt(1 - x^2)
        make_branch(sign * delta, 1 / (1 - x^2)^(3/2), valid)
      },
      "cohens_f2" = {
        valid <- x > 0
        delta <- sqrt(x)
        jacobian <- 1 / (2 * sqrt(x))
        if(one_sided){
          list(make_branch(delta, jacobian, valid))
        }else{
          list(
            make_branch( delta, jacobian, valid),
            make_branch(-delta, jacobian, valid)
          )
        }
      },
      "partial_r2" = {
        valid <- x > 0 & x < 1
        delta <- sqrt(x / (1 - x))
        jacobian <- 1 / (2 * sqrt(x) * (1 - x)^(3/2))
        if(one_sided){
          list(make_branch(delta, jacobian, valid))
        }else{
          list(
            make_branch( delta, jacobian, valid),
            make_branch(-delta, jacobian, valid)
          )
        }
      }
    ))
  }

  stop("Unknown BFF test type.")
}

.effect_size_bounded_local_upper <- function(effect_size){
  switch(
    effect_size,
    "risk_difference" = 1,
    "arcsine_h"       = pi,
    Inf
  )
}

.effect_size_density_normalizer <- function(test_type, effect_size = NULL, density, input = NULL, table_dim = NULL, table_margins = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)
  if(!.effect_size_chi2_family(test_type) || !effect_size %in% c("risk_difference", "arcsine_h")){
    return(1)
  }

  scale <- .effect_size_2x2_scale(
    effect_size   = effect_size,
    input         = input,
    table_dim     = table_dim,
    table_margins = table_margins
  )
  upper <- scale * .effect_size_bounded_local_upper(effect_size)

  if(!is.finite(upper) || upper <= 0){
    return(1)
  }

  normalizer <- tryCatch(
    stats::integrate(
      f = function(omega) density(omega),
      lower = 0,
      upper = upper,
      rel.tol = 1e-7,
      subdivisions = 1000
    )$value,
    error = function(e) NA_real_
  )

  if(!is.finite(normalizer) || normalizer <= 0){
    return(1)
  }

  normalizer
}

.effect_size_density <- function(test_type, effect_size = NULL, x, density, input = NULL, alternative = NULL, table_dim = NULL, table_margins = NULL, branch_sign = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)
  branches <- .effect_size_internal_branches(
    x           = x,
    test_type   = test_type,
    effect_size = effect_size,
    input       = input,
    alternative = alternative,
    table_dim   = table_dim,
    table_margins = table_margins,
    branch_sign = branch_sign
  )

  if(!is.list(branches[[1]])){
    branches <- list(branches)
  }

  out <- numeric(length(x))

  for(branch in branches){
    valid <- branch$valid & is.finite(branch$effect_size) & is.finite(branch$jacobian)
    contribution <- numeric(length(x))
    if(any(valid)){
      contribution[valid] <- density(branch$effect_size[valid]) * branch$jacobian[valid]
      contribution[!is.finite(contribution)] <- 0
    }
    out <- out + contribution
  }

  out[!is.finite(out)] <- 0
  out / .effect_size_density_normalizer(
    test_type     = test_type,
    effect_size   = effect_size,
    density       = density,
    input         = input,
    table_dim     = table_dim,
    table_margins = table_margins
  )
}

.effect_size_default_x_limit <- function(test_type, effect_size = NULL, alternative = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)

  if(.effect_size_chi2_family(test_type)){
    return(switch(
      effect_size,
      "contingency_coefficient" = c(0, 0.99),
      "logOR" = c(-4, 4),
      "OR" = exp(c(-4, 4)),
      "logRR" = c(-3, 3),
      "risk_ratio" = exp(c(-3, 3)),
      "risk_difference" = c(-0.99, 0.99),
      "arcsine_h" = c(-0.99 * pi, 0.99 * pi),
      c(0, 3)
    ))
  }

  if(test_type == "f_test"){
    return(switch(
      effect_size,
      "partial_eta2" = ,
      "partial_r2" = c(0, 0.99),
      "cohens_f2" = c(0, 9),
      c(0, 3)
    ))
  }

  if(test_type == "regression_test"){
    if(effect_size == "partial_r"){
      return(switch(
        alternative,
        "two.sided" = c(-0.99, 0.99),
        "greater"   = c(0, 0.99),
        "less"      = c(-0.99, 0)
      ))
    }
    if(effect_size == "partial_r2"){
      return(c(0, 0.99))
    }
    if(effect_size == "cohens_f2"){
      return(c(0, 9))
    }
  }

  NULL
}
