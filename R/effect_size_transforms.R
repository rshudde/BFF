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
    "t_test"          = "cohens_d",
    "z_test"          = "cohens_d",
    "chi2_test"       = "omega",
    "f_test"          = "omega",
    "regression_test" = "cohens_f",
    stop("Unknown BFF test type.")
  )
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
    "chi2_test" = switch(
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
      stop("Unsupported effect size for chi-square BFF objects.")
    ),
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
      "f_test" = "omega (RMSES)",
      "omega"
    ),
    "cohens_w" = "Cohen's w",
    "phi" = "phi coefficient",
    "cramers_v" = "Cramer's V",
    "tschuprow_t" = "Tschuprow's T",
    "contingency_coefficient" = "contingency coefficient",
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

.effect_size_check_interval <- function(value, effect_size, lower = 0, upper = Inf, upper_open = FALSE, allow_negative = FALSE){
  value <- as.numeric(value)
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

.effect_size_default_sequence <- function(test_type, effect_size = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)

  switch(
    effect_size,
    "partial_eta2" = ,
    "partial_r" = ,
    "partial_r2" = ,
    "contingency_coefficient" = seq(0.01, 0.99, by = 0.01),
    seq(0.01, 1, by = 0.01)
  )
}

.effect_size_minimum_cutoff <- function(test_type, effect_size = NULL, input = NULL, table_dim = NULL, default){
  cutpoints <- .get_effect_size_cutpoints(test_type, effect_size)

  if(length(cutpoints) == 0){
    return(default)
  }

  .effect_size_to_internal(
    value       = cutpoints[1],
    test_type   = test_type,
    effect_size = effect_size,
    input       = input,
    table_dim   = table_dim
  )
}

.effect_size_to_internal <- function(value, test_type, effect_size = NULL, input = NULL, table_dim = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)
  value <- as.numeric(value)

  switch(
    test_type,
    "t_test" = ,
    "z_test" = {
      .effect_size_check_interval(value, effect_size, allow_negative = TRUE, upper = Inf)
      abs(value)
    },
    "chi2_test" = {
      df <- input$df
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
        }
      )
    },
    "f_test" = {
      df1 <- input$df1
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

.effect_size_from_internal <- function(value, test_type, effect_size = NULL, input = NULL, table_dim = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)
  value <- as.numeric(value)

  switch(
    test_type,
    "t_test" = ,
    "z_test" = value,
    "chi2_test" = {
      df <- input$df
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
        }
      )
    },
    "f_test" = switch(
      effect_size,
      "omega" = value,
      "cohens_f" = sqrt(input$df1 / 2) * value,
      "cohens_f2" = input$df1 * value^2 / 2,
      "partial_eta2" = ,
      "partial_r2" = {
        f2 <- input$df1 * value^2 / 2
        f2 / (1 + f2)
      }
    ),
    "regression_test" = switch(
      effect_size,
      "cohens_f" = value,
      "partial_r" = value / sqrt(1 + value^2),
      "partial_r2" = value^2 / (1 + value^2),
      "cohens_f2" = value^2
    ),
    stop("Unknown BFF test type.")
  )
}

.effect_size_internal_branches <- function(x, test_type, effect_size = NULL, input = NULL, alternative = NULL, table_dim = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)
  x <- as.numeric(x)

  make_branch <- function(effect_size, jacobian, valid){
    list(effect_size = effect_size, jacobian = jacobian, valid = valid)
  }

  if(test_type == "chi2_test"){
    df <- input$df
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
      }
    ))
  }

  if(test_type == "f_test"){
    df1 <- input$df1
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

.effect_size_density <- function(test_type, effect_size = NULL, x, density, input = NULL, alternative = NULL, table_dim = NULL){
  branches <- .effect_size_internal_branches(
    x           = x,
    test_type   = test_type,
    effect_size = effect_size,
    input       = input,
    alternative = alternative,
    table_dim   = table_dim
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
  out
}

.effect_size_default_x_limit <- function(test_type, effect_size = NULL, alternative = NULL){
  effect_size <- .effect_size_normalize(test_type, effect_size)

  if(test_type %in% c("chi2_test", "f_test")){
    return(switch(
      effect_size,
      "contingency_coefficient" = c(0, 0.99),
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
