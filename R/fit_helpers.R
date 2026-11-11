# Helpers that fit one model and keep its treatment estimate. Every method
# returns a one-row data frame: the treatment coefficient (log HR, log OR, or
# mean difference), its SE, a status ("ok", "warning", "nonconvergence", or
# "error"), and the first warning or error message. A fit that fails still
# returns a row, with a missing estimate, so that the summaries can count it.

# `model_call` is a model call such as lm(y ~ x), passed as an argument. R does
# not run it here. It runs when capture_fit_with_status() first uses it, which
# is inside the handlers that record warnings and errors. For Cox models fitted
# with robust = TRUE, the reported SE is the robust SE.
fit_method <- function(model_call) {
  fit <- detect_nonconvergence(capture_fit_with_status(model_call))
  if (is.null(fit$model)) {
    return(failed_fit(fit$message))
  }
  coefficient <- tryCatch(
    treatment_coefficient(fit$model),
    error = function(e) list(estimate = NA_real_, se = NA_real_)
  )
  fit_result(coefficient$estimate, coefficient$se, fit$status, fit$message)
}

# The same, with a sandwich SE that is robust to unequal variances or, when
# `clusters` is given, to clustering within sites.
fit_with_sandwich_se <- function(model_call, clusters = NULL) {
  fit <- detect_nonconvergence(capture_fit_with_status(model_call))
  if (is.null(fit$model)) {
    return(failed_fit(fit$message))
  }
  fit_result(
    coef(fit$model)[["treatment"]], sandwich_se(fit$model, clusters),
    fit$status, fit$message
  )
}

# Continuous outcomes: estimate the treatment effect separately within each
# stratum (site or allocation ratio), then combine the estimates by
# inverse-variance weighting.
fit_pooled_within_strata <- function(data, strata) {
  stratum_fits <- lapply(split(data, strata), function(stratum) {
    treatment_coefficient(lm(continuous_outcome ~ treatment, data = stratum))
  })
  estimates <- vapply(stratum_fits, `[[`, numeric(1), "estimate")
  ses <- vapply(stratum_fits, `[[`, numeric(1), "se")
  weights <- 1 / ses^2
  fit_result(
    sum(weights * estimates) / sum(weights), sqrt(1 / sum(weights)),
    "ok", NA_character_
  )
}

# Evaluate a model call, recording warnings and errors instead of stopping.
# A fit that warns keeps its estimate and records the first warning.
capture_fit_with_status <- function(expr) {
  warnings_captured <- character()
  status <- "ok"
  message_text <- NA_character_

  model <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings_captured <<- c(warnings_captured, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      status <<- "error"
      message_text <<- conditionMessage(e)
      NULL
    }
  )

  if (!is.null(model) && length(warnings_captured) > 0) {
    status <- "warning"
    message_text <- warnings_captured[1]
  }
  list(model = model, status = status, message = message_text)
}

# Mark fits that the fitting software flags. lme4 lists convergence problems in
# its own messages, which also include singular fits (a variance estimated at
# its boundary), so status "nonconvergence" covers both for mixed models. glm
# sets a flag, and any warning it gave is kept as the message. coxph reports
# iteration problems by warning, so those fits have status "warning".
detect_nonconvergence <- function(fit) {
  model <- fit$model
  if (inherits(model, "merMod")) {
    lme4_messages <- unlist(model@optinfo$conv$lme4$messages)
    if (length(lme4_messages) > 0L) {
      fit$status <- "nonconvergence"
      fit$message <- lme4_messages[1]
    }
  } else if (inherits(model, "glm") && !model$converged) {
    fit$status <- "nonconvergence"
    if (is.na(fit$message)) fit$message <- "glm did not converge"
  }
  fit
}

# The treatment coefficient and its SE from the model summary. Cox models
# fitted with robust = TRUE report both SEs; use the robust one.
treatment_coefficient <- function(model) {
  row <- summary(model)$coefficients["treatment", ]
  if (!inherits(model, "coxph")) {
    return(list(estimate = row[["Estimate"]], se = row[["Std. Error"]]))
  }
  se_name <- if ("robust se" %in% names(row)) "robust se" else "se(coef)"
  list(estimate = row[["coef"]], se = row[[se_name]])
}

# HC0 sandwich SE of the treatment coefficient. With clusters, vcovCL() also
# applies its default small-sample factor G / (G - 1) for G sites. The robust
# SEs that coxph reports do not include that factor.
sandwich_se <- function(model, clusters = NULL) {
  covariance <- tryCatch(
    if (is.null(clusters)) {
      sandwich::vcovHC(model, type = "HC0")
    } else {
      sandwich::vcovCL(model, cluster = clusters, type = "HC0")
    },
    error = function(e) NULL
  )
  if (is.null(covariance)) {
    return(NA_real_)
  }
  sqrt(covariance["treatment", "treatment"])
}

# Nonfinite estimates or SEs are recorded as errors.
fit_result <- function(estimate, se, status, message) {
  if (!is.finite(estimate) || !is.finite(se) || se <= 0) {
    return(failed_fit("nonfinite coefficient or standard error"))
  }
  data.frame(estimate = estimate, se = se, status = status, message = message)
}

failed_fit <- function(message) {
  data.frame(
    estimate = NA_real_, se = NA_real_, status = "error", message = message
  )
}
