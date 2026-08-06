# ============================================================
#  gamlssCWM - Internal helper functions
#  These functions are not exported and are used internally
#  by the main fitting routine and its S3 methods.
# ============================================================

# Suppress R CMD check NOTE: 'weights_gamlss' is a data column
# created inside safe_gamlss(); it is not a global variable.
utils::globalVariables("weights_gamlss")

# ---- Family Function Finder ------------------------------------

#' Find a GAMLSS density function with fallback
#'
#' Searches for a GAMLSS density function (prefixed with "d") first within the
#' \code{gamlss.dist} namespace. If not found, it falls back to searching the
#' current search path and global environment (which is essential for custom or
#' truncated families generated at runtime, e.g., via \code{gamlss.tr}).
#' If the function cannot be located anywhere, it halts with an informative
#' error message.
#'
#' @param fam_name Character string specifying the GAMLSS family name
#'   (e.g., \code{"NO"}, \code{"GA"}, or a custom/truncated family name).
#'
#' @return The retrieved density function.
#'
#' @keywords internal
find_dfun <- function(fam_name) {
  fn <- paste0("d", fam_name)
  dFUN <- get0(fn, envir = asNamespace("gamlss.dist"), mode = "function")
  if (is.null(dFUN)) dFUN <- get0(fn, mode = "function")  # search path / global env
  if (is.null(dFUN)) {
    stop(sprintf(
      "Density function '%s' not found. Family '%s' is not supported (load the package defining it or use a gamlss.dist family).",
      fn, fam_name
    ), call. = FALSE)
  }
  dFUN
}

# ---- Initialisation ----------------------------------------

#' Initialize the Posterior Membership Matrix
#'
#' Produces an \eqn{n \times k} matrix of initial soft assignment weights
#' using *k*-means hard assignments, random soft assignments, or manual initialization.
#'
#' @param data A numeric `data.frame` (or matrix) with the observations.
#' @param mode Character string; either `"kmeans"` (default) for hard
#'   *k*-means initialization, `"rand.soft"` for random soft assignments,
#'   or `"manual"` to provide user-specified weights.
#' @param k Integer. Number of mixture components.
#' @param nstart Integer. Number of random starts (ignored for `mode = "manual"`).
#' @param start.z An \eqn{n \times k} numeric matrix of initial weights, required if `mode = "manual"`.
#' @param famX A `gamlss` family object for the marginal distribution of X.
#'   Used to check whether `kmeans` is applicable.
#' @param famY A `gamlss` family object for the conditional distribution of Y.
#'   Used to check whether `kmeans` is applicable.
#'
#' @return An \eqn{n \times k} numeric matrix of initial weights.
#'
#' @keywords internal
cwm_init <- function(data, mode, k, nstart = 1, start.z = NULL, famX = NULL, famY = NULL) {
  n <- nrow(data)

  if (!mode %in% c("kmeans", "rand.soft", "manual")) {
    stop("mode must be one of 'kmeans', 'rand.soft', or 'manual'")
  }

  if (mode == "kmeans") {
    if (!is.null(famX) && !is.null(famY) && famX$type == "Continuous" && famY$type == "Continuous") {
      res <- kmeans(data, centers = k, nstart = max(1, nstart))
      z <- as.matrix(model.matrix(~ factor(res$cluster) - 1))
      dimnames(z) <- NULL
    } else {
      warning("kmeans is only valid for continuous data -> switching to 'rand.soft'")
      mode <- "rand.soft"
    }
  }

  if (mode == "rand.soft") {
    best_score <- -Inf
    for (s in 1:max(1, nstart)) {
      temp <- matrix(runif(n * k), nrow = n, ncol = k)
      z_try <- temp / rowSums(temp)
      score <- sum(apply(z_try, 1, max))
      if (score > best_score) {
        best_score <- score
        z <- z_try
      }
    }
  } else if (mode == "manual") {
    if (is.null(start.z)) {
      stop("For 'manual' mode, 'start.z' must be provided")
    }
    if (!is.matrix(start.z)) {
      stop("'start.z' must be a matrix")
    }
    if (nrow(start.z) != n || ncol(start.z) != k) {
      stop(sprintf("'start.z' must be an %d x %d matrix", n, k))
    }
    if (any(start.z < 0)) {
      stop("'start.z' cannot contain negative values")
    }
    rs <- rowSums(start.z)
    if (any(!is.finite(rs))) {
      stop("Row sums of 'start.z' must be finite")
    }
    if (any(abs(rs - 1) > 1e-8)) {
      stop("Each row of 'start.z' must sum to 1 (within tolerance)")
    }
    z <- start.z
  }

  return(z)
}

# ---- Density evaluation (fitted values) --------------------

#' Evaluate a GAMLSS fitted density at observed values
#'
#' Given a fitted `gamlss` object and a numeric vector `x` with the same
#' length as the training data, evaluates the density of the fitted
#' distribution at `x` using the in-sample fitted parameter values.
#'
#' @param fit A fitted `gamlss` object, or `NULL`.  When `NULL`, a
#'   near-zero vector (`1e-10`) is returned so the EM algorithm can
#'   continue without crashing.
#' @param x Numeric vector at which to evaluate the density.
#'
#' @return Numeric vector of density values (same length as `x`).
#'
#' @keywords internal
dens_gamlss <- function(fit, x) {
  if (is.null(fit)) {
    return(rep(1e-10, length(x)))
  }

  fam_name <- fit$family[1]
  dFUN <- find_dfun(fam_name)

  fam_obj <- tryCatch(get(fam_name)(), error = function(e) NULL)
  if (!is.null(fam_obj$y.valid) && !is.null(x)) {
    if (!fam_obj$y.valid(x) && all(x >= 0, na.rm = TRUE)) {
      x <- pmax(x, 1e-5)
    }
  }

  args <- names(formals(dFUN))
  par_list <- list(x = x)

  if ("mu" %in% args && !is.null(fit$mu.fv)) {
    par_list$mu <- fit$mu.fv
  }
  if ("sigma" %in% args && !is.null(fit$sigma.fv)) {
    par_list$sigma <- fit$sigma.fv
  }
  if ("nu" %in% args && !is.null(fit$nu.fv)) {
    par_list$nu <- fit$nu.fv
  }
  if ("tau" %in% args && !is.null(fit$tau.fv)) {
    par_list$tau <- fit$tau.fv
  }


  suppressWarnings(do.call(dFUN, par_list))
}


# ---- Safe GAMLSS wrapper -----------------------------------

#' Fit a GAMLSS model, returning `NULL` on failure
#'
#' A lightweight wrapper around [gamlss::gamlss()] that catches errors
#' (via [try()]) and returns `NULL` instead of stopping the EM algorithm.
#' The training data are attached to the returned fit object so that
#' [stats::predict()] can use them when `newdata` is specified.
#'
#' @param formula A two-sided `formula` for the location parameter `mu`,
#'   e.g. `y ~ x`.
#' @param sigma_formula `formula` for the scale parameter `sigma`.
#'   Defaults to `~ 1` (intercept-only).
#' @param nu_formula `formula` for the third parameter `nu`.
#'   Defaults to `~ 1`.
#' @param tau_formula `formula` for the fourth parameter `tau`.
#'   Defaults to `~ 1`.
#' @param data A `data.frame` containing all variables in the formulae.
#' @param z Numeric vector of observation weights (soft assignments for
#'   the current EM component), floored at `1e-8`.
#' @param family A `gamlss` family object specifying the distribution.
#' @param ... Additional arguments passed to [gamlss::gamlss()].
#'
#' @return A fitted `gamlss` object (with `$data` attached), or `NULL`
#'   if fitting fails.
#'
#' @keywords internal
safe_gamlss <- function(formula,
                        sigma_formula = ~1,
                        nu_formula = ~1,
                        tau_formula = ~1,
                        data, z, family, prev_fit = NULL,
                        method = "RS", ...) {
  data$weights_gamlss <- pmax(as.numeric(z), 1e-8)

  # If the response violates family bounds, apply a minimal positive shift when safe.
  resp_var <- all.vars(formula)[1]
  if (length(resp_var) == 1 && resp_var %in% names(data)) {
    y <- data[[resp_var]]
    fam_obj <- if (is.character(family)) get(family)() else if (is.function(family)) family() else family
    if (!is.null(fam_obj$y.valid) && !fam_obj$y.valid(y)) {
      # Apply this fix only for non-negative data to avoid masking genuine data issues.
      if (all(!is.na(y)) && all(y >= 0) && fam_obj$y.valid(pmax(y, 1e-5))) {
        data[[resp_var]] <- pmax(y, 1e-5)
        warning(sprintf("Response '%s' contains values outside the support of family '%s'; shifted to 1e-5.",
                resp_var, fam_obj$family[1]), call. = FALSE)
      }
    }
  }

  start_args <- tryCatch({
    if (!is.null(prev_fit) && inherits(prev_fit, "gamlss")) {
      n_obs <- nrow(data)
      starts <- list()
      if (!is.null(prev_fit$mu.fv) && length(prev_fit$mu.fv) == n_obs) {
        starts$mu.start <- prev_fit$mu.fv
      }
      if (!is.null(prev_fit$sigma.fv) && length(prev_fit$sigma.fv) == n_obs) {
        starts$sigma.start <- prev_fit$sigma.fv
      }
      if (!is.null(prev_fit$nu.fv) && length(prev_fit$nu.fv) == n_obs) {
        starts$nu.start <- prev_fit$nu.fv
      }
      if (!is.null(prev_fit$tau.fv) && length(prev_fit$tau.fv) == n_obs) {
        starts$tau.start <- prev_fit$tau.fv
      }
      starts
    } else {
      list()
    }
  }, error = function(e) {
    list()
  })

  method_expr <- switch(method, 
                      RS    = quote(RS()), 
                      CG    = quote(CG()), 
                      mixed = quote(mixed()))

  gamlss_args <- c(list(
    formula = formula,
    sigma.formula = sigma_formula,
    nu.formula = nu_formula,
    tau.formula = tau_formula,
    data = data,
    weights = data$weights_gamlss,
    family = family,
    method = method_expr, 
    control = gamlss::gamlss.control(trace = FALSE, n.cyc = 500)
  ), start_args, list(...))

  fit <- try(
    suppressWarnings({
      do.call(gamlss::gamlss, gamlss_args)
    }),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(NULL)
  }

  fit$data <- data

  return(fit)
}


# ---- EM steps ----------------------------------------------

#' E-step of the EM algorithm
#'
#' Given the current M-step output (`M_out`) and the data, computes
#' posterior component-membership probabilities \eqn{z_{ij}} for each
#' observation \eqn{i} and component \eqn{j}.
#'
#' @param M_out List returned by `M_step()`, containing `$fits` (a list
#'   of per-component fit pairs) and `$pi` (mixing proportions).
#' @param data A `data.frame` with the observed variables.
#' @param k Integer. Number of mixture components.
#' @param formula The model `formula` used for fitting (used to extract
#'   variable names).
#'
#' @return An \eqn{n \times k} numeric matrix of posterior probabilities.
#'
#' @keywords internal
E_step <- function(M_out, data, k, formula) {
  n      <- nrow(data)
  pi     <- M_out$pi
  fits   <- M_out$fits
  y_name <- all.vars(formula)[1]
  x_name <- all.vars(formula)[2]

  log_dens <- matrix(0, nrow = n, ncol = k)
  for (j in 1:k) {
    densx <- pmax(dens_gamlss(fits[[j]]$fitX, data[[x_name]]), .Machine$double.xmin)
    densy <- pmax(dens_gamlss(fits[[j]]$fitY, data[[y_name]]), .Machine$double.xmin)
    log_dens[, j] <- log(pi[j]) + log(densx) + log(densy)
  }

  m     <- apply(log_dens, 1, max)
  lse   <- m + log(rowSums(exp(log_dens - m)))
  z_new <- exp(log_dens - lse)        
  return(z_new)
}


#' M-step of the EM algorithm
#'
#' Given the current soft assignments `z`, fits a GAMLSS model for `Y | X`
#' and a marginal GAMLSS model for `X` in each component using the supplied
#' formulas and family objects, then updates the mixing proportions.
#'
#' @param formula Two-sided `formula` for `Y`, e.g. `y ~ x`.
#' @param sigma_formula `formula` for the `sigma` parameter of `Y | X`.
#' @param nu_formula `formula` for the `nu` parameter of `Y | X`.
#' @param tau_formula `formula` for the `tau` parameter of `Y | X`.
#' @param z An \eqn{n \times k} matrix of current soft assignments.
#' @param famX A `gamlss` family object for the marginal distribution of `X`.
#' @param famY A `gamlss` family object for the conditional distribution of
#'   `Y | X`.
#' @param k Integer. Number of mixture components.
#' @param data A `data.frame` with the observed variables.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{`fits`}{A length-`k` list; each element is itself a list with
#'       `$fitX` (marginal fit for `X`) and `$fitY` (conditional fit for
#'       `Y | X`).}
#'     \item{`pi`}{Numeric vector of length `k` with updated mixing
#'       proportions.}
#'   }
#'
#' @keywords internal
M_step <- function(formula,
                   sigma_formula,
                   nu_formula,
                   tau_formula,
                   z, famX, famY, k, data, prev_fits = NULL, method = "RS") {
  fits <- vector(mode = "list", length = k)
  pi <- vector(mode = "numeric", length = k)

  for (j in 1:k) {
    prev_fitY <- NULL
    prev_fitX <- NULL
    if (!is.null(prev_fits) && length(prev_fits) >= j) {
      prev_fitY <- prev_fits[[j]]$fitY
      prev_fitX <- prev_fits[[j]]$fitX
    }

    # --- MODEL Y ---
    y_var <- data[[all.vars(formula)[1]]]

    fitY <- safe_gamlss(formula,
      sigma_formula,
      nu_formula,
      tau_formula,
      data, z[, j],
      family = famY,
      prev_fit = prev_fitY, 
      method = method 
    )

    # Fallback 1: cold-start
    if (is.null(fitY) && !is.null(prev_fitY)) {
      message(sprintf("  [EM] Component %d fitY failed with warm-start. Attempting cold-start...", j))
      fitY <- safe_gamlss(formula, sigma_formula, nu_formula, tau_formula,
                          data, z[, j], family = famY, prev_fit = NULL, method = method)
    }

    # Fallback 2: previous fit 
    if (!is.null(fitY) && !is.null(prev_fitY)) {
      if (keep_prev(fitY, prev_fitY, z[, j], y_var)) {
        message(sprintf("  [EM] Component %d fitY: log-likelihood decreased — rolling back to the previous fit", j))
        fitY <- prev_fitY
      }
    } else if (is.null(fitY) && !is.null(prev_fitY)) {
      fitY <- prev_fitY
    }

    if (is.null(fitY)) {
      stop(sprintf("Fit failed critically for component %d (Model Y). Stopping to avoid silent degeneration.", j))
    }

    # --- MODEL X ---
    x_name <- all.vars(formula)[2]
    formulaX <- as.formula(paste(x_name, "~ 1"))
    x_var <- data[[x_name]]

    fitX <- safe_gamlss(
      formula = formulaX,
      data = data,
      z = z[, j],
      family = famX,
      prev_fit = prev_fitX, 
      method = method
    )

    # Fallback 1: cold-start
    if (is.null(fitX) && !is.null(prev_fitX)) {
      fitX <- safe_gamlss(formulaX, data = data, z = z[, j], family = famX, prev_fit = NULL, method = method)
    }

    # Fallback 2: previous fit 
    if (!is.null(fitX) && !is.null(prev_fitX)) {
      if (keep_prev(fitX, prev_fitX, z[, j], x_var)) {
        warning(sprintf("  [EM] Component %d fitX: log-likelihood decreased — rolling back to the previous fit", j), call. = FALSE)
        fitX <- prev_fitX
      }
    } else if (is.null(fitX) && !is.null(prev_fitX)) {
      fitX <- prev_fitX
    }

    if (is.null(fitX)) {
      stop(sprintf("Fit failed critically for component %d (Model X).", j))
    }

    pi[j] <- sum(z[, j]) / nrow(data)
    # Check if pi collapse
    if (pi[j] < (k / nrow(data))) {
       warning(sprintf("Component %d has lost almost all of its mass (pi = %f). The algorithm is collapsing into a lower-dimensional mixture.", j, pi[j]), call. = FALSE)
       pi[j] <- 1e-8 
    }
    fits[[j]] <- list(fitX = fitX, fitY = fitY)
  }
  pi <- pi / sum(pi) # renorm pi
  return(list(fits = fits, pi = pi))
}


# ---- Log-likelihood ----------------------------------------

#' Compute the observed-data log-likelihood
#'
#' Evaluates the log-likelihood of the mixture model at the current
#' parameter estimates.
#'
#' @param M_out List returned by `M_step()`.
#' @param data A `data.frame` with the observed variables.
#' @param formula The model `formula` (used to extract variable names).
#' @param k Integer. Number of mixture components.
#'
#' @return A single numeric value: the observed-data log-likelihood.
#'
#' @keywords internal
compute_ll <- function(M_out, data, formula, k) {
  y_name <- all.vars(formula)[1]
  x_name <- all.vars(formula)[2]
  pi     <- M_out$pi
  fits   <- M_out$fits
  n      <- nrow(data)

  log_dens <- matrix(0, nrow = n, ncol = k)
  for (j in 1:k) {
    densx <- pmax(dens_gamlss(fits[[j]]$fitX, data[[x_name]]), .Machine$double.xmin)
    densy <- pmax(dens_gamlss(fits[[j]]$fitY, data[[y_name]]), .Machine$double.xmin)
    log_dens[, j] <- log(pi[j]) + log(densx) + log(densy)
  }

  m   <- apply(log_dens, 1, max)
  lse <- m + log(rowSums(exp(log_dens - m)))
  return(sum(lse))
}


# ---- Prediction helpers ------------------------------------

gamlss_predict_quiet <- function(object, what = "mu", ...) {
  withCallingHandlers(
    stats::predict(object, what = what, ...),
    warning = function(w) {
      msg <- conditionMessage(w)
      is_safe_refit_warning <- grepl(
        "discrepancy\\s+between the original and the re-fit\\s+used to achieve 'safe' predictions",
        msg,
        ignore.case = TRUE,
        perl = TRUE
      )
      if (is_safe_refit_warning) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Evaluate a GAMLSS density at new points
#'
#' Evaluates the density of a fitted `gamlss` distribution at `x_new`
#' using out-of-sample parameter predictions. Each distributional
#' parameter is predicted separately (avoiding the row-count mismatch
#' that can arise with `predictAll`).
#'
#' @param fit A fitted `gamlss` object, or `NULL`.  When `NULL`, a
#'   near-zero vector (`1e-10`) is returned.
#' @param x_new Numeric vector of points at which to evaluate the density.
#' @param newdata Optional `data.frame` of covariates for out-of-sample
#'   prediction.  If `NULL`, a one-column `data.frame` is constructed
#'   automatically from `x_new` using the response variable name of `fit`.
#'
#' @return Numeric vector of density values (same length as `x_new`).
#'
#' @importFrom stats formula setNames approx predict
#' @keywords internal
get_gamlss_density <- function(fit, x_new, newdata = NULL) {
  if (is.null(fit)) {
    return(rep(1e-10, length(x_new)))
  }

  fam_name <- fit$family[1]
  dFUN <- find_dfun(fam_name)
  d_args <- names(formals(dFUN))

  # If newdata is not supplied, construct it from x_new using the model response name.
  if (is.null(newdata)) {
    # Extract the response variable name (left-hand side of the fitted formula).
    resp_var <- deparse(formula(fit)[[2]])
    newdata <- setNames(data.frame(x_new), resp_var)
  }

  # predict() requires covariates only (exclude the response if present).
  pred_newdata <- as.data.frame(newdata)
  resp_var <- all.vars(formula(fit))[1]
  if (resp_var %in% names(pred_newdata)) {
    pred_newdata[[resp_var]] <- NULL
  }
  if (ncol(pred_newdata) == 0) {
    pred_newdata <- data.frame(.dummy = rep(0, length(x_new)))
  }

  # Predict parameters separately to avoid predictAll row-count mismatches
  # when newdata has a different size from the training data.
  par_list <- list(x = x_new)
  for (p in c("mu", "sigma", "nu", "tau")) {
    if (p %in% fit$parameters && p %in% d_args) {
      fv_name <- paste0(p, ".fv")
      fv <- fit[[fv_name]]

      pred_p <- try(
        gamlss_predict_quiet(
          fit,
          what = p,
          newdata = pred_newdata,
          data = fit$data,
          type = "response"
        ),
        silent = TRUE
      )

      if (inherits(pred_p, "try-error") || length(pred_p) != length(x_new)) {
        x_rhs_vars <- all.vars(formula(fit)[[3]])
        can_interpolate <-
          !is.null(fv) &&
          length(fv) == nrow(fit$data) &&
          length(x_rhs_vars) == 1 &&
          x_rhs_vars %in% names(pred_newdata) &&
          x_rhs_vars %in% names(fit$data)

        if (can_interpolate) {
          x_train <- fit$data[[x_rhs_vars]]
          x_out <- pred_newdata[[x_rhs_vars]]
          pred_p <- approx(x = x_train, y = fv, xout = x_out, rule = 2, ties = mean)$y
        } else if (!is.null(fv) && length(fv) > 0) {
          pred_p <- rep(fv[1], length(x_new))
        } else {
          pred_p <- rep(0, length(x_new))
        }
      }

      par_list[[p]] <- pred_p
    }
  }

  # Evaluate density values.
  dens <- suppressWarnings(do.call(dFUN, par_list))
  pmax(as.numeric(dens), 1e-10)
}

find_gamlss_obj <- function(name){
  obj <- get0(name, envir = globalenv(), inherits = TRUE, mode = "function")
  if(is.null(obj)){
    obj <- get0(name, envir = asNamespace("gamlss.dist"), mode = "function")
  }
  if(is.null(obj)){
    stop(sprintf("Unable to find ‘%s’. Please verify that the family/distribution is defined and loaded (e.g., via gamlss.tr::gen.trun() or library()).",
      name
    ))
  }
  obj
}

#' Compare Two GAMLSS Fits Using a Weighted Log-Likelihood
#'
#' Determines whether the previous fit should be retained based on the
#' weighted log-likelihood of the fitted densities.
#'
#' @param new_fit The newly obtained fitted model.
#' @param prev_fit The previously retained fitted model.
#' @param w A numeric vector of weights used to compute the weighted
#'   log-likelihood.
#' @param y The response values used to evaluate the fitted densities.
#'
#' @return A logical value indicating whether the previous fit should be
#'   retained. Returns `TRUE` if `new_fit` is `NULL`, or if the weighted
#'   log-likelihood of `prev_fit` is finite and greater than that of
#'   `new_fit`; otherwise returns `FALSE`.
#'
#' @keywords internal
keep_prev <- function(new_fit, prev_fit, w, y) {
  if (is.null(new_fit)) return(TRUE)
  if (is.null(prev_fit)) return(FALSE)

  ll_new <- sum(
    w * log(pmax(dens_gamlss(new_fit, y), .Machine$double.xmin))
  )
  ll_prev <- sum(
    w * log(pmax(dens_gamlss(prev_fit, y), .Machine$double.xmin))
  )

  return(is.finite(ll_prev) && ll_prev > ll_new)
}