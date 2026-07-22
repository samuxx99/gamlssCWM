# ============================================================
#  gamlssCWM — Exported functions and S3 methods
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b


# ---- Main fitting function ---------------------------------

#' Fit a Cluster-Weighted GAMLSS Model via EM
#'
#' Estimates a *k*-component bivariate Cluster-Weighted Model (CWM) in
#' which each component is characterized by:
#' \itemize{
#'   \item a marginal GAMLSS distribution for the covariate \eqn{X}, and
#'   \item a conditional GAMLSS distribution for the response \eqn{Y | X}.
#' }
#' Parameters are estimated via the Expectation-Maximization (EM)
#' algorithm.  The E-step updates the posterior component-membership
#' probabilities; the M-step re-fits a pair of weighted GAMLSS models
#' for each component.
#'
#' @param formula A two-sided `formula` specifying the conditional model
#'   for the location parameter (`mu`) of \eqn{Y | X}, e.g. `y ~ x`.
#' @param sigma_formula A one-sided `formula` for the scale parameter
#'   `sigma` of the \eqn{Y | X} model.  Defaults to `~ 1`
#'   (intercept-only / constant sigma).
#' @param nu_formula A one-sided `formula` for the third distributional
#'   parameter `nu`.  Defaults to `~ 1`.
#' @param tau_formula A one-sided `formula` for the fourth distributional
#'   parameter `tau`.  Defaults to `~ 1`.
#' @param data A `data.frame` (or coercible object) containing all
#'   variables referenced in the formulae.
#' @param famX A `gamlss` family object specifying the marginal
#'   distribution of \eqn{X} (e.g. [gamlss.dist::NO()]).
#' @param famY A `gamlss` family object specifying the conditional
#'   distribution of \eqn{Y | X} (e.g. [gamlss.dist::NO()]).
#' @param init Character string indicating the initialization strategy:
#'   `"kmeans"` for hard *k*-means assignments, or
#'   `"rand.soft"` (default) for random row-normalised weights.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance on the absolute change in
#'   log-likelihood between consecutive iterations.  Defaults to `1e-3`.
#' @param k Integer. Number of mixture components.
#' @param nstart Integer. Number of random starts for the initialization.
#' @param start.z Optional \eqn{n \times k} numeric matrix of initial
#'   component-membership probabilities.
#' @param seed Optional integer. If provided, the RNG state is set at the
#'   start of the fit and restored on exit for fully reproducible results.
#' @param verbose Logical. If `TRUE` (default), a progress bar is shown during EM iterations.
#'
#' @return An object of class `"gamlssCWM"`, which is a named list with
#'   the following elements:
#'   \describe{
#'     \item{`fits`}{A length-`k` list of per-component fit pairs, each
#'       containing `$fitX` and `$fitY` (`gamlss` objects or `NULL` if
#'       fitting failed).}
#'     \item{`pi`}{Numeric vector of length `k` with the estimated
#'       mixing proportions.}
#'     \item{`posteriors`}{An \eqn{n \times k} matrix of final posterior
#'       component-membership probabilities.}
#'     \item{`loglik`}{Scalar. Observed-data log-likelihood at convergence.}
#'     \item{`ll_history`}{Numeric vector of log-likelihood values across
#'       all EM iterations up to convergence.}
#'     \item {`converged`}{Logical; \code{TRUE} if the EM algorithm reached 
#'       the convergence criterion.} 
#'     \item{`call`}{The `formula` argument (stored for use by S3
#'       methods).}
#'   }
#'
#' @seealso `predict.gamlssCWM`, `plot.gamlssCWM`, `summary.gamlssCWM`
#'
#' @examples
#' library(gamlss.dist)
#'
#' set.seed(42)
#' n <- 200
#' x <- c(rnorm(n / 2, -2), rnorm(n / 2, 2))
#' y <- c(
#'   1 + 0.8 * x[1:(n / 2)] + rnorm(n / 2, sd = 0.5),
#'   -1 - 0.8 * x[(n / 2 + 1):n] + rnorm(n / 2, sd = 0.5)
#' )
#' df <- data.frame(x = x, y = y)
#'
#' fit <- gamlssCWM(
#'   formula  = y ~ x,
#'   data     = df,
#'   famX     = NO(),
#'   famY     = NO(),
#'   k        = 2,
#'   max_iter = 50
#' )
#'
#' summary(fit)
#'
#' @importFrom stats kmeans model.matrix runif as.formula
#' @importFrom stats logLik nobs AIC BIC
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom gamlss.dist NO
#'
#' @export
gamlssCWM <- function(formula,
                      sigma_formula = ~1,
                      nu_formula = ~1,
                      tau_formula = ~1,
                      data,
                      famX,
                      famY,
                      init = "rand.soft",
                      max_iter,
                      tol = 1e-3,
                      k,
                      nstart = 1,
                      start.z = NULL, 
                      seed = NULL,
                      verbose = TRUE) {
  # --- input validation ---
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("'formula' must be a two-sided formula, e.g. y ~ x.")
  }
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data.frame or matrix.")
  }
  data <- as.data.frame(data)
  vars <- all.vars(formula)
  if (!all(vars %in% names(data))) {
    stop(
      "Variables in 'formula' not found in 'data': ",
      paste(setdiff(vars, names(data)), collapse = ", ")
    )
  }
  if (!is.numeric(k) || length(k) != 1L || k < 1L || k != round(k)) {
    stop("'k' must be a positive integer (k >= 1).")
  }
  if (nrow(data) <= k) {
    stop("'data' must have more rows than 'k'.")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
    max_iter < 1L || max_iter != round(max_iter)) {
    stop("'max_iter' must be a positive integer.")
  }
  if (!is.numeric(nstart) || length(nstart) != 1L ||
    nstart < 1L || nstart != round(nstart)) {
    stop("'nstart' must be a positive integer.")
  }
  if (!inherits(famX, "gamlss.family") || !inherits(famY, "gamlss.family")) {
    stop("'famX' and 'famY' must be gamlss.family objects, e.g. NO().")
  }
  init <- match.arg(init, c("kmeans", "rand.soft", "manual"))
  # --- end validation ---

  formulas_health <- list(
    sigma = sigma_formula,
    nu = nu_formula,
    tau = tau_formula
  )

  for (param in names(formulas_health)) {
    has_formula <- deparse(formulas_health[[param]]) != "~1"
    param_exists <- param %in% names(famY$parameters)

    if (has_formula && !param_exists) {
      warning(sprintf(
        "Parameter '%s' is not in %s distribution and will be ignored",
        param,
        famY$family[1]
      ), call. = FALSE) # Nota: aggiunto call. = FALSE per un output più pulito
    }
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be a single finite numeric value.")
    }
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      stats::runif(1)
    }
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    set.seed(as.integer(seed))
  }
  loglik <- vector("numeric", length = max_iter)
  M_out <- vector("list", length = max_iter)

  if(k > 1){
    z <- matrix(NA, ncol = k, nrow = nrow(data))

    z <- cwm_init(
      data = data,
      mode = init,
      k = k,
      nstart = nstart,
      famX = famX,
      famY = famY,
      start.z = start.z
    )
  } else{
    z <- matrix(1, ncol = k, nrow = nrow(data))
  }

  converged <- FALSE
  if(verbose) pb <- txtProgressBar(min = 0, max = max_iter, style = 3)

  for (r in 1:max_iter) {
    if (r == 1) {
      M_out[[r]] <- M_step(
        formula = formula,
        sigma_formula = sigma_formula,
        nu_formula = nu_formula,
        tau_formula = tau_formula,
        z = z, famX = famX, famY = famY,
        k = k, data = data,
        prev_fits = NULL
      )

      loglik[r] <- compute_ll(M_out[[r]], data, formula, k)
    } else {
      z <- E_step(M_out[[r - 1]], data = data, k = k, formula = formula)

      M_out[[r]] <- M_step(
        formula = formula,
        sigma_formula = sigma_formula,
        nu_formula = nu_formula,
        tau_formula = tau_formula,
        z = z, famX = famX, famY = famY,
        k = k, data = data,
        prev_fits = M_out[[r - 1]]$fits
      )

      loglik[r] <- compute_ll(M_out[[r]], data, formula, k)
      if (r >= 3) M_out[r - 2] <- list(NULL)
    }

    converged_now <- if (r >= 3) {
      denom <- loglik[r - 1] - loglik[r - 2]
      if (abs(denom) < .Machine$double.eps^0.5) {
        abs(loglik[r] - loglik[r - 1]) < tol # fallback: flat LL
      } else {
        a_r <- loglik[r - 1] + (loglik[r] - loglik[r - 1]) / (1 - (loglik[r] - loglik[r - 1]) / denom)
        abs(a_r - loglik[r]) < tol
      }
    } else {
      FALSE
    }

    if (converged_now) {
      if (verbose) {
        setTxtProgressBar(pb, max_iter)
        close(pb)
      }


      if (verbose) {
        message(sprintf("EM converges at iteration %d", r))
      }

      converged <- TRUE
      break
    }
    if(verbose) setTxtProgressBar(pb, r)
  }

  last_M <- M_out[[r]]
  z <- E_step(last_M, data = data, k = k, formula = formula)

  res <- list(
    fits = last_M$fits,
    pi = last_M$pi,
    posteriors = z,
    loglik = loglik[r],
    ll_history = loglik[1:r],
    converged = converged,
    call = formula
  )

  class(res) <- "gamlssCWM"

  return(res)
}


# ---- S3 methods --------------------------------------------

#' Print a Cluster-Weighted GAMLSS Model
#'
#' Prints a concise summary of a fitted `gamlssCWM` object, including
#' the model formula, number of components, sample size, final
#' log-likelihood, number of EM iterations, estimated mixing
#' proportions, and per-component distribution families.
#'
#' @param x An object of class `"gamlssCWM"` returned by
#'   `gamlssCWM()`.
#' @param digits Integer indicating the number of significant digits to
#'   use for numeric values.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \dontshow{
#' library(gamlss.dist)
#' set.seed(42)
#' n <- 200
#' x <- c(rnorm(n / 2, -2), rnorm(n / 2, 2))
#' y <- c(
#'   1 + 0.8 * x[1:(n / 2)] + rnorm(n / 2, sd = 0.5),
#'   -1 - 0.8 * x[(n / 2 + 1):n] + rnorm(n / 2, sd = 0.5)
#' )
#' df <- data.frame(x = x, y = y)
#' fit <- gamlssCWM(
#'   formula  = y ~ x,
#'   data     = df,
#'   famX     = NO(),
#'   famY     = NO(),
#'   k        = 2,
#'   max_iter = 50
#' )
#' }
#' print(fit)
#'
#' @method print gamlssCWM
#' @export
print.gamlssCWM <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  object <- x

  k <- length(object$fits)
  n_obs <- if (!is.null(object$posteriors)) nrow(object$posteriors) else NA_integer_
  n_iter <- if (!is.null(object$ll_history)) length(object$ll_history) else NA_integer_

  cat("Family:  Cluster-Weighted GAMLSS Mixture\n")
  cat("Formula: ", paste(deparse(object$call), collapse = " "), "\n", sep = "")

  cat("Number of components: ", k, "\n", sep = "")
  if (!is.na(n_obs)) {
    cat("Number of observations: ", n_obs, "\n", sep = "")
  }
  if (!is.null(object$loglik)) {
    cat("Final log-likelihood: ", format(signif(object$loglik, digits = digits)), "\n", sep = "")
  }
  if (!is.na(n_iter)) {
    cat("EM iterations: ", n_iter, "\n", sep = "")
  }
  cat("Converged: ", object$converged, "\n", sep = "")

  if (!is.null(object$pi)) {
    pi_str <- paste(formatC(object$pi, format = "f", digits = digits), collapse = ", ")
    cat("Mixing proportions: ", pi_str, "\n", sep = "")
  }

  cat("Component families:\n")
  for (j in seq_len(k)) {
    fitX <- object$fits[[j]]$fitX
    fitY <- object$fits[[j]]$fitY

    famX <- if (is.null(fitX)) "failed" else fitX$family[1]
    famY <- if (is.null(fitY)) "failed" else fitY$family[1]

    cat(" ", j, ": X ~ ", famX, " | Y|X ~ ", famY, "\n", sep = "")
  }

  invisible(x)
}

#' Predict from a Cluster-Weighted GAMLSS Model
#'
#' Generates predictions from a fitted `gamlssCWM` object for new
#' observations.  Three types of output are supported.
#'
#' @param object An object of class `"gamlssCWM"` returned by
#'   `gamlssCWM()`.
#' @param newdata A `data.frame` containing the predictor and response
#'   variables for the new observations.
#' @param type Character string specifying the type of prediction:
#'   \describe{
#'     \item{`"response"`}{(default) The posterior-weighted expected value
#'       of \eqn{Y}, i.e. \eqn{\sum_j \hat{p}_{ij} \hat{\mu}_{ij}}.}
#'     \item{`"posterior"`}{An \eqn{n \times k} matrix of posterior
#'       component-membership probabilities based on the marginal density
#'       of \eqn{X} only.}
#'     \item{`"joint"`}{An \eqn{n \times k} matrix of weighted joint
#'       densities \eqn{\pi_j f_j(x) f_j(y | x)}.}
#'   }
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return Depending on `type`:
#'   \itemize{
#'     \item `"response"` — numeric vector of length `nrow(newdata)`.
#'     \item `"posterior"` — numeric matrix of dimension
#'       \eqn{n \times k}.
#'     \item `"joint"` — numeric matrix of dimension \eqn{n \times k}.
#'   }
#'
#' @examples
#' \dontshow{
#' library(gamlss.dist)
#' set.seed(42)
#' n <- 200
#' x <- c(rnorm(n / 2, -2), rnorm(n / 2, 2))
#' y <- c(
#'   1 + 0.8 * x[1:(n / 2)] + rnorm(n / 2, sd = 0.5),
#'   -1 - 0.8 * x[(n / 2 + 1):n] + rnorm(n / 2, sd = 0.5)
#' )
#' df <- data.frame(x = x, y = y)
#' fit <- gamlssCWM(
#'   formula  = y ~ x,
#'   data     = df,
#'   famX     = NO(),
#'   famY     = NO(),
#'   k        = 2,
#'   max_iter = 50
#' )
#' }
#' y_hat <- predict(fit, newdata = df, type = "response")
#' post <- predict(fit, newdata = df, type = "posterior")
#'
#' @importFrom stats predict
#'
#' @method predict gamlssCWM
#' @export
predict.gamlssCWM <- function(object, newdata, type = c("response", "posterior", "joint"), ...) {
  type <- match.arg(type)
  k <- length(object$fits)
  n_new <- nrow(newdata)

  x_name <- all.vars(object$call)[2]
  y_name <- all.vars(object$call)[1]

  dens_x <- matrix(0, nrow = n_new, ncol = k)
  dens_y <- matrix(0, nrow = n_new, ncol = k)
  mu_y <- NULL
  needs_mu <- identical(type, "response")
  if (needs_mu) {
    mu_y <- matrix(0, nrow = n_new, ncol = k)
  }

  for (j in 1:k) {
    dens_x[, j] <- get_gamlss_density(object$fits[[j]]$fitX,
      x_new = newdata[[x_name]]
    )

    if (type == "joint") {
      dens_y[, j] <- get_gamlss_density(object$fits[[j]]$fitY,
        x_new   = newdata[[y_name]],
        newdata = newdata
      )
    }

    if (needs_mu) {
      mu_y[, j] <- gamlss_predict_quiet(object$fits[[j]]$fitY,
        what = "mu", newdata = newdata,
        type = "response",
        data = object$fits[[j]]$fitY$data
      )
    }
  }

  joint_components <- sweep(dens_x * dens_y, 2, object$pi, "*")

  if (type == "joint") {
    return(joint_components)
  }

  log_dx <- sweep(log(pmax(dens_x, .Machine$double.xmin)), 2, log(object$pi), "+")
  m_dx <- apply(log_dx, 1, max)
  post_x <- exp(log_dx - (m_dx + log(rowSums(exp(log_dx - m_dx)))))
  if (type == "posterior") {
    return(post_x)
  }

  return(rowSums(post_x * mu_y))
}


#' Plot a Cluster-Weighted GAMLSS Model
#'
#' Produces a scatter plot of \eqn{Y} against \eqn{X} with observations
#' coloured by their most likely component (MAP assignment).  Component-
#' specific regression curves (fitted mu values) are overlaid.
#' Optionally, density contours of the full mixture and of each
#' individual component can be added.
#'
#' @param x An object of class `"gamlssCWM"` returned by
#'   `gamlssCWM()`.
#' @param ... Currently unused; included for S3 generic compatibility.
#' @param dataX Numeric vector of \eqn{X} values to plot (typically the
#'   training data).
#' @param dataY Numeric vector of \eqn{Y} values to plot.
#' @param col Character or integer vector of length \eqn{k} specifying
#'   the colours for the components. Defaults to `1:k`.
#' @param pch Point character; passed to [graphics::plot()]. Defaults to `1`.
#' @param lwd Line width; passed to [graphics::lines()] and [graphics::contour()]. Defaults to `2`.
#' @param contours Logical.  If `TRUE`, iso-density contours of the
#'   mixture (black solid) and each component (coloured dashed) are
#'   overlaid.  Defaults to `FALSE`. Contours are computed only when all
#'   component families are continuous; if any family is not continuous, a
#'   warning is issued and contours are skipped.
#' @param drawlabels Logical. If `TRUE`, contours are labeled with their
#'   corresponding density values. Defaults to `FALSE`.
#' @param main Character string for the plot title.  Defaults to `""`.
#' @param ylab Label for the \eqn{y}-axis.  Defaults to `"Y"`.
#' @param xlab Label for the \eqn{x}-axis.  Defaults to `"X"`.
#'
#' @return Invisibly returns `NULL`.  Called for its side-effect of
#'   producing a plot.
#'
#' @examples
#' \dontshow{
#' library(gamlss.dist)
#' set.seed(42)
#' n <- 200
#' x <- c(rnorm(n / 2, -2), rnorm(n / 2, 2))
#' y <- c(
#'   1 + 0.8 * x[1:(n / 2)] + rnorm(n / 2, sd = 0.5),
#'   -1 - 0.8 * x[(n / 2 + 1):n] + rnorm(n / 2, sd = 0.5)
#' )
#' df <- data.frame(x = x, y = y)
#' fit <- gamlssCWM(
#'   formula  = y ~ x,
#'   data     = df,
#'   famX     = NO(),
#'   famY     = NO(),
#'   k        = 2,
#'   max_iter = 50
#' )
#' }
#' plot(fit,
#'   dataX = df$x, dataY = df$y,
#'   col = c("steelblue", "tomato"), pch = 16,
#'   main = "Cluster-Weighted GAMLSS fit"
#' )
#'
#' @importFrom graphics plot lines contour
#'
#' @method plot gamlssCWM
#' @export
plot.gamlssCWM <- function(x, ..., dataX, dataY,
                           col = NULL, pch = 1, lwd = 2, contours = FALSE, drawlabels = FALSE,
                           main = "", ylab = NULL, xlab = NULL) {
  object <- x
  length(object$fits) -> k

  if (is.null(col)) {
    col <- seq_len(k)
  }

  cluster_ids <- apply(object$posteriors, 1, which.max)

  y_name <- all.vars(object$call)[1]
  x_name <- all.vars(object$call)[2]

  if (is.data.frame(dataX)) {
    data_df <- dataX
    if (!x_name %in% names(data_df) || !y_name %in% names(data_df)) {
      stop("dataX does not contain the variables in the fitted formula.")
    }
    dataX <- data_df[[x_name]]
    if (missing(dataY) || is.null(dataY)) {
      dataY <- data_df[[y_name]]
    } else if (is.character(dataY) && length(dataY) == 1L && dataY %in% names(data_df)) {
      dataY <- data_df[[dataY]]
    }
  }

  if (!is.numeric(dataX) || !is.numeric(dataY)) {
    stop("dataX and dataY must be numeric vectors (or dataX a data.frame with numeric columns).")
  }
  if (length(dataX) != length(dataY)) {
    stop("dataX and dataY must have the same length.")
  }

  if (is.null(xlab)) xlab <- x_name
  if (is.null(ylab)) ylab <- y_name

  plot(dataY ~ dataX,
    pch = pch,
    col = col[cluster_ids], main = main, ylab = ylab, xlab = xlab, ...
  )

  data <- data.frame(x = dataX, y = dataY, lab = cluster_ids)

  ## Component-specific regression lines
  for (j in 1:k) {
    if (is.null(object$fits[[j]]$fitY)) next
    if(!any(data$lab == j)) next
    x_grid <- seq(min(data$x[data$lab == j]), max(data$x[data$lab == j]), length.out = 200)
    nd <- setNames(data.frame(x_grid), x_name)
    mu_j <- try(
      gamlss_predict_quiet(object$fits[[j]]$fitY,
        what    = "mu",
        newdata = nd,
        type    = "response",
        data    = object$fits[[j]]$fitY$data
      ),
      silent = TRUE
    )
    if (inherits(mu_j, "try-error")) next
    lines(x_grid, mu_j, col = col[j], lwd = lwd)
  }

  ## Density contours
  if (contours) {

    family_type <- function(fit_obj) {
      if (is.null(fit_obj)) {
        return("continuous")
      }
      tolower(as.character(fit_obj$type %||% "continuous"))
    }

    has_non_continuous_family <- any(vapply(object$fits, function(comp_fit) {
      any(c(family_type(comp_fit$fitX), family_type(comp_fit$fitY)) != "continuous")
    }, logical(1)))

    if (has_non_continuous_family) {
      warning("Contours are computable only for continuous families.", call. = FALSE)
      return(invisible(NULL))
    }

    x_grid <- seq(min(data$x) - 0.05 * diff(range(data$x)), max(data$x) + 0.05 * diff(range(data$x)), length.out = 101)
    y_grid <- seq(min(data$y) - 0.05 * diff(range(data$y)), max(data$y) + 0.05 * diff(range(data$y)), length.out = 101)

    grid <- expand.grid(x = x_grid, y = y_grid) # Use names consistent with the fitted formula.
    colnames(grid) <- c(x_name, y_name)

    cont_mixt <- predict.gamlssCWM(object = object, newdata = grid, type = "joint")

    # cont_mixt is an n_grid x k matrix; sum across components for mixture density.
    mixt_dens <- rowSums(cont_mixt)
    # expand.grid varies x first, so this matrix layout preserves z(x_i, y_j).
    z_mat <- matrix(mixt_dens, nrow = length(x_grid), ncol = length(y_grid))

    contour(x_grid, y_grid, z_mat,
      add = TRUE,
      col = "black", lwd = lwd, drawlabels = drawlabels
    )

    for (j in 1:k) {
      # Marginal density of X; only x-values are required.
      dens_x_grid <- get_gamlss_density(object$fits[[j]]$fitX, grid[[x_name]])

      # Conditional density Y|X: pass the entire grid so that predictAll
      # can use X to predict mu and sigma, then evaluate the density at Y
      dens_y_grid <- get_gamlss_density(object$fits[[j]]$fitY,
        x_new = grid[[y_name]],
        newdata = grid
      )

      # Component-wise joint density: pi_j * f_j(x) * f_j(y | x).
      z_j <- dens_x_grid * dens_y_grid * object$pi[j]
      z_mat_j <- matrix(z_j, nrow = length(x_grid), ncol = length(y_grid))

      contour(x_grid, y_grid, z_mat_j, add = TRUE, col = col[j], lty = 2, drawlabels = drawlabels)
    }
  }
}

#' Perspective Plot for a Cluster-Weighted GAMLSS Model
#'
#' Creates a 3D perspective plot (`persp`) of the joint density
#' estimated by a fitted \code{gamlssCWM} object.
#'
#' @param x An object of class \code{gamlssCWM} returned by \code{gamlssCWM()}.
#' @param dataX Numeric vector with the original \eqn{X} data (used to set the grid range).
#' @param dataY Numeric vector with the original \eqn{Y} data (used to set the grid range).
#' @param theta Numeric specifying the azimuthal viewing angle. Default is 30.
#' @param phi Numeric specifying the colatitude viewing angle. Default is 30.
#' @param col Character string specifying the color of the surface. Default is \code{"white"}.
#' @param border Character string specifying the color of the lines on the surface. Default is \code{"black"}.
#' @param shade Numeric value specifying the shading of the surface. Default is \code{0.5}.
#' @param xlab Character string for the x-axis label.
#' @param ylab Character string for the y-axis label.
#' @param zlab Character string for the z-axis label.
#' @param main Character string specifying the main title.
#' @param ... Additional arguments passed to \code{\link[graphics]{persp}}.
#'
#' @seealso \code{plot.gamlssCWM}
#'
#' @importFrom graphics persp
#'
#' @method persp gamlssCWM
#' @export
persp.gamlssCWM <- function(x, dataX, dataY,
                            theta = 30, phi = 30,
                            col = "white", border = "black", shade = 0.5,
                            xlab = NULL, ylab = NULL, zlab = "Density", main = "", ...) {
  object <- x

  family_type <- function(fit_obj) {
    if (is.null(fit_obj)) {
      return("continuous")
    }
    tolower(as.character(fit_obj$type %||% "continuous"))
  }

  has_non_continuous_family <- any(vapply(object$fits, function(comp_fit) {
    any(c(family_type(comp_fit$fitX), family_type(comp_fit$fitY)) != "continuous")
  }, logical(1)))

  if (has_non_continuous_family) {
    stop("Perspective plots are computable only for continuous families.")
  }

  y_name <- all.vars(object$call)[1]
  x_name <- all.vars(object$call)[2]

  if (is.null(xlab)) xlab <- x_name
  if (is.null(ylab)) ylab <- y_name

  # Set up the grid
  x_grid <- seq(min(dataX) - 0.05 * diff(range(dataX)), max(dataX) + 0.05 * diff(range(dataX)), length.out = 101)
  y_grid <- seq(min(dataY) - 0.05 * diff(range(dataY)), max(dataY) + 0.05 * diff(range(dataY)), length.out = 101)

  grid <- expand.grid(x = x_grid, y = y_grid)
  colnames(grid) <- c(x_name, y_name)

  # Compute joint density using the existing predict method
  cont_mixt <- predict.gamlssCWM(object = object, newdata = grid, type = "joint")
  mixt_dens <- rowSums(cont_mixt)

  # Create the Z matrix
  z_mat <- matrix(mixt_dens, nrow = length(x_grid), ncol = length(y_grid))

  # Draw perspective plot
  persp(x_grid, y_grid, z_mat,
    theta = theta, phi = phi, col = col, border = border, shade = shade,
    xlab = xlab, ylab = ylab, zlab = zlab, main = main, ...
  )
}


#' Summarize a Cluster-Weighted GAMLSS Model
#'
#' Prints a formatted summary of a fitted `gamlssCWM` model, including
#' mixing proportions, the final log-likelihood and, for each component,
#' the estimated coefficients (with standard errors and significance
#' stars) of both the marginal \eqn{X} model and the conditional
#' \eqn{Y | X} model.
#'
#' @param object An object of class `"gamlssCWM"` returned by
#'   `gamlssCWM()`.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return Invisibly returns `object`.
#'
#' @examples
#' \dontshow{
#' library(gamlss.dist)
#' set.seed(42)
#' n <- 200
#' x <- c(rnorm(n / 2, -2), rnorm(n / 2, 2))
#' y <- c(
#'   1 + 0.8 * x[1:(n / 2)] + rnorm(n / 2, sd = 0.5),
#'   -1 - 0.8 * x[(n / 2 + 1):n] + rnorm(n / 2, sd = 0.5)
#' )
#' df <- data.frame(x = x, y = y)
#' fit <- gamlssCWM(
#'   formula  = y ~ x,
#'   data     = df,
#'   famX     = NO(),
#'   famY     = NO(),
#'   k        = 2,
#'   max_iter = 50
#' )
#' }
#' summary(fit)
#'
#' @importFrom stats printCoefmat
#' @importFrom utils capture.output
#'
#' @method summary gamlssCWM
#' @export
summary.gamlssCWM <- function(object, ...) {
  k <- length(object$fits)
  printed_any <- FALSE

  cat("==================================================================\n")
  cat("        Cluster-Weighted GAMLSS Model (EM Algorithm)              \n")
  cat("==================================================================\n")
  cat("\n")

  cat("Number of components : ", k, "\n", sep = "")
  cat("Final Log-Likelihood : ", round(object$loglik, 4), "\n", sep = "")

  pi_str <- paste(sprintf("%.4f", object$pi), collapse = ", ")
  cat("Mixing proportions   : ", pi_str, "\n", sep = "")
  cat("\n")

  for (j in 1:k) {
    cat("------------------------------------------------------------------\n")
    cat(sprintf(" COMPONENT %d  (pi = %.4f)\n", j, object$pi[j]))

    # --- MODEL X ---
    fitX <- object$fits[[j]]$fitX
    if (is.null(fitX)) {
      cat("-> MODEL X | fit failed\n\n")
    } else {
      cat("-> MODEL X | Distribution: ", fitX$family[1], "\n", sep = "")

      invisible(capture.output(sumX <- suppressWarnings(summary(fitX, save = TRUE))))

      for (p in fitX$parameters) {
        cat(sprintf("   * %s (link: %s)\n", p, fitX[[paste0(p, ".link")]]))
        printCoefmat(sumX[[paste0(p, ".coef.table")]], signif.stars = TRUE, signif.legend = FALSE)
        printed_any <- TRUE
        cat("\n")
      }
    }

    # --- MODEL Y ---
    fitY <- object$fits[[j]]$fitY
    if (is.null(fitY)) {
      cat("-> MODEL Y | fit failed\n\n")
    } else {
      cat("-> MODEL Y | Distribution: ", fitY$family[1], "\n", sep = "")

      invisible(capture.output(sumY <- suppressWarnings(summary(fitY, save = TRUE))))

      for (p in fitY$parameters) {
        cat(sprintf("   * %s (link: %s)\n", p, fitY[[paste0(p, ".link")]]))
        printCoefmat(sumY[[paste0(p, ".coef.table")]], signif.stars = TRUE, signif.legend = FALSE)
        printed_any <- TRUE
        cat("\n")
      }
    }
  }

  cat("==================================================================\n")
  if (printed_any) {
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }

  invisible(object)
}

#' @method logLik gamlssCWM
#' @export
logLik.gamlssCWM <- function(object, ...) {
  # Effective degrees of freedom: component df.fit totals plus (k - 1) mixing parameters.
  k <- length(object$fits)
  df_components <- sum(vapply(object$fits, function(comp) {
    (comp$fitX$df.fit %||% 0) + (comp$fitY$df.fit %||% 0)
  }, numeric(1)))
  df_total <- df_components + (k - 1L)
  structure(object$loglik, df = df_total, class = "logLik")
}

#' @method AIC gamlssCWM
#' @export
AIC.gamlssCWM <- function(object, ..., k = 2) {
  ll <- logLik(object)
  m <- attr(ll, "df")
  -2 * as.numeric(ll) + k * m
}

#' @method BIC gamlssCWM
#' @export
BIC.gamlssCWM <- function(object, ...) {
  ll <- logLik(object)
  m <- attr(ll, "df")
  n <- nobs(object)
  -2 * as.numeric(ll) + log(n) * m
}

#' @method nobs gamlssCWM
#' @export
nobs.gamlssCWM <- function(object, ...) {
  nrow(object$posteriors)
}

#' Information Criteria for Cluster-Weighted GAMLSS Models
#'
#' Computes several information criteria for model selection, including AIC,
#' BIC, and the Integrated Complete-data Likelihood (ICL), which is
#' particularly useful in mixture models.
#'
#' @param object A fitted model of class `gamlssCWM`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A named numeric vector containing AIC, AIC3, AICc, AICu, AWE, BIC, CAIC, and ICL.
#'
#' @export
getIC <- function(object, ...) {
  UseMethod("getIC")
}

#' @rdname getIC
#' @export
getIC.gamlssCWM <- function(object, ...) {
  ll <- logLik(object)
  m <- attr(logLik(object), "df")
  n <- nobs(object)

  # Information criteria
  aic <- -2 * ll + 2 * m
  aic3 <- -2 * ll + 3 * m

  if (n - m - 1 <= 0) {
    aicc <- NA_real_
    aicu <- NA_real_
  } else {
    aicc <- aic + 2 * (m * (m + 1)) / (n - m - 1)
    aicu <- aicc + n * log(n / (n - m - 1))
  }
  caic <- -2 * ll + m * (1 + log(n))
  bic <- -2 * ll + m * log(n)
  awe <- -2 * ll + 2 * m * (1.5 + log(n))

  # Compute ICL (Integrated Complete-data Likelihood) from posterior probabilities.
z <- object$posteriors
map_idx <- apply(z, 1, which.max)
entropy <- sum(log(pmax(z[cbind(seq_len(nrow(z)), map_idx)], .Machine$double.xmin)))
icl <- bic - 2 * entropy

  res <- c(
    AIC = aic, AIC3 = aic3, AICc = aicc, AICu = aicu,
    AWE = awe, BIC = bic, CAIC = caic, ICL = icl
  )
  return(res)
}

#' Extract Model Coefficients from a Cluster-Weighted GAMLSS Model
#'
#' Extracts the estimated coefficients from a fitted \code{gamlssCWM} model.
#' Because the model consists of multiple mixture components, each containing
#' a marginal and a conditional GAMLSS model with potentially multiple 
#' distributional parameters, the output is structured as a hierarchical list.
#'
#' @param object A fitted model of class \code{gamlssCWM}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{pi}}{Numeric vector of the estimated mixing proportions.}
#'     \item{\code{components}}{A list of length \code{k}. Each element corresponds 
#'       to a mixture component and contains sub-lists \code{modelX} and \code{modelY} 
#'       with the extracted coefficients for all active distributional parameters 
#'       (e.g., \code{mu}, \code{sigma}, \code{nu}, \code{tau}).}
#'   }
#'
#' @examples
#' \dontshow{
#' library(gamlss.dist)
#' set.seed(42)
#' df <- data.frame(x = rnorm(100), y = rnorm(100))
#' fit <- gamlssCWM(y ~ x, data = df, famX = NO(), famY = NO(), k = 2, max_iter = 5)
#' }
#' coef(fit)
#'
#' @importFrom stats coef
#' @method coef gamlssCWM
#' @export
coef.gamlssCWM <- function(object, ...) {
  k <- length(object$fits)
  res <- list(pi = object$pi)
  
  comp_list <- vector("list", k)
  names(comp_list) <- paste0("Component_", seq_len(k))
  
  for (j in seq_len(k)) {
    fitX <- object$fits[[j]]$fitX
    fitY <- object$fits[[j]]$fitY
    
    comp_coefs <- list()
    
    # Extract coefficients for the marginal model (X)
    if (!is.null(fitX)) {
      cfX <- list()
      for (p in fitX$parameters) {
        # In gamlss, coef(fit, what="mu") extracts coefficients for parameter 'p'
        cfX[[p]] <- stats::coef(fitX, what = p)
      }
      comp_coefs$modelX <- cfX
    } else {
      comp_coefs$modelX <- NA
    }
    
    # Extract coefficients for the conditional model (Y | X)
    if (!is.null(fitY)) {
      cfY <- list()
      for (p in fitY$parameters) {
        cfY[[p]] <- stats::coef(fitY, what = p)
      }
      comp_coefs$modelY <- cfY
    } else {
      comp_coefs$modelY <- NA
    }
    
    comp_list[[j]] <- comp_coefs
  }
  
  res$components <- comp_list
  return(res)
}