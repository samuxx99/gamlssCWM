library(gamlssCWM)
library(gamlss.dist)

# ---- Shared fixture ----------------------------------------
# Small but well-separated two-component dataset.
# n=80 keeps the EM fast; set.seed ensures reproducibility.

set.seed(123)
n  <- 80
x  <- c(stats::rnorm(n / 2, mean = -2), stats::rnorm(n / 2, mean = 2))
y  <- c(
  1  + 0.8 * x[seq_len(n / 2)]        + stats::rnorm(n / 2, sd = 0.4),
  -1 - 0.8 * x[seq(n / 2 + 1L, n)]   + stats::rnorm(n / 2, sd = 0.4)
)
df <- data.frame(x = x, y = y)

fit <- gamlssCWM(
  formula  = y ~ x,
  data     = df,
  famX     = NO(),
  famY     = NO(),
  k        = 2L,
  max_iter = 20L
)

# New dataframe for out-of-sample predictions
set.seed(456)
df_new <- data.frame(
  x = stats::rnorm(5, mean = 0),
  y = stats::rnorm(5, mean = 0)
)

# ---- 1. Return class and structure -------------------------
test_that("gamlssCWM returns correct class and structure", {
  expect_s3_class(fit, "gamlssCWM")
  expect_named(fit, c("fits", "pi", "posteriors", "loglik",
                      "ll_history", "converged", "call", "n_retained"),
               ignore.order = TRUE)
})

# ---- 2. Mixing proportions ---------------------------------
test_that("mixing proportions sum to 1", {
  expect_equal(sum(fit$pi), 1, tolerance = 1e-6)
  expect_true(all(fit$pi > 0))
})

# ---- 3. Posterior matrix -----------------------------------
test_that("posterior matrix has correct dimensions and rows sum to 1", {
  expect_equal(dim(fit$posteriors), c(nrow(df), 2L))
  expect_equal(rowSums(fit$posteriors), rep(1, nrow(df)), tolerance = 1e-6)
})

# ---- 4. Log-likelihood is finite and negative --------------
test_that("log-likelihood is finite", {
  expect_true(is.finite(fit$loglik))
})

# ---- 5. ll_history is non-decreasing -----------------------
test_that("ll_history is non-decreasing (EM ascent property)", {
  diffs <- diff(fit$ll_history)
  expect_true(all(diffs >= -1e-4))   # allow tiny numerical noise
})

# ---- 6. predict — response ---------------------------------
test_that("predict type='response' works on new data", {
  yhat <- predict(fit, newdata = df_new, type = "response")
  expect_type(yhat, "double")
  expect_length(yhat, nrow(df_new))
  expect_true(all(is.finite(yhat)))
})

# ---- 7. predict — posterior --------------------------------
test_that("predict type='posterior' works on new data", {
  post <- predict(fit, newdata = df_new, type = "posterior")
  expect_equal(dim(post), c(nrow(df_new), 2L))
  expect_equal(rowSums(post), rep(1, nrow(df_new)), tolerance = 1e-6)
})

# ---- 8. predict — joint ------------------------------------
test_that("predict type='joint' returns non-negative matrix for new data", {
  jnt <- predict(fit, newdata = df_new, type = "joint")
  expect_equal(dim(jnt), c(nrow(df_new), 2L))
  expect_true(all(jnt >= 0))
})

# ---- 9. logLik and nobs ------------------------------------
test_that("logLik and nobs work correctly", {
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_true(attr(ll, "df") > 0)

  expect_equal(nobs(fit), nrow(df))
})

# ---- 10. AIC / BIC -----------------------------------------
test_that("AIC and BIC return finite scalars", {
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
})

# ---- 11. print and summary content -------------------------
test_that("print and summary produce correct output phrases", {
  expect_output(print(fit), "Cluster-Weighted GAMLSS")
  expect_output(print(fit), "Number of components:\\s+2")
  expect_output(summary(fit), "Mixing proportions\\s+:")
  expect_output(summary(fit), "Cluster-Weighted GAMLSS")
})

# ---- 12. plot warning for discrete families ---------------
test_that("plot warns and skips contours with discrete families", {
  set.seed(456)
  n_d <- 60
  x_d <- stats::rnorm(n_d)
  lambda_d <- exp(0.4 + 0.2 * x_d)
  y_d <- stats::rpois(n_d, lambda = lambda_d)
  df_d <- data.frame(x = x_d, y = y_d)

  fit_d <- gamlssCWM(
    formula = y ~ x,
    data = df_d,
    famX = NO(),
    famY = PO(),
    k = 2L,
    max_iter = 10L
  )

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_warning(
    plot(fit_d,
      dataX = df_d$x,
      dataY = df_d$y,
      col = c("steelblue", "tomato"),
      pch = 16,
      contours = TRUE
    ),
    regexp = "Contours are computable only for continuous families"
  )
})

# ---- 12b. plot legend argument ------------------------------
test_that("plot draws a legend at a valid position without error", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_no_error(
    plot(fit,
      dataX = df$x,
      dataY = df$y,
      col = c("steelblue", "tomato"),
      pch = 16,
      legend = "topleft"
    )
  )
})

test_that("plot errors on an invalid legend argument", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_error(
    plot(fit,
      dataX = df$x,
      dataY = df$y,
      col = c("steelblue", "tomato"),
      pch = 16,
      legend = TRUE
    ),
    regexp = "'legend' must be NULL or a legend position string"
  )
})

# ---- 13. Input validation errors ---------------------------
test_that("gamlssCWM validates inputs correctly", {
  expect_error(
    gamlssCWM(y ~ x, data = df, famX = NO(), famY = NO(), k = 0L, max_iter = 5L),
    "'k' must be a positive integer"
  )
  expect_error(
    gamlssCWM(y ~ x, data = df, famX = NO(), famY = NO(), k = 2L, max_iter = 0L),
    "'max_iter' must be a positive integer"
  )
  expect_error(
    gamlssCWM(~ x, data = df, famX = NO(), famY = NO(), k = 2L, max_iter = 5L),
    "'formula' must be a two-sided formula"
  )
  expect_error(
    gamlssCWM(z ~ x, data = df, famX = NO(), famY = NO(), k = 2L, max_iter = 5L),
    "Variables in 'formula' not found"
  )
})

# ---- 14. Additional Formulas -------------------------------
test_that("gamlssCWM handles sigma_formula correctly", {
  fit_sigma <- gamlssCWM(
    formula = y ~ x,
    sigma_formula = ~ x,
    data = df,
    famX = NO(),
    famY = NO(),
    k = 2L,
    max_iter = 10L
  )
  expect_s3_class(fit_sigma, "gamlssCWM")
  
  # Ensure sigma sub-model of the first Y component contains x
  # fits[[1]]$fitY$sigma.coefficients should contain "x"
  expect_false(is.null(fit_sigma$fits[[1]]$fitY))
  expect_true("x" %in% names(fit_sigma$fits[[1]]$fitY$sigma.coefficients))
})

# ---- 15. Alternative Initialization ------------------------
test_that("gamlssCWM works with rand.soft initialization", {
  fit_rand <- gamlssCWM(
    formula = y ~ x,
    data = df,
    famX = NO(),
    famY = NO(),
    k = 2L,
    init = "rand.soft",
    max_iter = 10L
  )
  expect_s3_class(fit_rand, "gamlssCWM")
})

# ---- 16. safe_gamlss failure handling ----------------------
test_that("gamlssCWM stops with an informative error when a fit fails on iteration 1", {
  # LOGNO requires positive response, but df$x has negative values, so
  # safe_gamlss() returns NULL for fitX. With no previous fit to roll
  # back to (first iteration), M_step() stops rather than silently
  # returning a degenerate component.
  suppressWarnings(
    expect_error(
      gamlssCWM(
        formula = y ~ x,
        data = df,
        famX = LOGNO(),
        famY = NO(),
        k = 2L,
        max_iter = 5L
      ),
      "Fit failed critically"
    )
  )
})

# ---- 17. getIC function -----------------------------------
test_that("getIC returns a named vector indicating information criteria", {
  ic_vals <- getIC(fit)
  
  expect_type(ic_vals, "double")
  expect_length(ic_vals, 8)
  expect_named(ic_vals, c("AIC", "AIC3", "AICc", "AICu", "AWE", "BIC", "CAIC", "ICL"))
  expect_true(all(is.finite(ic_vals)))
})

# ---- 18. suppress safe predict warnings -------------------
test_that("predict and plot paths suppress safe refit warning, including splines", {
  skip_if_not(exists("pb", where = asNamespace("gamlss")), "gamlss::pb not available")

  fit_spline <- suppressWarnings(
    gamlssCWM(
      formula = y ~ gamlss::pb(x),
      data = df,
      famX = NO(),
      famY = NO(),
      k = 2L,
      max_iter = 10L
    )
  )

  grid <- data.frame(
    x = seq(min(df$x), max(df$x), length.out = 30),
    y = seq(min(df$y), max(df$y), length.out = 30)
  )

  expect_no_warning(
    predict(fit_spline, newdata = grid, type = "joint")
  )

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_no_warning(
    plot(fit_spline,
      dataX = df$x,
      dataY = df$y,
      col = c("steelblue", "tomato"),
      pch = 16,
      contours = TRUE
    )
  )
})
