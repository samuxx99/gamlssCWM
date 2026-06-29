# gamlssCWM

> Cluster-Weighted GAMLSS Models via the EM Algorithm

## Overview

`gamlssCWM` fits bivariate **Cluster-Weighted Models** (CWM) where each
mixture component is described by a pair of
[GAMLSS](https://www.gamlss.com/) distributions:

* a **marginal** GAMLSS model for the covariate **X**, and
* a **conditional** GAMLSS model for the response **Y | X**.

Parameters are estimated via the **Expectation-Maximisation (EM)**
algorithm. Any distribution from the `gamlss.dist` package can be used
independently for X and Y|X, and all four distributional parameters
(μ, σ, ν, τ) can be modelled as smooth functions of covariates.

## Installation

```r
# from a local source directory
devtools::install("/path/to/gamlssCWM")

# once on CRAN
install.packages("gamlssCWM")
```

## Quick start

```r
library(gamlssCWM)
library(gamlss.dist)

set.seed(42)
n <- 200
x <- c(rnorm(n/2, -2), rnorm(n/2, 2))
y <- c( 1 + 0.8 * x[1:(n/2)]     + rnorm(n/2, sd = 0.5),
       -1 - 0.8 * x[(n/2+1):n]   + rnorm(n/2, sd = 0.5))
df <- data.frame(x = x, y = y)

fit <- gamlssCWM(
  formula  = y ~ x,
  data     = df,
  famX     = NO(),
  famY     = NO(),
  k        = 2,
  max_iter = 50
)

summary(fit)

# Scatter plot coloured by component
plot(fit, dataX = df$x, dataY = df$y,
     col  = c("steelblue", "tomato"),
     pch  = 16,
     main = "gamlssCWM fit")

# Point predictions
y_hat <- predict(fit, newdata = df, type = "response")

# Posterior component probabilities
post <- predict(fit, newdata = df, type = "posterior")
```

## Main functions

| Function | Description |
|---|---|
| `gamlssCWM()` | Fit the CWM via EM |
| `predict.gamlssCWM()` | Predict responses, posteriors or joint densities |
| `plot.gamlssCWM()` | Scatter plot with regression curves (+ optional contours) |
| `summary.gamlssCWM()` | Formatted coefficient tables for each component |

## Dependencies

* [`gamlss`](https://CRAN.R-project.org/package=gamlss) (≥ 5.0.0)
* [`gamlss.dist`](https://CRAN.R-project.org/package=gamlss.dist) (≥ 5.0.0)
* `graphics`, `stats`, `utils` (base R)

## License

GPL-3
