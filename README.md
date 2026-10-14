# gamlssCWM

> Cluster-Weighted GAMLSS Models via the EM Algorithm

## Overview

`gamlssCWM` fits bivariate **Cluster-Weighted Models** (CWM) where each of
the *k* mixture components is described by a pair of
[GAMLSS](https://www.gamlss.com/) distributions:

* a **marginal** GAMLSS model for the covariate **X**, and
* a **conditional** GAMLSS model for the response **Y | X**.

Parameters are estimated via the **Expectation-Maximisation (EM)**
algorithm. Any continuous or discrete count distribution from the
`gamlss.dist` package can be used independently for X and Y|X, and all
four distributional parameters (μ, σ, ν, τ) of the conditional model can
be linked to the covariate through linear or smooth (penalized spline,
polynomial, ...) predictors. Distributions of type `"Mixed"` (e.g.
`BEINF()`) and discrete families with a fixed binomial denominator (e.g.
`BI()`, `BB()`) are not currently supported.

## Installation

```r
# from GitHub (not yet on CRAN)
remotes::install_github("samuxx99/gamlssCWM")

# from a local source directory
devtools::install("/path/to/gamlssCWM")
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

# Point predictions and posterior component probabilities
y_hat <- predict(fit, newdata = df, type = "response")
post  <- predict(fit, newdata = df, type = "posterior")
```

## Main functions and methods

| Function | Description |
|---|---|
| `gamlssCWM()` | Fit the *k*-component CWM via EM |
| `predict()` | Fitted response, posterior probabilities, or joint densities |
| `plot()` | Scatter plot with MAP colouring, regression curves, optional contours |
| `persp()` | 3D perspective plot of the estimated joint density (continuous families) |
| `summary()` / `print()` | Coefficient tables / concise fit overview |
| `coef()` | Estimated coefficients for each component's marginal and conditional model |
| `logLik()`, `AIC()`, `BIC()` | Standard likelihood-based methods, using effective degrees of freedom |
| `getIC()` | Additional information criteria (AIC3, AICc, AICu, AWE, CAIC, ICL) |

## Dependencies

* [`gamlss`](https://CRAN.R-project.org/package=gamlss) (≥ 5.0.0)
* [`gamlss.dist`](https://CRAN.R-project.org/package=gamlss.dist) (≥ 5.0.0)
* `graphics`, `stats`, `utils` (base R)

## Data

The package includes two datasets. `ENRON` is derived from the 'Rsafd' package (MIT licence). `ELECTION` is by Eric Sanders (Radboud University), with TweetGenie age estimates by Dong Nguyen (Utrecht University), redistributed under CC BY 4.0 with their permission. See `inst/COPYRIGHTS`.

## License

GPL (>= 2)
