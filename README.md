# xtcadfcoint

**Panel CADF Cointegration Test with Structural Breaks**

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/xtcadfcoint)](https://CRAN.R-project.org/package=xtcadfcoint)
<!-- badges: end -->

## Overview

`xtcadfcoint` tests the null hypothesis of no cointegration in balanced panel
data using cross-sectionally augmented Dickey-Fuller (CADF) regressions
following Banerjee and Carrion-i-Silvestre (2025, *Journal of Business &
Economic Statistics*).

Key features:

- **6 model specifications**: no constant, constant, trend, constant + breaks,
  trend + breaks, trend + level and slope breaks
- **0, 1, or 2 structural breaks** with endogenous break date estimation
- **CCE augmentation** for cross-sectional dependence (Pesaran, 2006)
- **Automatic lag selection**: AIC, BIC, MAIC, MBIC
- **Bootstrap critical value simulation**

## Installation

```r
install.packages("xtcadfcoint")
```

## Usage

```r
library(xtcadfcoint)

set.seed(42)
n <- 8; tt <- 30
dat <- data.frame(
  id   = rep(1:n, each = tt),
  time = rep(1:tt, times = n),
  y    = c(apply(matrix(rnorm(n * tt), tt, n), 2, cumsum)),
  x1   = c(apply(matrix(rnorm(n * tt), tt, n), 2, cumsum))
)

# No breaks, constant model
res <- xtcadfcoint(y ~ x1, data = dat, index = c("id", "time"),
                   model = 1, breaks = 0)
print(res)

# With one structural break
res_brk <- xtcadfcoint(y ~ x1, data = dat, index = c("id", "time"),
                        model = 3, breaks = 1)
print(res_brk)
```

## Reference

Banerjee, A. and Carrion-i-Silvestre, J.L. (2025).
Panel Data Cointegration Testing with Structural Instabilities.
*Journal of Business & Economic Statistics*, 43(2), 380–395.
<https://doi.org/10.1080/07350015.2024.2352073>
