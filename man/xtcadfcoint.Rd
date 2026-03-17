\name{xtcadfcoint}
\alias{xtcadfcoint}
\title{Panel CADF Cointegration Test with Structural Breaks}
\description{
  Tests the null hypothesis of no cointegration in panel data using the
  cross-sectionally augmented Dickey-Fuller (CADF) approach of Banerjee and
  Carrion-i-Silvestre (2025). Accounts for cross-sectional dependence via the
  Common Correlated Effects (CCE) estimator and allows for structural breaks.
}
\usage{
xtcadfcoint(formula, data, index, model = 1, breaks = 0,
            trimming = 0.15, maxlags = 4, lagselect = "bic",
            nfactors = 1, brk_slope = FALSE, brk_loadings = FALSE,
            cce = TRUE, simulate = 0, level = 95)
}
\arguments{
  \item{formula}{A formula of the form \code{y ~ x1 + x2 + ...}.}
  \item{data}{A data frame in long format.}
  \item{index}{Character vector of length 2: \code{c("id_var", "time_var")}.}
  \item{model}{Integer (0-5) specifying the deterministic component.}
  \item{breaks}{Integer (0, 1, or 2). Number of structural breaks.}
  \item{trimming}{Numeric. Trimming fraction for break date search. Default 0.15.}
  \item{maxlags}{Maximum lag order for ADF augmentation. Default 4.}
  \item{lagselect}{Lag selection: \code{"bic"}, \code{"aic"}, \code{"maic"},
    \code{"mbic"}, or \code{"fixed"}.}
  \item{nfactors}{Integer. Number of common factors for CCE. Default 1.}
  \item{brk_slope}{Logical. Allow breaks in cointegrating vector slopes.}
  \item{brk_loadings}{Logical. Allow breaks in factor loadings.}
  \item{cce}{Logical. Apply CCE augmentation. Default \code{TRUE}.}
  \item{simulate}{Integer. Bootstrap replications for critical values (0 = skip).}
  \item{level}{Confidence level in percent. Default 95.}
}
\value{
  An object of class \code{"xtcadfcoint"} containing panel CIPS statistics,
  individual CADF statistics, estimated break dates, and coefficient estimates.
}
\references{
Banerjee, A. and Carrion-i-Silvestre, J.L. (2025). Testing for panel cointegration
with structural breaks. \emph{Journal of Econometrics}.
}
\examples{
\donttest{
set.seed(42)
n <- 5; tt <- 30
uid  <- rep(1:n, each = tt)
tval <- rep(1:tt, times = n)
x_it <- cumsum(rnorm(n * tt))
y    <- x_it + rnorm(n * tt, sd = 0.5)
dat  <- data.frame(id = uid, time = tval, y = y, x = x_it)
res  <- xtcadfcoint(y ~ x, data = dat, index = c("id", "time"), model = 1)
print(res)
}
}
\keyword{htest}
