# xtcadfcoint 1.0.0

* Initial CRAN release.
* Implements the panel CADF cointegration test of
  Banerjee & Carrion-i-Silvestre (2025, Journal of Business & Economic Statistics).
* Supports models 0-5: no deterministic, constant, trend, and break variants.
* Allows 0, 1, or 2 structural breaks in cointegrating vectors, factor loadings.
* CCE augmentation for cross-sectional dependence (Pesaran 2006).
* Automatic lag selection: AIC, BIC, MAIC, MBIC.
* Bootstrap critical value simulation via `simulate` argument.
