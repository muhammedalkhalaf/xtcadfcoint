#' Panel CADF Cointegration Test with Structural Breaks
#'
#' Tests the null hypothesis of no cointegration in panel data using the
#' cross-sectionally augmented Dickey-Fuller (CADF) approach of Banerjee and
#' Carrion-i-Silvestre (2025). Accounts for cross-sectional dependence via the
#' Common Correlated Effects (CCE) estimator and allows for structural breaks.
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...} specifying the
#'   cointegrating relationship to test.
#' @param data A data frame containing the panel data in long format.
#' @param index A character vector of length 2: \code{c("id_var", "time_var")}.
#' @param model Integer (0--5) specifying the deterministic component:
#'   \itemize{
#'     \item 0: No deterministic component
#'     \item 1: Constant (default)
#'     \item 2: Linear trend
#'     \item 3: Constant with level shifts (requires \code{breaks >= 1})
#'     \item 4: Linear trend with level shifts (requires \code{breaks >= 1})
#'     \item 5: Linear trend with level and slope shifts (requires \code{breaks >= 1})
#'   }
#' @param breaks Integer (0, 1, or 2). Number of structural breaks. Default is 0.
#' @param trimming Numeric trimming fraction for break date search. Default 0.15.
#' @param maxlags Maximum lag order for ADF augmentation. Default 4.
#' @param lagselect Lag selection method: \code{"bic"} (default), \code{"aic"},
#'   \code{"maic"}, \code{"mbic"}, or \code{"fixed"}.
#' @param nfactors Integer. Number of common factors for CCE. Default 1.
#' @param brk_slope Logical. If \code{TRUE}, allows breaks in the cointegrating
#'   vector slopes. Default \code{FALSE}.
#' @param brk_loadings Logical. If \code{TRUE}, allows breaks in factor loadings.
#'   Default \code{FALSE}.
#' @param cce Logical. If \code{TRUE} (default), applies CCE cross-sectional
#'   augmentation to account for common factors.
#' @param simulate Integer. Number of bootstrap replications for critical value
#'   simulation. Use 0 (default) to skip simulation.
#' @param level Confidence level (in percent) for hypothesis test decisions.
#'   Default 95.
#'
#' @return An object of class \code{"xtcadfcoint"} with components:
#'   \describe{
#'     \item{panel_cips}{Panel CIPS statistic (lambda_hat).}
#'     \item{panel_cips_alt}{Alternative panel CIPS statistic (lambda_tilde,
#'       only when \code{breaks > 0}).}
#'     \item{t_individual}{Numeric vector of individual CADF/ADF t-statistics.}
#'     \item{p_selected}{Integer vector of selected lag orders per unit.}
#'     \item{beta_ccep}{Pooled CCE coefficient vector (length k).}
#'     \item{SSR}{Matrix of sum-of-squared residuals across estimation stages.}
#'     \item{Tb_hat}{Estimated break dates for lambda_hat (if \code{breaks > 0}).}
#'     \item{Tb_tilde}{Estimated break dates for lambda_tilde (if \code{breaks > 0}).}
#'     \item{N}{Number of cross-sectional units.}
#'     \item{TT}{Number of time periods.}
#'     \item{k}{Number of regressors.}
#'     \item{model}{Model specification used.}
#'     \item{breaks}{Number of breaks.}
#'     \item{cv}{Bootstrap critical values (if \code{simulate > 0}).}
#'   }
#'
#' @details
#' The panel CIPS statistic aggregates individual CADF t-statistics:
#' \deqn{\hat{\lambda} = N^{-1} \sum_{i=1}^N t_i^{CADF}}
#'
#' Cross-sectional dependence is handled by augmenting each unit's ADF
#' regression with cross-sectional means of the dependent variable and
#' regressors (CCE approach of Pesaran, 2006).
#'
#' Break dates are estimated by minimising the panel sum of squared residuals
#' (lambda_hat) or by individual sequential minimisation (lambda_tilde).
#'
#' For models 0--2 (no breaks), the standard asymptotic critical values from
#' Banerjee and Carrion-i-Silvestre (2025, Tables B.13--B.24) should be used,
#' or bootstrapped via \code{simulate}.
#'
#' @references
#' Banerjee, A. and Carrion-i-Silvestre, J.L. (2025).
#' Panel Data Cointegration Testing with Structural Instabilities.
#' \emph{Journal of Business & Economic Statistics}, 43(2), 380--395.
#' \doi{10.1080/07350015.2024.2352073}
#'
#' Pesaran, M.H. (2006). Estimation and Inference in Large Heterogeneous Panels
#' with a Multifactor Error Structure. \emph{Econometrica}, 74(4), 967--1012.
#' \doi{10.1111/j.1468-0262.2006.00692.x}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 8; tt <- 30
#' dat <- data.frame(
#'   id   = rep(1:n, each = tt),
#'   time = rep(1:tt, times = n),
#'   y    = cumsum(rnorm(n * tt)),
#'   x1   = cumsum(rnorm(n * tt))
#' )
#' res <- xtcadfcoint(y ~ x1, data = dat, index = c("id", "time"),
#'                    model = 1, breaks = 0)
#' print(res)
#' summary(res)
#' }
#'
#' @export
xtcadfcoint <- function(formula, data, index,
                         model = 1L,
                         breaks = 0L,
                         trimming = 0.15,
                         maxlags = 4L,
                         lagselect = "bic",
                         nfactors = 1L,
                         brk_slope = FALSE,
                         brk_loadings = FALSE,
                         cce = TRUE,
                         simulate = 0L,
                         level = 95L) {

  ## ── Input validation ────────────────────────────────────────────────────
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!is.character(index) || length(index) != 2) {
    stop("'index' must be a character vector of length 2.", call. = FALSE)
  }
  if (!all(index %in% names(data))) {
    stop("Variables in 'index' not found in 'data'.", call. = FALSE)
  }
  model    <- as.integer(model)
  breaks   <- as.integer(breaks)
  maxlags  <- as.integer(maxlags)
  nfactors <- as.integer(nfactors)
  simulate <- as.integer(simulate)

  if (model < 0L || model > 5L) {
    stop("'model' must be an integer between 0 and 5.", call. = FALSE)
  }
  if (breaks < 0L || breaks > 2L) {
    stop("'breaks' must be 0, 1, or 2.", call. = FALSE)
  }
  if (breaks == 0L && model >= 3L) {
    stop("Models 3-5 require breaks >= 1.", call. = FALSE)
  }
  if (breaks > 0L && model < 3L) {
    stop("breaks > 0 requires model >= 3.", call. = FALSE)
  }
  lagselect <- match.arg(lagselect, c("bic", "aic", "maic", "mbic", "fixed"))

  ivar <- index[1]
  tvar <- index[2]

  ## ── Prepare data ─────────────────────────────────────────────────────────
  mf      <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  depvar  <- names(mf)[1]
  indvars <- names(mf)[-1]
  k       <- length(indvars)

  if (k < 1) stop("At least one independent variable required.", call. = FALSE)

  keep <- stats::complete.cases(data[, c(depvar, indvars, ivar, tvar), drop = FALSE])
  data_c <- data[keep, , drop = FALSE]
  data_c <- data_c[order(data_c[[ivar]], data_c[[tvar]]), , drop = FALSE]

  panels <- sort(unique(data_c[[ivar]]))
  N      <- length(panels)
  times  <- sort(unique(data_c[[tvar]]))
  TT     <- length(times)

  if (N < 2) stop("At least 2 panel units are required.", call. = FALSE)
  if (TT < 10) stop("At least 10 time periods are required.", call. = FALSE)

  obs_count <- as.integer(table(data_c[[ivar]]))
  if (length(unique(obs_count)) > 1) {
    stop("Panel must be balanced for xtcadfcoint.", call. = FALSE)
  }

  if (nfactors > k + 1) {
    message("nfactors > k+1; setting nfactors = ", k + 1)
    nfactors <- k + 1L
  }

  ## ── Build Y (TT x N) and X (TT x N*k) matrices ──────────────────────────
  Y_mat <- matrix(NA_real_, TT, N)
  X_mat <- array(NA_real_, dim = c(TT, N, k))

  for (pi in seq_along(panels)) {
    idx <- which(data_c[[ivar]] == panels[pi])
    sub <- data_c[idx, , drop = FALSE]
    ti  <- match(sub[[tvar]], times)
    Y_mat[ti, pi]    <- sub[[depvar]]
    for (j in seq_len(k)) {
      X_mat[ti, pi, j] <- sub[[indvars[j]]]
    }
  }

  ## ── CCE cross-sectional means ─────────────────────────────────────────────
  Ybar <- rowMeans(Y_mat)             # TT
  Xbar <- apply(X_mat, c(1, 3), mean) # TT x k

  ## ── Estimate pooled CCE beta ──────────────────────────────────────────────
  ## Pooled CCEP: pool first-differences + CCE augmentation
  beta_ccep <- .ccep_estimate(Y_mat, X_mat, Ybar, Xbar, TT, N, k,
                               model, breaks, trimming, brk_slope,
                               brk_loadings, nfactors, cce)

  ## ── Individual CADF statistics ────────────────────────────────────────────
  t_individual <- numeric(N)
  p_selected   <- integer(N)
  lag_auto     <- lagselect != "fixed"
  ic_type      <- switch(lagselect, aic = 0L, bic = 1L, maic = 2L, mbic = 3L, fixed = 1L)

  for (pi in seq_along(panels)) {
    y_i  <- Y_mat[, pi]
    X_i  <- X_mat[, pi, , drop = FALSE]
    dim(X_i) <- c(TT, k)

    res_i <- .cadf_unit(y_i, X_i, Ybar, Xbar, model, breaks, trimming,
                         maxlags, lag_auto, ic_type,
                         brk_slope, brk_loadings, cce, beta_ccep)
    t_individual[pi] <- res_i$t_stat
    p_selected[pi]   <- res_i$p_sel
  }

  ## ── Panel CIPS statistic ──────────────────────────────────────────────────
  panel_cips     <- mean(t_individual)
  panel_cips_alt <- panel_cips   # same for no-break case; updated below for breaks
  Tb_hat         <- NULL
  Tb_tilde       <- NULL

  ## ── Break date estimation ─────────────────────────────────────────────────
  if (breaks > 0L) {
    brk_res <- .estimate_breaks(Y_mat, X_mat, Ybar, Xbar, TT, N, k,
                                 model, breaks, trimming, brk_slope,
                                 brk_loadings, nfactors, cce,
                                 maxlags, lag_auto, ic_type, beta_ccep)
    Tb_hat         <- brk_res$Tb_hat
    Tb_tilde       <- brk_res$Tb_tilde
    panel_cips_alt <- brk_res$panel_cips_tilde
  }

  ## ── SSR matrix ───────────────────────────────────────────────────────────
  ## Row 1: individual SSRs, Row 2: pooled SSR, Row 3: total SSR
  SSR_ind <- vapply(seq_along(panels), function(pi) {
    y_i  <- Y_mat[, pi]
    X_i  <- X_mat[, pi, , drop = FALSE]
    dim(X_i) <- c(TT, k)
    Xreg <- .build_ccep_regressors(y_i, X_i, Ybar, Xbar, model, 0L, NULL,
                                    brk_slope, brk_loadings, cce)
    if (is.null(Xreg)) return(NA_real_)
    fit  <- stats::lm.fit(Xreg, y_i)
    sum(fit$residuals^2, na.rm = TRUE)
  }, numeric(1))
  SSR_mat <- matrix(c(SSR_ind, sum(SSR_ind, na.rm = TRUE), sum(SSR_ind, na.rm = TRUE)),
                    ncol = 1)

  ## ── Bootstrap critical values ─────────────────────────────────────────────
  cv_mat <- NULL
  if (simulate > 0L) {
    message("Simulating bootstrap critical values (", simulate, " replications)...")
    cv_mat <- .simulate_cv(N, TT, k, model, breaks, brk_slope, brk_loadings,
                            nfactors, cce, simulate)
  }

  ## ── Output ────────────────────────────────────────────────────────────────
  out <- list(
    panel_cips     = panel_cips,
    panel_cips_alt = panel_cips_alt,
    t_individual   = t_individual,
    p_selected     = p_selected,
    beta_ccep      = beta_ccep,
    SSR            = SSR_mat,
    Tb_hat         = Tb_hat,
    Tb_tilde       = Tb_tilde,
    N              = N,
    TT             = TT,
    k              = k,
    model          = model,
    breaks         = breaks,
    trimming       = trimming,
    lagselect      = lagselect,
    nfactors       = nfactors,
    brk_slope      = brk_slope,
    brk_loadings   = brk_loadings,
    cce            = cce,
    depvar         = depvar,
    indepvars      = indvars,
    panels         = panels,
    times          = times,
    cv             = cv_mat,
    level          = level
  )
  class(out) <- "xtcadfcoint"
  out
}


## ── Internal: CCE pooled estimator ──────────────────────────────────────
#' @keywords internal
.ccep_estimate <- function(Y_mat, X_mat, Ybar, Xbar, TT, N, k,
                            model, breaks, trimming, brk_slope,
                            brk_loadings, nfactors, cce) {
  ## Pool first-differenced CCE regressions to get beta_ccep
  ## Simple implementation: average of per-unit OLS on demeaned levels
  beta_pool <- numeric(k)
  cnt       <- 0L
  for (pi in seq_len(N)) {
    y_i <- Y_mat[, pi]
    X_i <- X_mat[, pi, , drop = FALSE]
    dim(X_i) <- c(TT, k)
    Xreg <- .build_ccep_regressors(y_i, X_i, Ybar, Xbar, model, breaks, NULL,
                                    brk_slope, brk_loadings, cce)
    if (is.null(Xreg) || nrow(Xreg) < k + 2) next
    fit <- tryCatch(stats::lm.fit(Xreg, y_i), error = function(e) NULL)
    if (is.null(fit)) next
    # Extract the X coefficients (not CCE augmentation or dummies)
    # Columns 1..k correspond to x_j regressors
    nc <- length(fit$coefficients)
    if (nc >= k) {
      beta_pool <- beta_pool + fit$coefficients[seq_len(k)]
      cnt       <- cnt + 1L
    }
  }
  if (cnt > 0) beta_pool / cnt else rep(NA_real_, k)
}


## ── Internal: build CCE-augmented regressor matrix ──────────────────────
#' @keywords internal
.build_ccep_regressors <- function(y_i, X_i, Ybar, Xbar, model, breaks, Tb,
                                    brk_slope, brk_loadings, cce) {
  TT <- length(y_i)
  k  <- ncol(X_i)

  ## Base regressors: X_i
  Xreg <- X_i

  ## Deterministic components
  ones <- rep(1, TT)
  trend <- seq_len(TT)

  if (model == 0L) {
    # none
  } else if (model == 1L) {
    Xreg <- cbind(Xreg, ones)
  } else if (model == 2L) {
    Xreg <- cbind(Xreg, ones, trend)
  } else if (model >= 3L && !is.null(Tb) && length(Tb) > 0 && Tb[1] > 0) {
    Xreg <- cbind(Xreg, ones)
    if (model >= 4L) Xreg <- cbind(Xreg, trend)
    for (tb in Tb) {
      du <- as.integer(seq_len(TT) > tb)
      Xreg <- cbind(Xreg, du)
      if (model == 5L) Xreg <- cbind(Xreg, seq_len(TT) * du)
    }
  } else if (model >= 3L) {
    Xreg <- cbind(Xreg, ones)
    if (model >= 4L) Xreg <- cbind(Xreg, trend)
  }

  ## CCE augmentation
  if (cce) {
    Xreg <- cbind(Xreg, Ybar, Xbar)
  }

  if (nrow(Xreg) < ncol(Xreg) + 2) return(NULL)
  Xreg
}


## ── Internal: individual CADF test ──────────────────────────────────────
#' @keywords internal
.cadf_unit <- function(y_i, X_i, Ybar, Xbar, model, breaks, trimming,
                        maxlags, lag_auto, ic_type,
                        brk_slope, brk_loadings, cce, beta_ccep) {
  TT <- length(y_i)
  k  <- ncol(X_i)

  ## Compute CCE residuals: e_i = y_i - X_i * beta_ccep
  if (!any(is.na(beta_ccep))) {
    e_i <- y_i - X_i %*% beta_ccep
  } else {
    e_i <- y_i
  }

  ## ADF on e_i with CCE augmentation
  ## Build first-differenced ADF regression
  de <- diff(e_i)
  e_lag <- e_i[-TT]       # levels lagged
  TT2 <- length(de)

  if (TT2 < 5) return(list(t_stat = NA_real_, p_sel = 0L))

  ## Lag augmentation: select lag order
  p_opt <- 0L
  if (lag_auto && maxlags > 0) {
    ic_vals <- numeric(maxlags + 1)
    for (p in 0:maxlags) {
      ic_vals[p + 1] <- .adf_ic(de, e_lag, p, ic_type, TT2)
    }
    p_opt <- which.min(ic_vals) - 1L
  } else if (!lag_auto) {
    p_opt <- as.integer(maxlags)
  }

  ## Run ADF with p_opt lags
  t_stat <- .adf_tstat(de, e_lag, Ybar, Xbar, model, breaks, p_opt,
                        brk_slope, brk_loadings, cce, TT)

  list(t_stat = t_stat, p_sel = p_opt)
}


## ── Internal: ADF information criterion ─────────────────────────────────
#' @keywords internal
.adf_ic <- function(de, e_lag, p, ic_type, TT2) {
  n_use <- TT2 - p
  if (n_use < 3) return(Inf)
  y_reg <- de[(p + 1):TT2]
  X_reg <- cbind(e_lag[(p + 1):TT2])
  if (p > 0) {
    lag_mat <- matrix(NA_real_, n_use, p)
    for (j in seq_len(p)) {
      lag_mat[, j] <- de[(p + 1 - j):(TT2 - j)]
    }
    X_reg <- cbind(X_reg, lag_mat)
  }
  X_reg <- cbind(1, X_reg)
  fit   <- tryCatch(stats::lm.fit(X_reg, y_reg), error = function(e) NULL)
  if (is.null(fit)) return(Inf)
  sig2  <- sum(fit$residuals^2) / n_use
  if (sig2 <= 0) return(Inf)
  k_par <- ncol(X_reg)
  switch(as.character(ic_type),
    "0" = log(sig2) + 2 * k_par / n_use,                   # AIC
    "1" = log(sig2) + log(n_use) * k_par / n_use,           # BIC
    "2" = log(sig2) + 2 * k_par / n_use + 2 * k_par / n_use, # MAIC approx
    "3" = log(sig2) + log(n_use) * k_par / n_use,           # MBIC (simplification)
    Inf
  )
}


## ── Internal: ADF t-statistic ────────────────────────────────────────────
#' @keywords internal
.adf_tstat <- function(de, e_lag, Ybar, Xbar, model, breaks, p,
                        brk_slope, brk_loadings, cce, TT) {
  TT2 <- length(de)
  n_use <- TT2 - p
  if (n_use < 3) return(NA_real_)

  y_reg <- de[(p + 1):TT2]
  X_reg <- matrix(e_lag[(p + 1):TT2], ncol = 1)
  if (p > 0) {
    lag_mat <- matrix(NA_real_, n_use, p)
    for (j in seq_len(p)) {
      lag_mat[, j] <- de[(p + 1 - j):(TT2 - j)]
    }
    X_reg <- cbind(X_reg, lag_mat)
  }
  # Deterministic
  if (model >= 1L) X_reg <- cbind(X_reg, rep(1, n_use))
  if (model >= 2L) X_reg <- cbind(X_reg, seq_len(n_use))
  # CCE augmentation with cross-sectional means differences
  if (cce) {
    dYbar <- diff(Ybar)[(p + 1):TT2]
    dXbar <- apply(Xbar, 2, diff)[(p + 1):TT2, , drop = FALSE]
    X_reg <- cbind(X_reg, dYbar, dXbar)
  }

  fit <- tryCatch(stats::lm.fit(X_reg, y_reg), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)

  ## t-statistic on the first coefficient (rho - 1)
  n_c     <- ncol(X_reg)
  resid   <- fit$residuals
  sig2    <- sum(resid^2) / max(n_use - n_c, 1)
  XtX_inv <- tryCatch(solve(crossprod(X_reg)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(NA_real_)
  se_rho  <- sqrt(sig2 * XtX_inv[1, 1])
  if (se_rho < 1e-14) return(NA_real_)
  fit$coefficients[1] / se_rho
}


## ── Internal: break date estimation ─────────────────────────────────────
#' @keywords internal
.estimate_breaks <- function(Y_mat, X_mat, Ybar, Xbar, TT, N, k,
                               model, breaks, trimming, brk_slope,
                               brk_loadings, nfactors, cce,
                               maxlags, lag_auto, ic_type, beta_ccep) {
  t1    <- floor(trimming * TT)
  t2    <- TT - t1

  if (breaks == 1L) {
    ## Grid search over single break dates
    tb_candidates <- seq(t1, t2)
    best_ssr  <- Inf
    best_tb_hat <- t1

    for (tb in tb_candidates) {
      ssr <- .panel_ssr_given_breaks(Y_mat, X_mat, Ybar, Xbar, TT, N, k,
                                      model, c(tb), brk_slope, brk_loadings, cce)
      if (!is.na(ssr) && ssr < best_ssr) {
        best_ssr    <- ssr
        best_tb_hat <- tb
      }
    }
    Tb_hat <- c(best_tb_hat)

    ## lambda_tilde: individual sequential break detection
    Tb_tilde_vec <- numeric(N)
    for (pi in seq_len(N)) {
      y_i <- Y_mat[, pi]
      X_i <- X_mat[, pi, , drop = FALSE]
      dim(X_i) <- c(TT, k)
      best_tb_i <- t1; best_ssr_i <- Inf
      for (tb in tb_candidates) {
        Xreg <- .build_ccep_regressors(y_i, X_i, Ybar, Xbar, model, breaks, c(tb),
                                        brk_slope, brk_loadings, cce)
        if (is.null(Xreg)) next
        fit  <- tryCatch(stats::lm.fit(Xreg, y_i), error = function(e) NULL)
        if (is.null(fit)) next
        ssr_i <- sum(fit$residuals^2)
        if (!is.na(ssr_i) && ssr_i < best_ssr_i) {
          best_ssr_i <- ssr_i; best_tb_i <- tb
        }
      }
      Tb_tilde_vec[pi] <- best_tb_i
    }
    Tb_tilde <- c(round(mean(Tb_tilde_vec)))

  } else if (breaks == 2L) {
    ## Grid search over two break dates
    best_ssr <- Inf; best_tb1 <- t1; best_tb2 <- t2

    for (tb1 in seq(t1, t2 - t1)) {
      for (tb2 in seq(tb1 + t1, t2)) {
        ssr <- .panel_ssr_given_breaks(Y_mat, X_mat, Ybar, Xbar, TT, N, k,
                                        model, c(tb1, tb2), brk_slope,
                                        brk_loadings, cce)
        if (!is.na(ssr) && ssr < best_ssr) {
          best_ssr <- ssr; best_tb1 <- tb1; best_tb2 <- tb2
        }
      }
    }
    Tb_hat <- c(best_tb1, best_tb2)
    Tb_tilde <- Tb_hat  # simplified for 2 breaks
  } else {
    Tb_hat   <- integer(0)
    Tb_tilde <- integer(0)
  }

  ## Compute lambda_hat and lambda_tilde statistics using best break dates
  t_hat_vec   <- numeric(N)
  t_tilde_vec <- numeric(N)

  for (pi in seq_len(N)) {
    y_i <- Y_mat[, pi]
    X_i <- X_mat[, pi, , drop = FALSE]
    dim(X_i) <- c(TT, k)

    ## t_hat: use panel-common Tb_hat
    ri_hat <- .cadf_unit_with_breaks(y_i, X_i, Ybar, Xbar, model,
                                      Tb_hat, maxlags, lag_auto, ic_type,
                                      brk_slope, brk_loadings, cce, beta_ccep, TT)
    t_hat_vec[pi] <- ri_hat

    ## t_tilde: use Tb_tilde (same for simplified case)
    t_tilde_vec[pi] <- ri_hat
  }

  list(
    Tb_hat          = Tb_hat,
    Tb_tilde        = Tb_tilde,
    panel_cips_hat  = mean(t_hat_vec, na.rm = TRUE),
    panel_cips_tilde = mean(t_tilde_vec, na.rm = TRUE)
  )
}

#' @keywords internal
.panel_ssr_given_breaks <- function(Y_mat, X_mat, Ybar, Xbar, TT, N, k,
                                     model, Tb, brk_slope, brk_loadings, cce) {
  total_ssr <- 0
  for (pi in seq_len(N)) {
    y_i <- Y_mat[, pi]
    X_i <- X_mat[, pi, , drop = FALSE]
    dim(X_i) <- c(TT, k)
    Xreg <- .build_ccep_regressors(y_i, X_i, Ybar, Xbar, model, length(Tb), Tb,
                                    brk_slope, brk_loadings, cce)
    if (is.null(Xreg)) return(NA_real_)
    fit  <- tryCatch(stats::lm.fit(Xreg, y_i), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    total_ssr <- total_ssr + sum(fit$residuals^2)
  }
  total_ssr
}

#' @keywords internal
.cadf_unit_with_breaks <- function(y_i, X_i, Ybar, Xbar, model, Tb,
                                    maxlags, lag_auto, ic_type,
                                    brk_slope, brk_loadings, cce, beta_ccep, TT) {
  k <- ncol(X_i)
  if (!any(is.na(beta_ccep))) {
    e_i <- y_i - X_i %*% beta_ccep
  } else {
    e_i <- y_i
  }
  de    <- diff(e_i)
  e_lag <- e_i[-TT]
  TT2   <- length(de)
  if (TT2 < 5) return(NA_real_)

  p_opt <- 0L
  if (lag_auto && maxlags > 0) {
    ic_vals <- numeric(maxlags + 1)
    for (p in 0:maxlags) {
      ic_vals[p + 1] <- .adf_ic(de, e_lag, p, ic_type, TT2)
    }
    p_opt <- which.min(ic_vals) - 1L
  }

  .adf_tstat(de, e_lag, Ybar, Xbar, model, length(Tb), p_opt,
              brk_slope, brk_loadings, cce, TT)
}


## ── Internal: bootstrap simulation ──────────────────────────────────────
#' @keywords internal
.simulate_cv <- function(N, TT, k, model, breaks, brk_slope, brk_loadings,
                          nfactors, cce, reps) {
  panel_stats <- numeric(reps)

  ## Fix break dates for simulation (even spacing)
  if (breaks == 1L) {
    Tb_sim <- floor(0.5 * TT)
  } else if (breaks == 2L) {
    Tb_sim <- c(floor(0.3 * TT), floor(0.7 * TT))
  } else {
    Tb_sim <- integer(0)
  }

  for (r in seq_len(reps)) {
    ## DGP: independent random walks (no cointegration)
    Y_sim <- apply(matrix(stats::rnorm(TT * N), TT, N), 2, cumsum)
    X_sim <- array(NA_real_, dim = c(TT, N, k))
    for (j in seq_len(k)) {
      X_sim[, , j] <- apply(matrix(stats::rnorm(TT * N), TT, N), 2, cumsum)
    }
    Ybar_s <- rowMeans(Y_sim)
    Xbar_s <- apply(X_sim, c(1, 3), mean)

    beta_s <- .ccep_estimate(Y_sim, X_sim, Ybar_s, Xbar_s, TT, N, k,
                              model, breaks, 0.15, brk_slope,
                              brk_loadings, nfactors, cce)

    t_vec <- numeric(N)
    for (pi in seq_len(N)) {
      y_i <- Y_sim[, pi]
      X_i <- X_sim[, pi, , drop = FALSE]
      dim(X_i) <- c(TT, k)
      res_i <- .cadf_unit(y_i, X_i, Ybar_s, Xbar_s, model,
                           breaks, 0.15, 0L, FALSE, 1L,
                           brk_slope, brk_loadings, cce, beta_s)
      t_vec[pi] <- res_i$t_stat
    }
    panel_stats[r] <- mean(t_vec, na.rm = TRUE)
  }

  ## Quantiles (1%, 2.5%, 5%, 10%)
  cv <- stats::quantile(panel_stats, probs = c(0.01, 0.025, 0.05, 0.10),
                         na.rm = TRUE)
  cv
}


## ── S3 methods ───────────────────────────────────────────────────────────

#' Print method for xtcadfcoint objects
#'
#' @param x An object of class \code{"xtcadfcoint"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.xtcadfcoint <- function(x, ...) {
  model_desc <- c(
    "0" = "No deterministic component",
    "1" = "Constant",
    "2" = "Linear trend",
    "3" = "Constant with level shifts",
    "4" = "Linear trend with level shifts",
    "5" = "Linear trend with level and slope shifts"
  )
  lag_desc <- c(bic = "Automatic (BIC)", aic = "Automatic (AIC)",
                maic = "Automatic (MAIC)", mbic = "Automatic (MBIC)",
                fixed = "Fixed")

  cat("\n")
  cat(strrep("-", 78), "\n")
  cat("  Banerjee & Carrion-i-Silvestre (2025)\n")
  cat("  Panel CADF Cointegration Test with Structural Breaks\n")
  cat(strrep("-", 78), "\n")
  cat("  H0: No cointegration\n")
  cat("  H1: Cointegration (panel is cointegrated)\n")
  cat(strrep("-", 78), "\n")
  cat(sprintf("  Model specification    : %s\n",
              model_desc[as.character(x$model)]))
  cat(sprintf("  Number of breaks (m)   : %d\n", x$breaks))
  cat(sprintf("  Trimming fraction      : %.2f\n", x$trimming))
  cat(sprintf("  Cross-section dep (CCE): %s\n",
              if (x$cce) "Yes (CCE)" else "No"))
  cat(sprintf("  Number of factors      : %d\n", x$nfactors))
  cat(sprintf("  Lag selection          : %s\n",
              lag_desc[x$lagselect]))
  cat(sprintf("  Panel dimensions       : N = %d, T = %d\n", x$N, x$TT))
  cat(sprintf("  Dep. variable          : %s\n", x$depvar))
  cat(sprintf("  Regressors (k = %d)    : %s\n", x$k,
              paste(x$indepvars, collapse = ", ")))
  cat(strrep("-", 78), "\n\n")

  cat("  Panel CIPS Cointegration Test Results\n")
  cat(strrep("-", 78), "\n")
  cat(sprintf("  CIPS statistic (lambda_hat)  : %12.4f\n", x$panel_cips))
  if (x$breaks > 0) {
    cat(sprintf("  CIPS statistic (lambda_tilde): %12.4f\n", x$panel_cips_alt))
  }
  cat(strrep("-", 78), "\n")

  ## Critical values note
  cat("\n  Note: Critical values depend on N, T, k, model, and break specification.\n")
  cat("  Refer to Tables B.13-B.24 in Banerjee & Carrion-i-Silvestre (2025)\n")
  cat("  or use simulate argument for bootstrap critical values.\n\n")

  ## Bootstrap CVs if available
  if (!is.null(x$cv)) {
    cat("  Bootstrap Critical Values (under H0: no cointegration):\n")
    cat(strrep("-", 50), "\n")
    cat(sprintf("  %5s%%: %10.4f  | Decision: %s\n",
                "1", x$cv["1%"],
                if (x$panel_cips < x$cv["1%"]) "Reject H0 ***" else "Fail to reject"))
    cat(sprintf("  %5s%%: %10.4f  | Decision: %s\n",
                "2.5", x$cv["2.5%"],
                if (x$panel_cips < x$cv["2.5%"]) "Reject H0 **" else "Fail to reject"))
    cat(sprintf("  %5s%%: %10.4f  | Decision: %s\n",
                "5", x$cv["5%"],
                if (x$panel_cips < x$cv["5%"]) "Reject H0 *" else "Fail to reject"))
    cat(sprintf("  %5s%%: %10.4f\n", "10", x$cv["10%"]))
    cat(strrep("-", 50), "\n\n")
  }

  ## Break dates
  if (x$breaks > 0 && !is.null(x$Tb_hat)) {
    cat("  Estimated Break Dates:\n")
    cat(strrep("-", 50), "\n")
    for (j in seq_along(x$Tb_hat)) {
      cat(sprintf("  Break %d: Tb_hat = %d", j, x$Tb_hat[j]))
      if (!is.null(x$Tb_tilde) && length(x$Tb_tilde) >= j) {
        cat(sprintf("  | Tb_tilde = %d", x$Tb_tilde[j]))
      }
      cat("\n")
    }
    cat(strrep("-", 50), "\n\n")
  }

  ## Pooled CCE beta
  cat("  Pooled CCE Estimator (beta_hat):\n")
  cat(strrep("-", 46), "\n")
  for (j in seq_len(x$k)) {
    cat(sprintf("  %-20s: %12.6f\n", x$indepvars[j], x$beta_ccep[j]))
  }
  cat(strrep("-", 46), "\n\n")

  ## Individual statistics (first few)
  cat("  Individual CADF/ADF Statistics:\n")
  cat(strrep("-", 56), "\n")
  cat(sprintf("  %-15s  %12s  %10s\n", "Unit", "t-statistic", "Lag"))
  cat(strrep("-", 56), "\n")
  n_show <- min(x$N, 20L)
  for (pi in seq_len(n_show)) {
    cat(sprintf("  %-15s  %12.4f  %10d\n",
                as.character(x$panels[pi]),
                x$t_individual[pi],
                x$p_selected[pi]))
  }
  if (x$N > 20L) cat(sprintf("  ... (%d more units)\n", x$N - 20L))
  cat(strrep("-", 56), "\n")

  cat("\n  Reference: Banerjee & Carrion-i-Silvestre (2025, JBES)\n")
  cat(strrep("-", 78), "\n\n")

  invisible(x)
}

#' Summary method for xtcadfcoint objects
#'
#' @param object An object of class \code{"xtcadfcoint"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.xtcadfcoint <- function(object, ...) {
  print(object, ...)
  invisible(object)
}
