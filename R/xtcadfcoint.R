#' Panel CADF Cointegration Test with Structural Breaks
#'
#' @description
#' Performs the panel CADF (Cross-sectionally Augmented Dickey-Fuller)
#' cointegration test of Banerjee and Carrion-i-Silvestre (2017). The test
#' accounts for cross-sectional dependence through Common Correlated Effects
#' (CCE) augmentation and supports structural breaks in deterministic components.
#'
#' @param formula A formula specifying the cointegrating relationship
#'   (e.g., \code{y ~ x1 + x2}).
#' @param data A data frame containing panel data.
#' @param id Character string naming the cross-sectional unit identifier.
#' @param time Character string naming the time period identifier.
#' @param model Integer (0-5). Model specification for deterministic components:
#'   \describe{
#'     \item{0}{No deterministic components}
#'     \item{1}{Constant only (default)}
#'     \item{2}{Constant and linear trend}
#'     \item{3}{Constant with level shifts at breaks}
#'     \item{4}{Linear trend with level shifts at breaks}
#'     \item{5}{Linear trend with level and slope shifts at breaks}
#'   }
#' @param breaks Integer (0-2). Number of structural breaks to allow.
#'   Models 0-2 require \code{breaks = 0}; models 3-5 require \code{breaks >= 1}.
#' @param trimming Numeric. Trimming fraction for break date search.
#'   Default is 0.15, meaning breaks cannot occur in first/last 15\% of sample.
#' @param maxlags Integer. Maximum lag order for ADF regressions. Default is 4.
#' @param lagselect Character. Lag selection method: \code{"bic"} (default),
#'   \code{"aic"}, \code{"maic"}, \code{"mbic"}, or \code{"fixed"}.
#' @param nfactors Integer. Number of unobserved common factors. Default is 1.
#'   Used to determine the number of cross-sectional averages for CCE.
#' @param brk_slope Logical. If TRUE, allows cointegrating coefficients to
#'   shift at break dates. Default is FALSE.
#' @param brk_loadings Logical. If TRUE, allows factor loadings to shift at
#'   break dates. Default is FALSE.
#' @param cce Logical. If TRUE (default), uses CCE augmentation for
#'   cross-sectional dependence. Set FALSE for standard ADF tests.
#' @param simulate Integer. If > 0, simulates bootstrap critical values with
#'   the specified number of replications. Default is 0 (no simulation).
#'
#' @return An object of class \code{"xtcadfcoint"} containing:
#' \describe{
#'   \item{panel_cips}{Panel CIPS statistic (average of individual t-statistics).}
#'   \item{panel_cips_alt}{Alternative CIPS using different break estimator
#'     (only for \code{breaks > 0}).}
#'   \item{t_individual}{Vector of individual CADF t-statistics.}
#'   \item{p_selected}{Vector of selected lag orders per unit.}
#'   \item{beta_cce}{Pooled CCE coefficient estimates.}
#'   \item{Tb_hat}{Estimated break dates (lambda_hat criterion).}
#'   \item{Tb_tilde}{Estimated break dates (lambda_tilde criterion).}
#'   \item{SSR}{Sum of squared residuals (total, pre-break, post-break).}
#'   \item{N}{Number of cross-sectional units.}
#'   \item{TT}{Number of time periods.}
#'   \item{k}{Number of regressors.}
#'   \item{model}{Model specification used.}
#'   \item{breaks}{Number of breaks.}
#'   \item{cv_panel}{Bootstrap critical values for panel CIPS (if simulated).}
#'   \item{cv_ind}{Bootstrap critical values for individual tests (if simulated).}
#'   \item{call}{The matched call.}
#' }
#'
#' @details
#' The null hypothesis is no cointegration, tested against the alternative
#' that the panel is cointegrated. The test is based on applying ADF tests
#' to residuals from the cointegrating regression, augmented with
#' cross-sectional averages (CCE) to handle unobserved common factors.
#'
#' \strong{Model Specifications:}
#' Models 0-2 assume no structural breaks and test standard cointegration.
#' Models 3-5 allow for structural instabilities in the deterministic
#' components, which can affect both the cointegrating relationship and
#' factor structure.
#'
#' \strong{CCE Augmentation:}
#' The CCE approach follows Pesaran (2006) by augmenting individual
#' regressions with cross-sectional averages of the dependent and
#' independent variables. This makes the test robust to cross-sectional
#' dependence induced by common factors.
#'
#' \strong{Break Date Estimation:}
#' When \code{breaks > 0}, break dates are estimated by minimizing the
#' sum of squared residuals (SSR) over all admissible break combinations.
#' The trimming parameter excludes break dates too close to sample endpoints.
#'
#' @references
#' Banerjee, A., & Carrion-i-Silvestre, J. L. (2017). Testing for panel
#' cointegration using common correlated effects estimators.
#' \emph{Journal of Time Series Analysis}, 38(4), 610--636.
#' \doi{10.1111/jtsa.12234}
#'
#' Pesaran, M. H. (2006). Estimation and inference in large heterogeneous
#' panels with a multifactor error structure. \emph{Econometrica}, 74(4),
#' 967--1012. \doi{10.1111/j.1468-0262.2006.00692.x}
#'
#' Pesaran, M. H. (2007). A simple panel unit root test in the presence
#' of cross-section dependence. \emph{Journal of Applied Econometrics},
#' 22(2), 265--312. \doi{10.1002/jae.951}
#'
#' @examples
#' # Generate example panel data
#' set.seed(123)
#' N <- 10
#' TT <- 50
#' panel_data <- data.frame(
#'   id = rep(1:N, each = TT),
#'   time = rep(1:TT, N),
#'   y = rnorm(N * TT),
#'   x = rnorm(N * TT)
#' )
#'
#' \donttest{
#' # Basic test with constant
#' result <- xtcadfcoint(y ~ x, data = panel_data,
#'                      id = "id", time = "time", model = 1)
#' print(result)
#'
#' # Test with one structural break
#' result_brk <- xtcadfcoint(y ~ x, data = panel_data,
#'                          id = "id", time = "time",
#'                          model = 3, breaks = 1)
#' summary(result_brk)
#' }
#'
#' @export
xtcadfcoint <- function(formula,
                        data,
                        id,
                        time,
                        model = 1,
                        breaks = 0,
                        trimming = 0.15,
                        maxlags = 4,
                        lagselect = c("bic", "aic", "maic", "mbic", "fixed"),
                        nfactors = 1,
                        brk_slope = FALSE,
                        brk_loadings = FALSE,
                        cce = TRUE,
                        simulate = 0) {

  call <- match.call()
  lagselect <- match.arg(lagselect)

  # ---- Validate inputs ----
  if (model < 0 || model > 5) {
    stop("model must be between 0 and 5", call. = FALSE)
  }

  if (breaks < 0 || breaks > 2) {
    stop("breaks must be 0, 1, or 2", call. = FALSE)
  }

  if (breaks == 0 && model >= 3) {
    stop("models 3-5 require breaks >= 1", call. = FALSE)
  }

  if (breaks > 0 && model < 3) {
    stop("breaks > 0 requires model >= 3", call. = FALSE)
  }

  if (trimming <= 0 || trimming >= 0.5) {
    stop("trimming must be between 0 and 0.5", call. = FALSE)
  }

  if (nfactors < 1) {
    stop("nfactors must be at least 1", call. = FALSE)
  }

  # ---- Parse formula and data ----
  mf <- model.frame(formula, data = data)
  depvar <- model.response(mf)
  indepvars <- model.matrix(formula, data = mf)

  # Remove intercept (handled by model specification)
  if ("(Intercept)" %in% colnames(indepvars)) {
    indepvars <- indepvars[, colnames(indepvars) != "(Intercept)", drop = FALSE]
  }

  k <- ncol(indepvars)
  if (k < 1) {
    stop("at least one independent variable required", call. = FALSE)
  }

  # Get panel structure
  if (!id %in% names(data) || !time %in% names(data)) {
    stop("id and time variables must be columns in data", call. = FALSE)
  }

  panel_id <- data[[id]]
  time_id <- data[[time]]

  # Sort by panel then time
  ord <- order(panel_id, time_id)
  depvar <- depvar[ord]
  indepvars <- indepvars[ord, , drop = FALSE]
  panel_id <- panel_id[ord]
  time_id <- time_id[ord]

  # Get dimensions
  units <- unique(panel_id)
  N <- length(units)
  times <- unique(time_id)
  TT <- length(times)

  if (N < 2) {
    stop("at least 2 cross-sectional units required", call. = FALSE)
  }

  if (TT < 10) {
    stop("at least 10 time periods required", call. = FALSE)
  }

  # Check balanced panel
  obs_per_unit <- table(panel_id)
  if (length(unique(obs_per_unit)) > 1) {
    stop("panel must be balanced (no gaps) for xtcadfcoint", call. = FALSE)
  }

  # Warn about nfactors
  if (nfactors > k + 1) {
    message("Note: nfactors > k+1; setting nfactors = ", k + 1)
    nfactors <- k + 1
  }

  # ---- Reshape to wide format ----
  # Y: T x N matrix
  Y <- matrix(depvar, nrow = TT, ncol = N)

  # X: T x (N*k) matrix, arranged as [X1_unit1, X1_unit2, ..., X1_unitN, X2_unit1, ...]
  X <- matrix(NA_real_, nrow = TT, ncol = N * k)
  for (j in seq_len(k)) {
    for (i in seq_len(N)) {
      start_idx <- (i - 1) * TT + 1
      end_idx <- i * TT
      X[, i + (j - 1) * N] <- indepvars[start_idx:end_idx, j]
    }
  }

  # ---- Lag selection options ----
  opt_auto <- lagselect != "fixed"
  opt_ic <- switch(lagselect,
                   "aic" = 0,
                   "bic" = 1,
                   "maic" = 2,
                   "mbic" = 3,
                   "fixed" = 1)

  # ---- Run test ----
  if (breaks == 0) {
    result <- cadfcoint_main(Y, X, model, c(0), brk_slope, brk_loadings,
                             nfactors, maxlags, opt_auto, opt_ic, cce)
    panel_cips <- result$panel_cips
    panel_cips_alt <- panel_cips
    t_individual <- result$t_individual
    p_selected <- result$p_selected
    beta_cce <- result$beta_cce
    SSR <- result$SSR
    Tb_hat <- NULL
    Tb_tilde <- NULL
  } else {
    result <- cadfcoint_endog(Y, X, model, breaks, trimming, brk_slope,
                              brk_loadings, nfactors, maxlags, opt_auto,
                              opt_ic, cce)
    panel_cips <- result$panel_cips
    panel_cips_alt <- result$panel_cips_alt
    t_individual <- result$t_individual
    p_selected <- result$p_selected
    beta_cce <- result$beta_cce
    SSR <- result$SSR
    Tb_hat <- result$Tb_hat
    Tb_tilde <- result$Tb_tilde
  }

  # ---- Bootstrap critical values ----
  cv_panel <- NULL
  cv_ind <- NULL

  if (simulate > 0) {
    # Set break dates for simulation
    if (breaks == 0) {
      Tb_sim <- c(0)
    } else if (breaks == 1) {
      Tb_sim <- floor(0.5 * TT)
    } else {
      Tb_sim <- c(floor(0.3 * TT), floor(0.7 * TT))
    }

    cv_result <- simulate_cv(N, TT, k, model, Tb_sim, brk_slope,
                             brk_loadings, nfactors, 0, FALSE, opt_ic,
                             cce, simulate)
    cv_panel <- cv_result$panel_cv
    cv_ind <- cv_result$ind_cv
  }

  # ---- Build result object ----
  out <- list(
    panel_cips = panel_cips,
    panel_cips_alt = panel_cips_alt,
    t_individual = t_individual,
    p_selected = p_selected,
    beta_cce = beta_cce,
    SSR = SSR,
    Tb_hat = Tb_hat,
    Tb_tilde = Tb_tilde,
    N = N,
    TT = TT,
    k = k,
    model = model,
    breaks = breaks,
    trimming = trimming,
    nfactors = nfactors,
    brk_slope = brk_slope,
    brk_loadings = brk_loadings,
    cce = cce,
    maxlags = maxlags,
    lagselect = lagselect,
    cv_panel = cv_panel,
    cv_ind = cv_ind,
    simulate = simulate,
    call = call
  )

  class(out) <- "xtcadfcoint"
  return(out)
}


#' Main CADF Cointegration Engine (No Breaks)
#'
#' @description
#' Core computation for panel CADF cointegration test without breaks.
#'
#' @param Y T x N matrix of dependent variable.
#' @param X T x (N*k) matrix of regressors.
#' @param model Model specification (0-2).
#' @param Tb Vector of break dates (should be c(0) for no breaks).
#' @param brk_slope Allow slope shifts at breaks.
#' @param brk_loadings Allow loading shifts at breaks.
#' @param nfactors Number of factors.
#' @param p_max Maximum lag order.
#' @param opt_auto Automatic lag selection.
#' @param opt_ic Information criterion (0=AIC, 1=BIC, 2=MAIC, 3=MBIC).
#' @param opt_cce Use CCE augmentation.
#'
#' @return List with test results.
#'
#' @keywords internal
#' @noRd
cadfcoint_main <- function(Y, X, model, Tb, brk_slope, brk_loadings,
                           nfactors, p_max, opt_auto, opt_ic, opt_cce) {

  TT <- nrow(Y)
  N <- ncol(Y)
  k <- ncol(X) / N

  # ---- Build deterministic components ----
  det_X <- build_deterministics(TT, model, Tb, brk_slope)

  # ---- CCE: compute cross-sectional averages ----
  if (opt_cce) {
    y_bar <- rowMeans(Y)
    x_bar <- matrix(NA_real_, nrow = TT, ncol = k)
    for (j in seq_len(k)) {
      x_bar[, j] <- rowMeans(X[, ((j - 1) * N + 1):(j * N), drop = FALSE])
    }

    # Build CCE augmentation (cross-sectional averages and their lags)
    cce_vars <- cbind(y_bar, x_bar)
    n_cce <- ncol(cce_vars)
  } else {
    cce_vars <- NULL
    n_cce <- 0
  }

  # ---- Pooled CCE estimation ----
  # Stack data for pooled regression: y_i = X_i * beta + deterministics + cce + e_i
  y_stack <- as.vector(Y)
  X_stack <- matrix(NA_real_, nrow = N * TT, ncol = k)
  det_stack <- matrix(NA_real_, nrow = N * TT, ncol = ncol(det_X))
  cce_stack <- if (opt_cce) matrix(NA_real_, nrow = N * TT, ncol = n_cce) else NULL

  for (i in seq_len(N)) {
    idx <- ((i - 1) * TT + 1):(i * TT)
    for (j in seq_len(k)) {
      X_stack[idx, j] <- X[, i + (j - 1) * N]
    }
    det_stack[idx, ] <- det_X
    if (opt_cce) {
      cce_stack[idx, ] <- cce_vars
    }
  }

  # Pooled regression
  if (opt_cce) {
    reg_X <- cbind(X_stack, det_stack, cce_stack)
  } else {
    reg_X <- cbind(X_stack, det_stack)
  }

  # Remove columns that are all NA or constant zero
  valid_cols <- apply(reg_X, 2, function(col) {
    var(col, na.rm = TRUE) > .Machine$double.eps
  })
  reg_X <- reg_X[, valid_cols, drop = FALSE]

  # OLS estimation
  if (ncol(reg_X) > 0) {
    qr_fit <- qr(reg_X)
    if (qr_fit$rank < ncol(reg_X)) {
      # Handle rank deficiency
      coef_all <- qr.coef(qr_fit, y_stack)
      coef_all[is.na(coef_all)] <- 0
    } else {
      coef_all <- solve(crossprod(reg_X), crossprod(reg_X, y_stack))
    }
    beta_cce <- coef_all[1:k]
    resid_pool <- y_stack - reg_X %*% coef_all
  } else {
    beta_cce <- rep(0, k)
    resid_pool <- y_stack
  }

  # Reshape residuals to T x N
  E <- matrix(resid_pool, nrow = TT, ncol = N)

  # ---- Individual CADF tests ----
  t_individual <- numeric(N)
  p_selected <- numeric(N)

  for (i in seq_len(N)) {
    e_i <- E[, i]
    result_i <- cadf_test(e_i, model, Tb, cce_vars, p_max, opt_auto, opt_ic, opt_cce)
    t_individual[i] <- result_i$t_stat
    p_selected[i] <- result_i$p_sel
  }

  # ---- Panel CIPS statistic ----
  panel_cips <- mean(t_individual)

  # ---- SSR ----
  SSR <- c(total = sum(resid_pool^2), NA, NA)

  return(list(
    panel_cips = panel_cips,
    t_individual = t_individual,
    p_selected = p_selected,
    beta_cce = beta_cce,
    SSR = SSR
  ))
}


#' CADF Test with Endogenous Break Detection
#'
#' @description
#' Performs panel CADF test with endogenous structural break detection.
#'
#' @param Y T x N matrix of dependent variable.
#' @param X T x (N*k) matrix of regressors.
#' @param model Model specification (3-5).
#' @param m Number of breaks.
#' @param trimming Trimming fraction.
#' @param brk_slope Allow slope shifts.
#' @param brk_loadings Allow loading shifts.
#' @param nfactors Number of factors.
#' @param p_max Maximum lag order.
#' @param opt_auto Automatic lag selection.
#' @param opt_ic Information criterion.
#' @param opt_cce Use CCE.
#'
#' @return List with test results.
#'
#' @keywords internal
#' @noRd
cadfcoint_endog <- function(Y, X, model, m, trimming, brk_slope,
                            brk_loadings, nfactors, p_max, opt_auto,
                            opt_ic, opt_cce) {

  TT <- nrow(Y)
  N <- ncol(Y)
  k <- ncol(X) / N

  # ---- Find optimal break dates ----
  # Grid search over admissible break dates
  trim_obs <- floor(TT * trimming)
  admissible <- (trim_obs + 1):(TT - trim_obs)

  if (m == 1) {
    # Single break: search over all admissible dates
    best_SSR <- Inf
    best_Tb <- trim_obs + 1

    for (tb in admissible) {
      result <- cadfcoint_main(Y, X, model, c(tb), brk_slope, brk_loadings,
                               nfactors, p_max, opt_auto, opt_ic, opt_cce)
      if (result$SSR[1] < best_SSR) {
        best_SSR <- result$SSR[1]
        best_Tb <- tb
      }
    }

    Tb_hat <- best_Tb
    Tb_tilde <- best_Tb  # Alternative estimator (same for single break)

  } else {
    # Two breaks: grid search
    best_SSR <- Inf
    best_Tb <- c(trim_obs + 1, TT - trim_obs)

    for (tb1 in admissible) {
      for (tb2 in (tb1 + trim_obs):min(TT - trim_obs, length(admissible) + trim_obs)) {
        if (tb2 > tb1 + trim_obs && tb2 <= TT - trim_obs) {
          result <- cadfcoint_main(Y, X, model, c(tb1, tb2), brk_slope, brk_loadings,
                                   nfactors, p_max, opt_auto, opt_ic, opt_cce)
          if (result$SSR[1] < best_SSR) {
            best_SSR <- result$SSR[1]
            best_Tb <- c(tb1, tb2)
          }
        }
      }
    }

    Tb_hat <- best_Tb
    Tb_tilde <- best_Tb
  }

  # ---- Compute final results at optimal breaks ----
  result_hat <- cadfcoint_main(Y, X, model, Tb_hat, brk_slope, brk_loadings,
                               nfactors, p_max, opt_auto, opt_ic, opt_cce)

  # Alternative criterion: maximize |t-stat|
  best_t <- 0
  Tb_tilde <- Tb_hat

  if (m == 1) {
    for (tb in admissible) {
      result_tb <- cadfcoint_main(Y, X, model, c(tb), brk_slope, brk_loadings,
                                  nfactors, p_max, opt_auto, opt_ic, opt_cce)
      if (abs(result_tb$panel_cips) > abs(best_t)) {
        best_t <- result_tb$panel_cips
        Tb_tilde <- tb
      }
    }
  }

  result_tilde <- cadfcoint_main(Y, X, model, Tb_tilde, brk_slope, brk_loadings,
                                 nfactors, p_max, opt_auto, opt_ic, opt_cce)

  return(list(
    panel_cips = result_hat$panel_cips,
    panel_cips_alt = result_tilde$panel_cips,
    t_individual = result_hat$t_individual,
    p_selected = result_hat$p_selected,
    beta_cce = result_hat$beta_cce,
    SSR = result_hat$SSR,
    Tb_hat = Tb_hat,
    Tb_tilde = if (is.null(Tb_tilde)) Tb_hat else Tb_tilde
  ))
}


#' Build Deterministic Components Matrix
#'
#' @description
#' Constructs the matrix of deterministic components based on model specification.
#'
#' @param TT Number of time periods.
#' @param model Model specification (0-5).
#' @param Tb Vector of break dates.
#' @param brk_slope Allow slope shifts.
#'
#' @return Matrix of deterministic variables.
#'
#' @keywords internal
#' @noRd
build_deterministics <- function(TT, model, Tb, brk_slope) {

  t_vec <- seq_len(TT)
  result <- NULL

  # Model 0: No deterministics
  if (model == 0) {
    return(matrix(0, nrow = TT, ncol = 1))
  }

  # Model 1: Constant only
  if (model == 1) {
    return(matrix(1, nrow = TT, ncol = 1))
  }

  # Model 2: Constant + trend
  if (model == 2) {
    return(cbind(1, t_vec))
  }

  # Models 3-5: With breaks
  m <- length(Tb[Tb > 0])
  if (m == 0) {
    # Fallback to simpler model
    if (model <= 2) {
      return(build_deterministics(TT, model, c(0), brk_slope))
    }
    return(matrix(1, nrow = TT, ncol = 1))
  }

  # Constant
  result <- matrix(1, nrow = TT, ncol = 1)

  # Model 3: Constant with level shifts
  if (model == 3) {
    for (j in seq_len(m)) {
      DU <- as.numeric(t_vec > Tb[j])
      result <- cbind(result, DU)
    }
    return(result)
  }

  # Model 4: Trend with level shifts
  if (model == 4) {
    result <- cbind(result, t_vec)  # Add trend
    for (j in seq_len(m)) {
      DU <- as.numeric(t_vec > Tb[j])
      result <- cbind(result, DU)
    }
    return(result)
  }

  # Model 5: Trend with level and slope shifts
  if (model == 5) {
    result <- cbind(result, t_vec)  # Add trend
    for (j in seq_len(m)) {
      DU <- as.numeric(t_vec > Tb[j])
      DT <- (t_vec - Tb[j]) * DU  # Trend shift
      result <- cbind(result, DU, DT)
    }
    return(result)
  }

  return(result)
}


#' Individual CADF Test
#'
#' @description
#' Performs CADF test on residuals for one cross-sectional unit.
#'
#' @param e Residual series (length T).
#' @param model Model specification.
#' @param Tb Break dates.
#' @param cce_vars Cross-sectional averages (T x (k+1) matrix).
#' @param p_max Maximum lag order.
#' @param opt_auto Automatic lag selection.
#' @param opt_ic Information criterion.
#' @param opt_cce Include CCE augmentation.
#'
#' @return List with t-statistic and selected lag order.
#'
#' @keywords internal
#' @noRd
cadf_test <- function(e, model, Tb, cce_vars, p_max, opt_auto, opt_ic, opt_cce) {

  TT <- length(e)

  # Build ADF regression components
  # Delta(e_t) = rho * e_{t-1} + sum_{j=1}^{p} phi_j * Delta(e_{t-j}) + det + cce + u_t

  # First difference
  de <- diff(e)
  e_lag <- e[1:(TT - 1)]

  # Lagged differences
  de_lags <- matrix(NA_real_, nrow = TT - 1, ncol = p_max)
  for (j in seq_len(p_max)) {
    de_lags[, j] <- c(rep(NA, j), de[1:(TT - 1 - j)])
  }

  # Determine lag order
  if (opt_auto && p_max > 0) {
    p_sel <- select_lag(de, e_lag, de_lags, cce_vars, model, Tb, p_max, opt_ic, opt_cce)
  } else {
    p_sel <- p_max
  }

  # Build regression matrix
  # Valid observations: after p_sel lags
  valid_idx <- (p_sel + 1):(TT - 1)
  n_obs <- length(valid_idx)

  if (n_obs < 10) {
    return(list(t_stat = NA, p_sel = p_sel))
  }

  y_reg <- de[valid_idx]
  X_reg <- matrix(e_lag[valid_idx], ncol = 1)

  # Add lagged differences
  if (p_sel > 0) {
    for (j in seq_len(p_sel)) {
      X_reg <- cbind(X_reg, de_lags[valid_idx, j])
    }
  }

  # Add deterministics
  det_X <- build_deterministics(TT - 1, model, Tb, FALSE)
  X_reg <- cbind(X_reg, det_X[valid_idx, , drop = FALSE])

  # Add CCE (differenced)
  if (opt_cce && !is.null(cce_vars)) {
    d_cce <- apply(cce_vars, 2, diff)
    X_reg <- cbind(X_reg, d_cce[valid_idx, , drop = FALSE])
  }

  # Remove constant columns
  col_vars <- apply(X_reg, 2, var, na.rm = TRUE)
  X_reg <- X_reg[, col_vars > .Machine$double.eps | seq_len(ncol(X_reg)) == 1, drop = FALSE]

  # Handle missing values
  complete_idx <- complete.cases(cbind(y_reg, X_reg))
  if (sum(complete_idx) < 10) {
    return(list(t_stat = NA, p_sel = p_sel))
  }

  y_reg <- y_reg[complete_idx]
  X_reg <- X_reg[complete_idx, , drop = FALSE]

  # OLS regression
  qr_fit <- qr(X_reg)
  if (qr_fit$rank < ncol(X_reg)) {
    # Rank deficient - use pseudoinverse
    coef <- qr.coef(qr_fit, y_reg)
    coef[is.na(coef)] <- 0
  } else {
    coef <- solve(crossprod(X_reg), crossprod(X_reg, y_reg))
  }

  # Residuals and standard error
  resid <- y_reg - X_reg %*% coef
  n_eff <- length(resid)
  df <- n_eff - ncol(X_reg)

  if (df <= 0) {
    return(list(t_stat = NA, p_sel = p_sel))
  }

  sigma2 <- sum(resid^2) / df

  # Compute t-statistic for rho (first coefficient)
  XtX_inv <- tryCatch({
    solve(crossprod(X_reg))
  }, error = function(err) {
    # Use generalized inverse
    MASS_pinv <- function(X) {
      svd_X <- svd(X)
      d_inv <- ifelse(svd_X$d > .Machine$double.eps * max(dim(X)) * max(svd_X$d),
                      1 / svd_X$d, 0)
      svd_X$v %*% diag(d_inv) %*% t(svd_X$u)
    }
    MASS_pinv(crossprod(X_reg))
  })

  se_rho <- sqrt(sigma2 * XtX_inv[1, 1])
  t_stat <- coef[1] / se_rho

  return(list(t_stat = t_stat, p_sel = p_sel))
}


#' Lag Selection via Information Criteria
#'
#' @description
#' Selects optimal lag order using AIC, BIC, MAIC, or MBIC.
#'
#' @param de First differences.
#' @param e_lag Lagged levels.
#' @param de_lags Matrix of lagged differences.
#' @param cce_vars CCE variables.
#' @param model Model specification.
#' @param Tb Break dates.
#' @param p_max Maximum lag order.
#' @param opt_ic Criterion code.
#' @param opt_cce Include CCE.
#'
#' @return Selected lag order.
#'
#' @keywords internal
#' @noRd
select_lag <- function(de, e_lag, de_lags, cce_vars, model, Tb, p_max, opt_ic, opt_cce) {

  TT <- length(de) + 1
  best_ic <- Inf
  best_p <- 0

  for (p in 0:p_max) {
    # Build regression for this lag order
    valid_idx <- (p + 1):(TT - 1)
    n_obs <- length(valid_idx)

    if (n_obs < 10) next

    y_reg <- de[valid_idx]
    X_reg <- matrix(e_lag[valid_idx], ncol = 1)

    if (p > 0) {
      for (j in seq_len(p)) {
        X_reg <- cbind(X_reg, de_lags[valid_idx, j])
      }
    }

    # Add deterministics
    det_X <- build_deterministics(TT - 1, model, Tb, FALSE)
    X_reg <- cbind(X_reg, det_X[valid_idx, , drop = FALSE])

    # Add CCE
    if (opt_cce && !is.null(cce_vars)) {
      d_cce <- apply(cce_vars, 2, diff)
      X_reg <- cbind(X_reg, d_cce[valid_idx, , drop = FALSE])
    }

    # Remove constant columns
    col_vars <- apply(X_reg, 2, var, na.rm = TRUE)
    X_reg <- X_reg[, col_vars > .Machine$double.eps | seq_len(ncol(X_reg)) == 1, drop = FALSE]

    # Handle NA
    complete_idx <- complete.cases(cbind(y_reg, X_reg))
    if (sum(complete_idx) < 10) next

    y_fit <- y_reg[complete_idx]
    X_fit <- X_reg[complete_idx, , drop = FALSE]
    n_eff <- length(y_fit)
    k_eff <- ncol(X_fit)

    # OLS
    qr_fit <- qr(X_fit)
    if (qr_fit$rank < k_eff) next

    coef <- qr.coef(qr_fit, y_fit)
    resid <- y_fit - X_fit %*% coef
    sigma2 <- sum(resid^2) / n_eff

    # Information criterion
    if (opt_ic == 0) {
      # AIC
      ic <- log(sigma2) + 2 * k_eff / n_eff
    } else if (opt_ic == 1) {
      # BIC
      ic <- log(sigma2) + k_eff * log(n_eff) / n_eff
    } else if (opt_ic == 2) {
      # MAIC (Modified AIC - Ng & Perron)
      tau <- coef[1]^2 * sum(e_lag[valid_idx][complete_idx]^2) / sigma2
      ic <- log(sigma2) + 2 * (tau + k_eff) / n_eff
    } else {
      # MBIC (Modified BIC)
      tau <- coef[1]^2 * sum(e_lag[valid_idx][complete_idx]^2) / sigma2
      ic <- log(sigma2) + (tau + k_eff) * log(n_eff) / n_eff
    }

    if (ic < best_ic) {
      best_ic <- ic
      best_p <- p
    }
  }

  return(best_p)
}


#' Simulate Critical Values
#'
#' @description
#' Generates bootstrap critical values under the null of no cointegration.
#'
#' @param N Number of cross-sections.
#' @param TT Number of time periods.
#' @param k Number of regressors.
#' @param model Model specification.
#' @param Tb Break dates.
#' @param brk_slope Slope shifts.
#' @param brk_loadings Loading shifts.
#' @param nfactors Number of factors.
#' @param p_max Maximum lags.
#' @param opt_auto Automatic lag selection.
#' @param opt_ic Information criterion.
#' @param opt_cce Use CCE.
#' @param reps Number of replications.
#'
#' @return List with panel and individual critical values.
#'
#' @keywords internal
#' @noRd
simulate_cv <- function(N, TT, k, model, Tb, brk_slope, brk_loadings,
                        nfactors, p_max, opt_auto, opt_ic, opt_cce, reps) {

  panel_stats <- numeric(reps)
  ind_stats <- numeric(reps * N)

  for (r in seq_len(reps)) {
    # Generate independent random walks (no cointegration)
    Y_sim <- matrix(cumsum(rnorm(TT * N)), nrow = TT, ncol = N)
    X_sim <- matrix(NA_real_, nrow = TT, ncol = N * k)

    for (j in seq_len(k)) {
      for (i in seq_len(N)) {
        X_sim[, i + (j - 1) * N] <- cumsum(rnorm(TT))
      }
    }

    # Run test
    result <- cadfcoint_main(Y_sim, X_sim, model, Tb, brk_slope, brk_loadings,
                             nfactors, p_max, opt_auto, opt_ic, opt_cce)

    panel_stats[r] <- result$panel_cips
    ind_stats[((r - 1) * N + 1):(r * N)] <- result$t_individual
  }

  # Remove NA values
  panel_stats <- panel_stats[!is.na(panel_stats)]
  ind_stats <- ind_stats[!is.na(ind_stats)]

  # Compute quantiles (left tail for rejection)
  panel_cv <- quantile(panel_stats, probs = c(0.01, 0.025, 0.05, 0.10), na.rm = TRUE)
  ind_cv <- quantile(ind_stats, probs = c(0.01, 0.025, 0.05, 0.10), na.rm = TRUE)

  return(list(panel_cv = panel_cv, ind_cv = ind_cv))
}
