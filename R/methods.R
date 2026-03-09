#' @export
print.xtcadfcoint <- function(x, digits = 4, ...) {

  cat("\n")
  cat(rep("=", 70), "\n", sep = "")
  cat(" Banerjee & Carrion-i-Silvestre Panel CADF Cointegration Test\n")
  cat(rep("=", 70), "\n", sep = "")
  cat("\n")

  cat("  H0: No cointegration\n")
  cat("  H1: Panel is cointegrated\n")
  cat("\n")

  # Model description
  model_desc <- switch(as.character(x$model),
                       "0" = "No deterministic components",
                       "1" = "Constant",
                       "2" = "Constant and linear trend",
                       "3" = "Constant with level shifts",
                       "4" = "Linear trend with level shifts",
                       "5" = "Linear trend with level and slope shifts")

  cat(rep("-", 70), "\n", sep = "")
  cat("  Model specification     :", model_desc, "\n")
  cat("  Number of breaks (m)    :", x$breaks, "\n")
  cat("  Cross-section dep (CCE) :", ifelse(x$cce, "Yes", "No"), "\n")
  cat("  Number of factors       :", x$nfactors, "\n")
  cat("  Cointegrating vec shift :", ifelse(x$brk_slope, "Yes", "No"), "\n")
  cat("  Factor loading shift    :", ifelse(x$brk_loadings, "Yes", "No"), "\n")
  cat("  Lag selection           :", x$lagselect, "(max =", x$maxlags, ")\n")
  cat(rep("-", 70), "\n", sep = "")
  cat("  Panel dimensions        : N =", x$N, ", T =", x$TT, "\n")
  cat("  Number of regressors    : k =", x$k, "\n")
  cat(rep("-", 70), "\n", sep = "")
  cat("\n")

  # Test statistics
  cat(" Panel CIPS Cointegration Test Results\n")
  cat(rep("-", 70), "\n", sep = "")

  if (x$breaks == 0) {
    cat("  CIPS statistic          :", format(x$panel_cips, digits = digits), "\n")
  } else {
    cat("  CIPS statistic (SSR)    :", format(x$panel_cips, digits = digits), "\n")
    cat("  CIPS statistic (|t|)    :", format(x$panel_cips_alt, digits = digits), "\n")
  }

  cat(rep("-", 70), "\n", sep = "")

  # Critical values if available
  if (!is.null(x$cv_panel)) {
    cat("\n")
    cat(" Simulated Critical Values (", x$simulate, " replications)\n", sep = "")
    cat(rep("-", 70), "\n", sep = "")
    cat("  Significance     Critical Value     Test Stat     Decision\n")
    cat(rep("-", 70), "\n", sep = "")

    levels <- c(1, 2.5, 5, 10)
    cv_names <- c("1%", "2.5%", "5%", "10%")

    for (i in seq_along(levels)) {
      cv <- x$cv_panel[i]
      decision <- if (x$panel_cips < cv) "Reject H0" else "Fail to reject"
      stars <- if (x$panel_cips < cv) {
        if (i == 1) " ***" else if (i == 2) " ***" else if (i == 3) " **" else " *"
      } else ""

      cat(sprintf("  %5s          %12.4f       %12.4f    %s%s\n",
                  cv_names[i], cv, x$panel_cips, decision, stars))
    }
    cat(rep("-", 70), "\n", sep = "")
  } else {
    cat("\n  Note: Critical values depend on N, T, k, model specification.\n")
    cat("  Use simulate argument to generate bootstrap critical values.\n")
  }

  # Break dates
  if (x$breaks > 0 && !is.null(x$Tb_hat)) {
    cat("\n")
    cat(" Estimated Break Dates\n")
    cat(rep("-", 50), "\n", sep = "")
    for (j in seq_along(x$Tb_hat)) {
      cat(sprintf("  Break %d: Tb_hat = %d, Tb_tilde = %d\n",
                  j, x$Tb_hat[j],
                  if (length(x$Tb_tilde) >= j) x$Tb_tilde[j] else x$Tb_hat[j]))
    }
    cat(rep("-", 50), "\n", sep = "")
  }

  cat("\n")
  cat(rep("=", 70), "\n", sep = "")
  cat("  Reference: Banerjee & Carrion-i-Silvestre (2017, JTSA)\n")
  cat(rep("=", 70), "\n", sep = "")
  cat("\n")

  invisible(x)
}


#' @export
summary.xtcadfcoint <- function(object, ...) {

  # Print main results
  print(object, ...)

  # Additional details
  cat("\n")
  cat(" Pooled CCE Coefficients (beta_hat)\n")
  cat(rep("-", 40), "\n", sep = "")

  if (!is.null(object$beta_cce)) {
    for (j in seq_along(object$beta_cce)) {
      cat(sprintf("  beta[%d]  = %12.6f\n", j, object$beta_cce[j]))
    }
  }
  cat(rep("-", 40), "\n", sep = "")

  # Individual statistics
  cat("\n")
  cat(" Individual CADF Statistics\n")
  cat(rep("-", 50), "\n", sep = "")
  cat("  Unit     t-statistic    Selected lag\n")
  cat(rep("-", 50), "\n", sep = "")

  for (i in seq_len(object$N)) {
    cat(sprintf("  %4d    %12.4f     %6d\n",
                i, object$t_individual[i], object$p_selected[i]))
  }
  cat(rep("-", 50), "\n", sep = "")

  cat("\n")

  invisible(object)
}


#' Extract Coefficients from xtcadfcoint Object
#'
#' @description
#' Extracts the pooled CCE coefficient estimates.
#'
#' @param object An object of class \code{"xtcadfcoint"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A named numeric vector of coefficient estimates.
#'
#' @export
coef.xtcadfcoint <- function(object, ...) {
  beta <- object$beta_cce
  names(beta) <- paste0("beta", seq_along(beta))
  return(beta)
}


#' Check Cointegration Rejection
#'
#' @description
#' Determines whether the null of no cointegration is rejected at a
#' given significance level.
#'
#' @param object An object of class \code{"xtcadfcoint"}.
#' @param level Significance level (e.g., 0.05 for 5\%). Default is 0.05.
#'
#' @return Logical indicating whether H0 is rejected (TRUE = cointegrated).
#'
#' @details
#' Requires that bootstrap critical values were simulated (simulate > 0).
#' The test rejects the null of no cointegration if the CIPS statistic
#' is below the critical value (left-tail test).
#'
#' @examples
#' \donttest{
#' # Example with simulated CVs
#' set.seed(123)
#' N <- 10; TT <- 50
#' panel_data <- data.frame(
#'   id = rep(1:N, each = TT),
#'   time = rep(1:TT, N),
#'   y = rnorm(N * TT),
#'   x = rnorm(N * TT)
#' )
#' result <- xtcadfcoint(y ~ x, data = panel_data,
#'                      id = "id", time = "time",
#'                      simulate = 100)
#' is_cointegrated(result, level = 0.05)
#' }
#'
#' @export
is_cointegrated <- function(object, level = 0.05) {

  if (!inherits(object, "xtcadfcoint")) {
    stop("object must be of class 'xtcadfcoint'", call. = FALSE)
  }

  if (is.null(object$cv_panel)) {
    stop("Bootstrap critical values not available. ",
         "Rerun with simulate > 0.", call. = FALSE)
  }

  # Map level to quantile index
  levels <- c(0.01, 0.025, 0.05, 0.10)
  idx <- which.min(abs(levels - level))

  # Reject if test stat < critical value
  return(object$panel_cips < object$cv_panel[idx])
}
