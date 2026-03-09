#' Fisher Effect Panel Data
#'
#' @description
#' Simulated panel data for testing the Fisher effect relationship between
#' nominal interest rates and inflation. Contains 10 countries over 50 time
#' periods, designed to exhibit cointegration under the full Fisher hypothesis.
#'
#' @format A data frame with 500 rows and 4 variables:
#' \describe{
#'   \item{country}{Integer. Country identifier (1 to 10).}
#'   \item{year}{Integer. Time period (1 to 50).}
#'   \item{interest}{Numeric. Nominal interest rate (percent).}
#'   \item{inflation}{Numeric. Inflation rate (percent).}
#' }
#'
#' @details
#' The data is simulated from a panel cointegration model with common factors:
#' \deqn{interest_{it} = \alpha_i + \beta \cdot inflation_{it} + \gamma_i f_t + e_{it}}
#' where \eqn{f_t} is a common factor, \eqn{\gamma_i} are heterogeneous
#' factor loadings, and \eqn{e_{it}} is a stationary error.
#'
#' The data exhibits:
#' \itemize{
#'   \item Cross-sectional dependence through the common factor
#'   \item Heterogeneous intercepts across countries
#'   \item Cointegration between interest and inflation
#' }
#'
#' @examples
#' data(fisher_panel)
#' head(fisher_panel)
#'
#' # Test for cointegration
#' \donttest{
#' result <- xtcadfcoint(interest ~ inflation, data = fisher_panel,
#'                      id = "country", time = "year", model = 1)
#' print(result)
#' }
#'
#' @references
#' Fisher, I. (1930). \emph{The Theory of Interest}. Macmillan.
#'
#' @source Simulated data for illustration purposes.
#'
#' @docType data
#' @keywords datasets
#' @name fisher_panel
#' @usage data(fisher_panel)
NULL
