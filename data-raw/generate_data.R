# Generate example dataset for xtcadfcoint package

set.seed(42)

N <- 10  # countries
TT <- 50  # time periods

# Common factor (global shock)
f_t <- cumsum(rnorm(TT, 0, 0.3))

# Generate cointegrated panel
fisher_panel <- data.frame(
  country = rep(1:N, each = TT),
  year = rep(1:TT, N)
)

# Heterogeneous parameters
alpha_i <- rnorm(N, 2, 0.5)  # intercepts
gamma_i <- rnorm(N, 1, 0.3)  # factor loadings
beta <- 1.0  # cointegrating coefficient (Fisher effect)

# Generate data
inflation <- numeric(N * TT)
interest <- numeric(N * TT)

for (i in 1:N) {
  idx <- ((i - 1) * TT + 1):(i * TT)

  # Inflation follows a random walk
  infl_i <- cumsum(rnorm(TT, 0.02, 0.5))
  inflation[idx] <- infl_i

  # Interest rate is cointegrated with inflation
  # interest = alpha_i + beta * inflation + gamma_i * f_t + stationary error
  e_it <- arima.sim(list(ar = 0.5), n = TT, sd = 0.3)
  interest[idx] <- alpha_i[i] + beta * infl_i + gamma_i[i] * f_t + e_it
}

fisher_panel$inflation <- inflation
fisher_panel$interest <- interest

# Save data
usethis::use_data(fisher_panel, overwrite = TRUE)
