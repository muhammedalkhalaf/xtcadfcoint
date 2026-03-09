test_that("xtcadfcoint runs with basic input", {
  set.seed(123)
  N <- 10
  TT <- 50

  panel_data <- data.frame(
    id = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y = rnorm(N * TT),
    x = rnorm(N * TT)
  )

  result <- xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time")

  expect_s3_class(result, "xtcadfcoint")
  expect_equal(result$N, N)
  expect_equal(result$TT, TT)
  expect_equal(result$k, 1)
  expect_length(result$t_individual, N)
  expect_length(result$p_selected, N)
})


test_that("model specification validation works", {
  set.seed(456)
  N <- 5
  TT <- 30

  panel_data <- data.frame(
    id = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y = rnorm(N * TT),
    x = rnorm(N * TT)
  )

  # Models 0-2 should work without breaks
  for (m in 0:2) {
    result <- xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time",
                          model = m, breaks = 0)
    expect_equal(result$model, m)
  }

  # Model 3+ requires breaks
  expect_error(
    xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time",
                model = 3, breaks = 0),
    "models 3-5 require breaks >= 1"
  )

  # Breaks > 0 requires model >= 3
  expect_error(
    xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time",
                model = 1, breaks = 1),
    "breaks > 0 requires model >= 3"
  )
})


test_that("CCE option works", {
  set.seed(789)
  N <- 8
  TT <- 40

  panel_data <- data.frame(
    id = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y = rnorm(N * TT),
    x = rnorm(N * TT)
  )

  # With CCE
  result_cce <- xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time",
                            cce = TRUE)
  expect_true(result_cce$cce)

  # Without CCE
  result_no_cce <- xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time",
                               cce = FALSE)
  expect_false(result_no_cce$cce)

  # Results should differ
  expect_false(isTRUE(all.equal(result_cce$panel_cips, result_no_cce$panel_cips)))
})


test_that("lag selection methods work", {
  set.seed(111)
  N <- 6
  TT <- 35

  panel_data <- data.frame(
    id = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y = rnorm(N * TT),
    x = rnorm(N * TT)
  )

  methods <- c("aic", "bic", "maic", "mbic", "fixed")

  for (method in methods) {
    result <- xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time",
                          lagselect = method, maxlags = 2)
    expect_equal(result$lagselect, method)
  }
})


test_that("multiple regressors work", {
  set.seed(222)
  N <- 8
  TT <- 40

  panel_data <- data.frame(
    id = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y = rnorm(N * TT),
    x1 = rnorm(N * TT),
    x2 = rnorm(N * TT)
  )

  result <- xtcadfcoint(y ~ x1 + x2, data = panel_data, id = "id", time = "time")

  expect_equal(result$k, 2)
  expect_length(result$beta_cce, 2)
})


test_that("structural breaks are estimated", {
  set.seed(333)
  N <- 10
  TT <- 60

  panel_data <- data.frame(
    id = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y = rnorm(N * TT),
    x = rnorm(N * TT)
  )

  # One break
  result1 <- xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time",
                         model = 3, breaks = 1)

  expect_equal(result1$breaks, 1)
  expect_length(result1$Tb_hat, 1)
  expect_true(result1$Tb_hat > 0 && result1$Tb_hat < TT)

  # Two breaks
  result2 <- xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time",
                         model = 3, breaks = 2, trimming = 0.15)

  expect_equal(result2$breaks, 2)
  expect_length(result2$Tb_hat, 2)
  expect_true(all(result2$Tb_hat > 0 & result2$Tb_hat < TT))
})


test_that("print method works", {
  set.seed(444)
  N <- 5
  TT <- 30

  panel_data <- data.frame(
    id = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y = rnorm(N * TT),
    x = rnorm(N * TT)
  )

  result <- xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time")

  # Capture output
  output <- capture.output(print(result))
  expect_true(length(output) > 0)
  expect_true(any(grepl("CIPS", output)))
})


test_that("coef method returns coefficients", {
  set.seed(555)
  N <- 5
  TT <- 30

  panel_data <- data.frame(
    id = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y = rnorm(N * TT),
    x1 = rnorm(N * TT),
    x2 = rnorm(N * TT)
  )

  result <- xtcadfcoint(y ~ x1 + x2, data = panel_data, id = "id", time = "time")
  coeffs <- coef(result)

  expect_named(coeffs)
  expect_length(coeffs, 2)
})


test_that("input validation catches errors", {
  set.seed(666)
  N <- 3
  TT <- 20

  panel_data <- data.frame(
    id = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y = rnorm(N * TT),
    x = rnorm(N * TT)
  )

  # Invalid model
  expect_error(
    xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time", model = 6),
    "model must be between 0 and 5"
  )

  # Invalid breaks
  expect_error(
    xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time", breaks = 3),
    "breaks must be 0, 1, or 2"
  )

  # Invalid trimming
  expect_error(
    xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time", trimming = 0.6),
    "trimming must be between 0 and 0.5"
  )

  # Missing id/time columns
  expect_error(
    xtcadfcoint(y ~ x, data = panel_data, id = "missing", time = "time"),
    "id and time variables must be columns in data"
  )
})


test_that("simulation of critical values works", {
  skip_on_cran()  # Skip on CRAN due to time

  set.seed(777)
  N <- 5
  TT <- 25

  panel_data <- data.frame(
    id = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y = rnorm(N * TT),
    x = rnorm(N * TT)
  )

  result <- xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time",
                        simulate = 50)  # Small number for speed

  expect_true(!is.null(result$cv_panel))
  expect_length(result$cv_panel, 4)  # 1%, 2.5%, 5%, 10%
  expect_true(!is.null(result$cv_ind))
})


test_that("is_cointegrated function works", {
  skip_on_cran()

  set.seed(888)
  N <- 5
  TT <- 25

  panel_data <- data.frame(
    id = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y = rnorm(N * TT),
    x = rnorm(N * TT)
  )

  result <- xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time",
                        simulate = 50)

  # Should return logical
  coint <- is_cointegrated(result, level = 0.05)
  expect_type(coint, "logical")

  # Error without CVs
  result_no_cv <- xtcadfcoint(y ~ x, data = panel_data, id = "id", time = "time")
  expect_error(is_cointegrated(result_no_cv), "Bootstrap critical values not available")
})
