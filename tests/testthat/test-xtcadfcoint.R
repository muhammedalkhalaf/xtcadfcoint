test_that("xtcadfcoint returns correct structure (no breaks)", {
  set.seed(123)
  n <- 6; tt <- 30
  dat <- data.frame(
    id   = rep(1:n, each = tt),
    time = rep(1:tt, times = n),
    y    = c(apply(matrix(rnorm(n * tt), tt, n), 2, cumsum)),
    x1   = c(apply(matrix(rnorm(n * tt), tt, n), 2, cumsum))
  )
  res <- xtcadfcoint(y ~ x1, data = dat, index = c("id", "time"),
                     model = 1, breaks = 0)
  expect_s3_class(res, "xtcadfcoint")
  expect_equal(res$N, n)
  expect_equal(res$TT, tt)
  expect_equal(res$k, 1)
  expect_length(res$t_individual, n)
  expect_false(is.na(res$panel_cips))
})

test_that("xtcadfcoint model 2 works", {
  set.seed(42)
  n <- 5; tt <- 25
  dat <- data.frame(
    id   = rep(1:n, each = tt),
    time = rep(1:tt, times = n),
    y    = cumsum(rnorm(n * tt)),
    x1   = cumsum(rnorm(n * tt))
  )
  res <- xtcadfcoint(y ~ x1, data = dat, index = c("id", "time"),
                     model = 2, breaks = 0)
  expect_s3_class(res, "xtcadfcoint")
  expect_equal(res$model, 2)
})

test_that("xtcadfcoint rejects wrong model/break combos", {
  set.seed(1)
  n <- 4; tt <- 20
  dat <- data.frame(
    id   = rep(1:n, each = tt),
    time = rep(1:tt, times = n),
    y    = cumsum(rnorm(n * tt)),
    x1   = cumsum(rnorm(n * tt))
  )
  expect_error(
    xtcadfcoint(y ~ x1, data = dat, index = c("id", "time"),
                model = 3, breaks = 0),
    "breaks"
  )
  expect_error(
    xtcadfcoint(y ~ x1, data = dat, index = c("id", "time"),
                model = 1, breaks = 1),
    "breaks"
  )
})

test_that("print and summary do not error", {
  set.seed(7)
  n <- 4; tt <- 20
  dat <- data.frame(
    id   = rep(1:n, each = tt),
    time = rep(1:tt, times = n),
    y    = cumsum(rnorm(n * tt)),
    x1   = cumsum(rnorm(n * tt))
  )
  res <- xtcadfcoint(y ~ x1, data = dat, index = c("id", "time"),
                     model = 1, breaks = 0)
  expect_output(print(res))
  expect_output(summary(res))
})
