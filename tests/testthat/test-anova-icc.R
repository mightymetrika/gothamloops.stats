test_that("anova_icc matches a hand-calculated balanced example", {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = 1:9
  )

  # Group means are 2, 5, and 8; grand mean is 5.
  # SS_between = 54, MS_between = 27.
  # SS_within = 6, MS_within = 1.
  # ICC = (27 - 1) / (27 + 2 * 1) = 26 / 29.
  expect_equal(
    anova_icc(dat, outcome = "y", group = "group"),
    26 / 29,
    tolerance = 1e-12
  )
})


test_that("internal ANOVA components match known values", {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = 1:9
  )

  out <- .anova_icc_components(
    dat,
    outcome = "y",
    group = "group"
  )

  expect_equal(out$n_groups, 3L)
  expect_equal(out$group_size, 3L)
  expect_equal(out$grand_mean, 5)
  expect_equal(out$group_means, c(2, 5, 8))
  expect_equal(out$ss_between, 54)
  expect_equal(out$ss_within, 6)
  expect_equal(out$df_between, 2L)
  expect_equal(out$df_within, 6L)
  expect_equal(out$ms_between, 27)
  expect_equal(out$ms_within, 1)
})


test_that("anova_icc is invariant to location and positive scale", {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 4),
    y = c(
      1.0, 1.5, 2.0, 2.5,
      4.0, 4.5, 5.0, 5.5,
      7.0, 7.5, 8.0, 8.5
    )
  )

  original <- anova_icc(dat, outcome = "y", group = "group")

  shifted <- dat
  shifted$y <- shifted$y + 100

  scaled <- dat
  scaled$y <- 7 * scaled$y

  expect_equal(
    anova_icc(shifted, outcome = "y", group = "group"),
    original,
    tolerance = 1e-12
  )

  expect_equal(
    anova_icc(scaled, outcome = "y", group = "group"),
    original,
    tolerance = 1e-12
  )
})


test_that("anova_icc retains negative ANOVA ICC estimates", {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = rep(c(0, 1, 2), times = 3)
  )

  expect_equal(
    anova_icc(dat, outcome = "y", group = "group"),
    -0.5,
    tolerance = 1e-12
  )
})


test_that("anova_icc agrees with ANOVA mean squares", {
  dat <- data.frame(
    group = factor(rep(letters[1:4], each = 5)),
    y = c(
      2, 3, 4, 5, 6,
      4, 5, 5, 6, 7,
      8, 8, 9, 10, 10,
      10, 11, 12, 12, 13
    )
  )

  fit <- stats::aov(y ~ group, data = dat)
  tab <- summary(fit)[[1L]]
  ms_between <- tab["group", "Mean Sq"]
  ms_within <- tab["Residuals", "Mean Sq"]
  n <- 5

  expected <- (ms_between - ms_within) /
    (ms_between + (n - 1) * ms_within)

  expect_equal(
    anova_icc(dat, outcome = "y", group = "group"),
    unname(expected),
    tolerance = 1e-12
  )
})


test_that("anova_icc rejects unsupported or invalid inputs", {
  balanced <- data.frame(
    group = rep(c("A", "B"), each = 3),
    y = 1:6
  )

  unbalanced <- data.frame(
    group = c("A", "A", "A", "B", "B"),
    y = 1:5
  )

  expect_error(
    anova_icc(unbalanced, outcome = "y", group = "group"),
    "balanced design"
  )

  missing_y <- balanced
  missing_y$y[1] <- NA_real_
  expect_error(
    anova_icc(missing_y, outcome = "y", group = "group"),
    "finite, non-missing"
  )

  nonnumeric <- balanced
  nonnumeric$y <- as.character(nonnumeric$y)
  expect_error(
    anova_icc(nonnumeric, outcome = "y", group = "group"),
    "must be numeric"
  )

  one_group <- data.frame(group = "A", y = 1:4)
  expect_error(
    anova_icc(one_group, outcome = "y", group = "group"),
    "At least two groups"
  )

  singleton_groups <- data.frame(
    group = c("A", "B", "C"),
    y = c(1, 2, 3)
  )
  expect_error(
    anova_icc(singleton_groups, outcome = "y", group = "group"),
    "at least two observations"
  )

  constant <- data.frame(
    group = rep(c("A", "B"), each = 3),
    y = rep(1, 6)
  )
  expect_error(
    anova_icc(constant, outcome = "y", group = "group"),
    "undefined"
  )
})
