test_that("gotham_state matches a hand-calculated balanced example", {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = 1:9
  )

  out <- gotham_state(
    dat,
    outcome = "y",
    group = "group"
  )

  expect_equal(out$n_groups, 3L)
  expect_equal(out$group_size, 3L)
  expect_equal(out$grand_mean, 5)
  expect_equal(out$group_means, c(A = 2, B = 5, C = 8))
  expect_equal(out$group_offsets, c(A = -3, B = 0, C = 3))
  expect_equal(out$group_offset_ss, 18)
  expect_equal(out$ss_between, 54)
  expect_equal(out$ss_within, 6)
  expect_equal(out$ms_between, 27)
  expect_equal(out$ms_within, 1)
  expect_equal(out$icc, 26 / 29, tolerance = 1e-12)
})


test_that("gotham_state obeys the balanced between-group identity", {
  dat <- data.frame(
    school = factor(rep(c("north", "south", "west", "east"), each = 5)),
    score = c(
      2, 3, 4, 5, 6,
      4, 5, 5, 6, 7,
      8, 8, 9, 10, 10,
      10, 11, 12, 12, 13
    )
  )

  out <- gotham_state(
    dat,
    outcome = "score",
    group = "school"
  )

  expect_equal(
    out$ss_between,
    out$group_size * out$group_offset_ss,
    tolerance = 1e-12
  )
  expect_equal(
    sum(out$group_offsets),
    0,
    tolerance = 1e-12
  )
})


test_that("gotham_state agrees with stats::aov and anova_icc", {
  dat <- data.frame(
    group = factor(rep(letters[1:4], each = 5)),
    y = c(
      2, 3, 4, 5, 6,
      4, 5, 5, 6, 7,
      8, 8, 9, 10, 10,
      10, 11, 12, 12, 13
    )
  )

  out <- gotham_state(
    dat,
    outcome = "y",
    group = "group"
  )

  fit <- stats::aov(y ~ group, data = dat)
  tab <- summary(fit)[[1L]]

  expect_equal(
    out$ss_between,
    unname(tab["group", "Sum Sq"]),
    tolerance = 1e-12
  )
  expect_equal(
    out$ss_within,
    unname(tab["Residuals", "Sum Sq"]),
    tolerance = 1e-12
  )
  expect_equal(
    out$ms_between,
    unname(tab["group", "Mean Sq"]),
    tolerance = 1e-12
  )
  expect_equal(
    out$ms_within,
    unname(tab["Residuals", "Mean Sq"]),
    tolerance = 1e-12
  )
  expect_equal(
    out$icc,
    anova_icc(dat, outcome = "y", group = "group"),
    tolerance = 1e-12
  )
})


test_that("gotham_state preserves factor group order", {
  dat <- data.frame(
    group = factor(
      rep(c("B", "A", "C"), each = 3),
      levels = c("B", "A", "C")
    ),
    y = 1:9
  )

  out <- gotham_state(
    dat,
    outcome = "y",
    group = "group"
  )

  expect_identical(names(out$group_means), c("B", "A", "C"))
  expect_equal(out$group_means, c(B = 2, A = 5, C = 8))
})


test_that("gotham_state uses the same input restrictions as anova_icc", {
  unbalanced <- data.frame(
    group = c("A", "A", "A", "B", "B"),
    y = 1:5
  )

  expect_error(
    gotham_state(unbalanced, outcome = "y", group = "group"),
    "balanced design"
  )

  missing_y <- data.frame(
    group = rep(c("A", "B"), each = 3),
    y = c(1, 2, NA, 4, 5, 6)
  )

  expect_error(
    gotham_state(missing_y, outcome = "y", group = "group"),
    "finite, non-missing"
  )
})
