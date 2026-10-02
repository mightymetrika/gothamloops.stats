test_that("gotham_update matches brute-force data augmentation", {
  dat <- data.frame(
    group = factor(rep(c("A", "B", "C", "D"), each = 4)),
    y = c(
      1, 2, 3, 4,
      4, 5, 6, 7,
      8, 9, 10, 11,
      12, 13, 14, 15
    )
  )

  state <- gotham_state(dat, outcome = "y", group = "group")
  deviations <- c(A = -2, B = 1.5, C = -0.5, D = 3)

  out <- gotham_update(state, deviations)

  added <- data.frame(
    group = factor(
      names(deviations),
      levels = levels(dat$group)
    ),
    y = state$group_means[names(deviations)] + deviations
  )

  brute <- gotham_state(
    rbind(dat, added),
    outcome = "y",
    group = "group"
  )

  fields <- c(
    "n_groups",
    "group_size",
    "grand_mean",
    "group_means",
    "group_offsets",
    "group_offset_ss",
    "ss_between",
    "ss_within",
    "ms_between",
    "ms_within",
    "icc"
  )

  for (field in fields) {
    expect_equal(out[[field]], brute[[field]], tolerance = 1e-12)
  }
})


test_that("gotham_update obeys the analytic within-group update identity", {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = 1:9
  )

  state <- gotham_state(dat, outcome = "y", group = "group")
  deviations <- c(A = -1, B = 2, C = 4)

  out <- gotham_update(state, deviations)

  expected_ss_within <- state$ss_within +
    state$group_size / (state$group_size + 1) * sum(deviations^2)

  expect_equal(out$ss_within, expected_ss_within, tolerance = 1e-12)
  expect_equal(out$proposal$deviation_ss, sum(deviations^2), tolerance = 1e-12)
})


test_that("gotham_update obeys the analytic between-group update identity", {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = 1:9
  )

  state <- gotham_state(dat, outcome = "y", group = "group")
  deviations <- c(A = -1, B = 2, C = 4)
  centered <- deviations - mean(deviations)
  n_new <- state$group_size + 1

  expected_ss_between <-
    n_new * state$group_offset_ss +
    2 * sum(state$group_offsets * centered) +
    sum(centered^2) / n_new

  out <- gotham_update(state, deviations)

  expect_equal(out$ss_between, expected_ss_between, tolerance = 1e-12)
  expect_equal(
    out$proposal$centered_deviations,
    centered,
    tolerance = 1e-12
  )
})


test_that("gotham_update records the expected proposal geometry", {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = 1:9
  )

  state <- gotham_state(dat, outcome = "y", group = "group")

  reinforcing <- gotham_update(
    state,
    deviations = c(A = -2, B = 0, C = 2)
  )
  opposing <- gotham_update(
    state,
    deviations = c(A = 2, B = 0, C = -2)
  )
  common_shift <- gotham_update(
    state,
    deviations = c(A = 2, B = 2, C = 2)
  )

  expect_equal(reinforcing$proposal$alignment, 1, tolerance = 1e-12)
  expect_equal(opposing$proposal$alignment, -1, tolerance = 1e-12)
  expect_true(is.na(common_shift$proposal$alignment))
  expect_equal(common_shift$proposal$centered_deviation_ss, 0)
})


test_that("gotham_update reorders named deviations to the state group order", {
  dat <- data.frame(
    group = factor(
      rep(c("B", "A", "C"), each = 3),
      levels = c("B", "A", "C")
    ),
    y = 1:9
  )

  state <- gotham_state(dat, outcome = "y", group = "group")

  out <- gotham_update(
    state,
    deviations = c(C = 3, B = 1, A = 2)
  )

  expect_identical(names(out$proposal$deviations), c("B", "A", "C"))
  expect_equal(out$proposal$deviations, c(B = 1, A = 2, C = 3))
})


test_that("gotham_update supports sequential analytic updates", {
  dat <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 3)),
    y = 1:9
  )

  state0 <- gotham_state(dat, outcome = "y", group = "group")

  d1 <- c(A = -1, B = 0.5, C = 2)
  state1 <- gotham_update(state0, d1)

  d2 <- c(A = 1.25, B = -0.75, C = 0.5)
  state2 <- gotham_update(state1, d2)

  added1 <- data.frame(
    group = factor(names(d1), levels = levels(dat$group)),
    y = state0$group_means[names(d1)] + d1
  )
  dat1 <- rbind(dat, added1)
  brute1 <- gotham_state(dat1, outcome = "y", group = "group")

  added2 <- data.frame(
    group = factor(names(d2), levels = levels(dat$group)),
    y = brute1$group_means[names(d2)] + d2
  )
  brute2 <- gotham_state(
    rbind(dat1, added2),
    outcome = "y",
    group = "group"
  )

  expect_equal(state2$group_size, 5L)
  expect_equal(state2$grand_mean, brute2$grand_mean, tolerance = 1e-12)
  expect_equal(state2$group_means, brute2$group_means, tolerance = 1e-12)
  expect_equal(state2$ss_between, brute2$ss_between, tolerance = 1e-12)
  expect_equal(state2$ss_within, brute2$ss_within, tolerance = 1e-12)
  expect_equal(state2$icc, brute2$icc, tolerance = 1e-12)
})


test_that("gotham_update rejects invalid proposals and invalid states", {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = 1:9
  )

  state <- gotham_state(dat, outcome = "y", group = "group")

  expect_error(
    gotham_update(state, c(1, 2)),
    "one value per group"
  )
  expect_error(
    gotham_update(state, c(A = 1, B = 2, D = 3)),
    "group labels exactly"
  )
  expect_error(
    gotham_update(state, c(A = 1, B = NA_real_, C = 3)),
    "finite, non-missing"
  )
  expect_error(
    gotham_update(list(), c(1, 2, 3)),
    "returned by gotham_state"
  )
})
