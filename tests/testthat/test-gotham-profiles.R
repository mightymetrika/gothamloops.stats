.make_profile_state <- function() {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = 1:9
  )

  gotham_state(dat, outcome = "y", group = "group")
}


test_that("gotham_profile reproduces requested centered geometry", {
  state <- .make_profile_state()

  for (current_alignment in c(-1, -0.5, 0, 0.5, 1)) {
    out <- gotham_profile(
      state,
      radius = 2.75,
      alignment = current_alignment
    )

    expect_equal(mean(out$deviations), 0, tolerance = 1e-12)
    expect_equal(
      sqrt(sum(out$deviations^2)),
      2.75,
      tolerance = 1e-10
    )
    expect_equal(
      out$achieved$alignment,
      current_alignment,
      tolerance = 1e-10
    )
    expect_equal(out$achieved$radius, 2.75, tolerance = 1e-10)
    expect_equal(
      out$proposed_observations,
      state$group_means + out$deviations,
      tolerance = 1e-12
    )
  }
})


test_that("gotham_profile mean deviation shifts all proposal deviations equally", {
  state <- .make_profile_state()

  centered <- gotham_profile(
    state,
    radius = 2,
    alignment = 0,
    mean_deviation = 0
  )

  shifted <- gotham_profile(
    state,
    radius = 2,
    alignment = 0,
    mean_deviation = 3
  )

  expect_equal(
    unname(shifted$deviations - centered$deviations),
    rep(3, state$n_groups),
    tolerance = 1e-12
  )
  expect_equal(shifted$achieved$radius, centered$achieved$radius, tolerance = 1e-12)
  expect_equal(
    shifted$achieved$alignment,
    centered$achieved$alignment,
    tolerance = 1e-12
  )
  expect_equal(shifted$achieved$mean_deviation, 3, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(shifted$achieved$icc, centered$achieved$icc)))
})


test_that("orientation changes microstate but not centered macrostate", {
  state <- .make_profile_state()

  positive <- gotham_profile(
    state,
    radius = 3,
    alignment = 0.25,
    orientation = 1
  )

  negative <- gotham_profile(
    state,
    radius = 3,
    alignment = 0.25,
    orientation = -1
  )

  expect_false(isTRUE(all.equal(positive$deviations, negative$deviations)))
  expect_equal(positive$achieved$radius, negative$achieved$radius, tolerance = 1e-12)
  expect_equal(
    positive$achieved$alignment,
    negative$achieved$alignment,
    tolerance = 1e-12
  )
  expect_equal(positive$achieved$icc, negative$achieved$icc, tolerance = 1e-12)
})


test_that("gotham_profile preserves state group order", {
  dat <- data.frame(
    group = factor(
      rep(c("B", "A", "C"), each = 3),
      levels = c("B", "A", "C")
    ),
    y = 1:9
  )
  state <- gotham_state(dat, outcome = "y", group = "group")

  out <- gotham_profile(state, radius = 2, alignment = 0)

  expect_identical(names(out$deviations), c("B", "A", "C"))
  expect_identical(names(out$proposed_observations), c("B", "A", "C"))
})


test_that("zero-radius profile is a common shift", {
  state <- .make_profile_state()

  out <- gotham_profile(
    state,
    radius = 0,
    alignment = 0,
    mean_deviation = 2
  )

  expect_equal(out$deviations, c(A = 2, B = 2, C = 2))
  expect_equal(out$achieved$radius, 0)
  expect_true(is.na(out$achieved$alignment))
})


test_that("gotham_profile_set hits threshold and epsilon targets exactly", {
  state <- .make_profile_state()

  out <- gotham_profile_set(
    state,
    threshold = 0.80,
    epsilon = 0.02,
    alignment = c(-1, 0, 1)
  )

  expect_equal(nrow(out$summary), 9L)
  expect_true(all(out$summary$status == "generated"))
  expect_equal(
    out$summary$achieved_icc,
    out$summary$target_icc,
    tolerance = 1e-10
  )
  expect_equal(
    out$summary$achieved_alignment,
    out$summary$alignment,
    tolerance = 1e-10
  )

  expect_equal(nrow(out$profiles), 27L)
  expect_identical(
    unique(out$summary$target),
    c(
      "threshold_plus_epsilon",
      "threshold",
      "threshold_minus_epsilon"
    )
  )
})


test_that("profile-set radii order naturally in the monotone example", {
  state <- .make_profile_state()

  out <- gotham_profile_set(
    state,
    threshold = 0.80,
    epsilon = 0.02,
    alignment = c(-1, 0, 1)
  )

  for (current_alignment in c(-1, 0, 1)) {
    rows <- out$summary[out$summary$alignment == current_alignment, ]

    expect_lt(
      rows$radius[rows$relation == "above"],
      rows$radius[rows$relation == "equal"]
    )
    expect_lt(
      rows$radius[rows$relation == "equal"],
      rows$radius[rows$relation == "below"]
    )
  }
})


test_that("generated profile rows reproduce their requested ICCs", {
  state <- .make_profile_state()

  out <- gotham_profile_set(
    state,
    threshold = 0.80,
    epsilon = 0.02,
    alignment = c(-1, 0, 1)
  )

  keys <- unique(out$profiles[c("target", "alignment")])

  for (i in seq_len(nrow(keys))) {
    rows <- out$profiles[
      out$profiles$target == keys$target[i] &
        out$profiles$alignment == keys$alignment[i],
      ,
      drop = FALSE
    ]

    deviations <- stats::setNames(rows$deviation, rows$group)
    updated <- gotham_update(state, deviations)

    expect_equal(
      updated$icc,
      unique(rows$target_icc),
      tolerance = 1e-10
    )
    expect_equal(
      updated$proposal$alignment,
      unique(rows$alignment),
      tolerance = 1e-10
    )
  }
})


test_that("gotham_profile_set supports exact-threshold profiles without epsilon", {
  state <- .make_profile_state()

  out <- gotham_profile_set(
    state,
    threshold = 0.80,
    epsilon = 0,
    alignment = c(-1, 0, 1)
  )

  expect_equal(nrow(out$summary), 3L)
  expect_true(all(out$summary$target == "threshold"))
  expect_equal(out$summary$achieved_icc, rep(0.80, 3), tolerance = 1e-10)
})


test_that("gotham_profile validates geometric constraints", {
  state <- .make_profile_state()

  expect_error(
    gotham_profile(state, radius = -1, alignment = 0),
    "non-negative"
  )
  expect_error(
    gotham_profile(state, radius = 1, alignment = 1.1),
    "interval \\[-1, 1\\]"
  )
  expect_error(
    gotham_profile(state, radius = 1, alignment = 0, orientation = 0),
    "either -1 or 1"
  )

  two_group <- data.frame(
    group = rep(c("A", "B"), each = 3),
    y = c(1, 2, 3, 5, 6, 7)
  )
  two_state <- gotham_state(two_group, outcome = "y", group = "group")

  expect_error(
    gotham_profile(two_state, radius = 1, alignment = 0),
    "With two groups"
  )
})


test_that("gotham_profile_set validates threshold band and root", {
  state <- .make_profile_state()

  expect_error(
    gotham_profile_set(state, threshold = 0.01, epsilon = 0.02),
    "must both lie"
  )
  expect_error(
    gotham_profile_set(state, threshold = 0.99, epsilon = 0.02),
    "must both lie"
  )
  expect_error(
    gotham_profile_set(state, threshold = 0.80, root = 3),
    "either 1 or 2"
  )
})
