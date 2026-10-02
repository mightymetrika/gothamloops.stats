.make_band_state <- function() {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = 1:9
  )

  gotham_state(dat, outcome = "y", group = "group")
}


test_that("gotham_band reproduces known frontier radii", {
  state <- .make_band_state()

  out <- gotham_band(
    state,
    threshold = 0.80,
    epsilon = 0.02,
    alignment = c(-1, 0, 1)
  )

  expect_equal(out$radius_above, c(2.778148, 3.962204, 5.650908), tolerance = 1e-6)
  expect_equal(out$radius_threshold, c(3.026088, 4.369986, 6.310713), tolerance = 1e-6)
  expect_equal(out$radius_below, c(3.254062, 4.764208, 6.975182), tolerance = 1e-6)
  expect_true(all(out$status == "complete"))
})


test_that("gotham_band distances equal profile-space Euclidean distances", {
  state <- .make_band_state()

  band <- gotham_band(
    state,
    threshold = 0.80,
    epsilon = 0.02,
    alignment = c(-1, 0, 1)
  )

  profiles <- gotham_profile_set(
    state,
    threshold = 0.80,
    epsilon = 0.02,
    alignment = c(-1, 0, 1)
  )$profiles

  for (current_alignment in c(-1, 0, 1)) {
    current_band <- band[band$alignment == current_alignment, , drop = FALSE]

    above <- profiles[
      profiles$alignment == current_alignment & profiles$relation == "above",
      "deviation"
    ]
    threshold_profile <- profiles[
      profiles$alignment == current_alignment & profiles$relation == "equal",
      "deviation"
    ]
    below <- profiles[
      profiles$alignment == current_alignment & profiles$relation == "below",
      "deviation"
    ]

    expect_equal(
      sqrt(sum((threshold_profile - above)^2)),
      current_band$distance_above_to_threshold,
      tolerance = 1e-10
    )
    expect_equal(
      sqrt(sum((below - threshold_profile)^2)),
      current_band$distance_threshold_to_below,
      tolerance = 1e-10
    )
    expect_equal(
      sqrt(sum((below - above)^2)),
      current_band$band_width,
      tolerance = 1e-10
    )
  }
})


test_that("the example band is wider in reinforcing directions", {
  state <- .make_band_state()

  out <- gotham_band(
    state,
    threshold = 0.80,
    epsilon = 0.02,
    alignment = c(-1, 0, 1)
  )

  expect_lt(out$band_width[1], out$band_width[2])
  expect_lt(out$band_width[2], out$band_width[3])
})


test_that("gotham_band returns requested ICC targets", {
  state <- .make_band_state()

  out <- gotham_band(
    state,
    threshold = 0.80,
    epsilon = 0.02,
    alignment = c(-1, 0, 1)
  )

  expect_equal(out$icc_above, rep(0.82, 3))
  expect_equal(out$icc_threshold, rep(0.80, 3))
  expect_equal(out$icc_below, rep(0.78, 3))
})


test_that("gotham_band validates inputs", {
  state <- .make_band_state()

  expect_error(
    gotham_band(state, threshold = 0.80, epsilon = 0),
    "positive"
  )
  expect_error(
    gotham_band(state, threshold = 0.01, epsilon = 0.02),
    "must both lie"
  )
  expect_error(
    gotham_band(state, threshold = 0.99, epsilon = 0.02),
    "must both lie"
  )
  expect_error(
    gotham_band(state, threshold = 0.80, epsilon = 0.02, alignment = 2),
    "interval \\[-1, 1\\]"
  )
  expect_error(
    gotham_band(state, threshold = 0.80, epsilon = 0.02, root = 3),
    "either 1 or 2"
  )
})
