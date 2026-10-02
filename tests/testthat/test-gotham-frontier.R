.make_frontier_state <- function() {
  dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = 1:9
  )

  gotham_state(dat, outcome = "y", group = "group")
}


.make_centered_frontier_deviation <- function(state, radius, alignment) {
  u <- state$group_offsets / sqrt(sum(state$group_offsets^2))

  # For the three-group validation state, this vector is centered and
  # orthogonal to the group-offset direction (-3, 0, 3).
  v <- c(A = 1, B = -2, C = 1)
  v <- v / sqrt(sum(v^2))

  stopifnot(abs(sum(v)) < 1e-12)
  stopifnot(abs(sum(u * v)) < 1e-12)

  radius * (
    alignment * u +
      sqrt(max(0, 1 - alignment^2)) * v
  )
}


test_that("gotham_frontier finds the known first crossings", {
  state <- .make_frontier_state()

  out <- gotham_frontier(
    state,
    threshold = 0.80,
    alignment = c(-1, 0, 1)
  )

  expect_equal(out$frontier$radius_1[1], 3.0260881250, tolerance = 1e-9)
  expect_equal(out$frontier$radius_1[2], 4.3699856056, tolerance = 1e-9)
  expect_equal(out$frontier$radius_1[3], 6.3107131731, tolerance = 1e-9)
  expect_true(all(out$frontier$n_roots == 1))
  expect_true(all(out$frontier$status == "one_root"))
})


test_that("frontier radii reproduce the requested ICC through gotham_update", {
  state <- .make_frontier_state()
  threshold <- 0.80
  alignments <- c(-1, -0.5, 0, 0.5, 1)

  frontier <- gotham_frontier(
    state,
    threshold = threshold,
    alignment = alignments
  )

  for (i in seq_len(nrow(frontier$frontier))) {
    radius <- frontier$frontier$radius_1[i]
    current_alignment <- frontier$frontier$alignment[i]

    deviation <- .make_centered_frontier_deviation(
      state,
      radius = radius,
      alignment = current_alignment
    )

    updated <- gotham_update(state, deviation)

    expect_equal(
      updated$proposal$mean_deviation,
      0,
      tolerance = 1e-12
    )
    expect_equal(
      updated$proposal$alignment,
      current_alignment,
      tolerance = 1e-10
    )
    expect_equal(updated$icc, threshold, tolerance = 1e-10)
    expect_equal(
      updated$ss_between / updated$ss_within,
      frontier$boundary_ss_ratio,
      tolerance = 1e-10
    )
  }
})


test_that("gotham_frontier can identify two crossings and no crossing", {
  state <- .make_frontier_state()

  out <- gotham_frontier(
    state,
    threshold = 0.05,
    alignment = c(-1, 0, 1)
  )

  expect_equal(out$frontier$n_roots, c(2, 0, 0))
  expect_identical(
    out$frontier$status,
    c("two_roots", "no_root", "no_root")
  )
  expect_true(is.finite(out$frontier$radius_1[1]))
  expect_true(is.finite(out$frontier$radius_2[1]))
  expect_true(out$frontier$radius_2[1] > out$frontier$radius_1[1])
  expect_true(is.na(out$frontier$radius_1[2]))
  expect_true(is.na(out$frontier$radius_1[3]))
})


test_that("both roots of a two-crossing ray reproduce the threshold ICC", {
  state <- .make_frontier_state()
  threshold <- 0.05

  out <- gotham_frontier(
    state,
    threshold = threshold,
    alignment = -1
  )

  radii <- c(
    out$frontier$radius_1,
    out$frontier$radius_2
  )

  expect_equal(length(radii), 2L)

  for (radius in radii) {
    deviation <- .make_centered_frontier_deviation(
      state,
      radius = radius,
      alignment = -1
    )

    updated <- gotham_update(state, deviation)
    expect_equal(updated$icc, threshold, tolerance = 1e-9)
  }
})


test_that("frontier summaries have the expected interpretation", {
  state <- .make_frontier_state()

  out <- gotham_frontier(
    state,
    threshold = 0.80,
    alignment = 0
  )

  zero_update <- gotham_update(
    state,
    deviations = c(A = 0, B = 0, C = 0)
  )

  expect_equal(out$icc_at_zero, zero_update$icc, tolerance = 1e-12)
  expect_equal(out$icc_at_infinity, 1 / 9, tolerance = 1e-12)
  expect_true(out$icc_at_zero > out$threshold)
})


test_that("gotham_frontier handles the asymptotic-threshold linear case", {
  state <- .make_frontier_state()
  asymptotic_threshold <- 1 / 9

  out <- gotham_frontier(
    state,
    threshold = asymptotic_threshold,
    alignment = c(-1, 0, 1)
  )

  expect_true(all(abs(out$frontier$quadratic) < 1e-12))
  expect_equal(out$frontier$n_roots[c(1, 3)], c(1, 0))
  expect_true(is.na(out$frontier$discriminant[2]))
})


test_that("gotham_frontier validates its inputs", {
  state <- .make_frontier_state()

  expect_error(
    gotham_frontier(state, threshold = -0.01),
    "interval \\[0, 1\\)"
  )
  expect_error(
    gotham_frontier(state, threshold = 1),
    "interval \\[0, 1\\)"
  )
  expect_error(
    gotham_frontier(state, threshold = 0.5, alignment = 1.1),
    "interval \\[-1, 1\\]"
  )
  expect_error(
    gotham_frontier(state, threshold = 0.5, alignment = NA_real_),
    "finite numeric"
  )
  expect_error(
    gotham_frontier(state, threshold = 0.5, tolerance = 0),
    "positive finite"
  )

  flat_dat <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = rep(c(-1, 0, 1), times = 3)
  )
  flat_state <- gotham_state(
    flat_dat,
    outcome = "y",
    group = "group"
  )

  expect_error(
    gotham_frontier(flat_state, threshold = 0.5),
    "nonzero between-group offsets"
  )
})
