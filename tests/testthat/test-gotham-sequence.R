.make_sequence_data <- function() {
  data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    y = 1:9
  )
}


.make_sequence_state <- function() {
  dat <- .make_sequence_data()
  gotham_state(dat, outcome = "y", group = "group")
}


test_that("gotham_sequence hugs the requested ICC boundary", {
  state <- .make_sequence_state()

  for (current_alignment in c(-1, 0, 1)) {
    out <- gotham_sequence(
      state,
      steps = 5,
      threshold = 0.80,
      alignment = current_alignment
    )

    expect_equal(nrow(out$trajectory), 5L)
    expect_equal(out$trajectory$step, 1:5)
    expect_equal(out$trajectory$group_size_before, 3:7)
    expect_equal(out$trajectory$group_size_after, 4:8)
    expect_true(all(is.finite(out$trajectory$radius)))
    expect_true(all(out$trajectory$radius >= 0))
    expect_equal(
      out$trajectory$alignment_achieved,
      rep(current_alignment, 5),
      tolerance = 1e-10
    )
    expect_equal(
      out$trajectory$icc_after,
      rep(0.80, 5),
      tolerance = 1e-10
    )
    expect_true(all(out$trajectory$icc_at_zero >= 0.80 - 1e-10))
  }
})


test_that("gotham_sequence final state matches brute-force sequential augmentation", {
  dat <- .make_sequence_data()
  state <- gotham_state(dat, outcome = "y", group = "group")

  out <- gotham_sequence(
    state,
    steps = 4,
    threshold = 0.80,
    alignment = c(-1, 0, 0.5, 1)
  )

  brute_data <- dat

  for (step in seq_len(out$steps)) {
    rows <- out$profiles[out$profiles$step == step, , drop = FALSE]
    added <- data.frame(
      group = rows$group,
      y = rows$proposed_observation,
      stringsAsFactors = FALSE
    )
    brute_data <- rbind(brute_data, added)
  }

  brute <- gotham_state(
    brute_data,
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
    expect_equal(
      out$final_state[[field]],
      brute[[field]],
      tolerance = 1e-9
    )
  }
})


test_that("each sequential radius is the first frontier crossing", {
  state <- .make_sequence_state()

  out <- gotham_sequence(
    state,
    steps = 3,
    threshold = 0.80,
    alignment = 0.5
  )

  for (step in seq_len(out$steps)) {
    current_state <- out$states[[step]]
    frontier <- gotham_frontier(
      current_state,
      threshold = 0.80,
      alignment = 0.5
    )

    expect_equal(
      out$trajectory$radius[step],
      frontier$frontier$radius_1,
      tolerance = 1e-10
    )

    midpoint <- gotham_profile(
      current_state,
      radius = out$trajectory$radius[step] / 2,
      alignment = 0.5
    )

    expect_gte(midpoint$achieved$icc, 0.80 - 1e-10)
  }
})


test_that("step-specific alignment and orientation schedules are honored", {
  state <- .make_sequence_state()

  out <- gotham_sequence(
    state,
    steps = 4,
    threshold = 0.80,
    alignment = c(-1, -0.25, 0.25, 1),
    orientation = c(1, -1, 1, -1)
  )

  expect_equal(
    out$trajectory$alignment_requested,
    c(-1, -0.25, 0.25, 1)
  )
  expect_equal(
    out$trajectory$alignment_achieved,
    c(-1, -0.25, 0.25, 1),
    tolerance = 1e-10
  )
  expect_equal(out$trajectory$orientation, c(1, -1, 1, -1))
})


test_that("sequence profile rows agree with trajectory and updated states", {
  state <- .make_sequence_state()

  out <- gotham_sequence(
    state,
    steps = 3,
    threshold = 0.80,
    alignment = 0
  )

  expect_equal(nrow(out$profiles), 3L * state$n_groups)

  for (step in seq_len(out$steps)) {
    rows <- out$profiles[out$profiles$step == step, , drop = FALSE]
    current_state <- out$states[[step]]
    updated_state <- out$states[[step + 1L]]

    expect_equal(
      rows$current_group_mean,
      unname(current_state$group_means),
      tolerance = 1e-12
    )
    expect_equal(
      rows$updated_group_mean,
      unname(updated_state$group_means),
      tolerance = 1e-12
    )
    expect_equal(
      rows$achieved_icc,
      rep(updated_state$icc, state$n_groups),
      tolerance = 1e-12
    )
  }
})


test_that("gotham_sequence validates schedule inputs", {
  state <- .make_sequence_state()

  expect_error(
    gotham_sequence(state, steps = 0, threshold = 0.80, alignment = 0),
    "positive integer"
  )
  expect_error(
    gotham_sequence(state, steps = 2, threshold = 1, alignment = 0),
    "interval \\[0, 1\\)"
  )
  expect_error(
    gotham_sequence(
      state,
      steps = 3,
      threshold = 0.80,
      alignment = c(-1, 1)
    ),
    "length 1 or length equal to steps"
  )
  expect_error(
    gotham_sequence(state, steps = 2, threshold = 0.80, alignment = 1.1),
    "interval \\[-1, 1\\]"
  )
  expect_error(
    gotham_sequence(
      state,
      steps = 2,
      threshold = 0.80,
      alignment = 0,
      orientation = c(1, 0)
    ),
    "only -1 or 1"
  )
})
