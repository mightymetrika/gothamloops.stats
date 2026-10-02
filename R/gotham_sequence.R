#' Generate a sequential boundary-hugging Gotham process
#'
#' Repeatedly adds one observation per group while moving as far as possible
#' along a requested centered proposal direction without leaving the connected
#' ICC-admissible region that contains zero proposal magnitude.
#'
#' At each step, [gotham_frontier()] is evaluated for the current state and the
#' requested alignment. The first finite non-negative frontier crossing is used
#' as the step radius. A representative profile at that radius is constructed
#' with [gotham_profile()], added analytically, and the resulting state becomes
#' the starting state for the next step.
#'
#' Thus, when the zero-radius proposal has ICC at least `threshold`, each
#' completed step ends on
#' \deqn{ICC = \tau,}
#' where \eqn{\tau} is `threshold`. This first sequential process fixes the mean
#' proposed deviation at zero, so each step changes the relative group
#' structure without imposing a common shift across all groups.
#'
#' `alignment` can be a scalar, which is reused at every step, or a vector with
#' one value per step. Likewise, `orientation` can be fixed or supplied as a
#' step-specific vector. For alignments strictly between -1 and 1, orientation
#' selects between reflected representative microstates with the same centered
#' Gotham geometry.
#'
#' @param state A state returned by [gotham_state()] or [gotham_update()].
#' @param steps Positive integer number of sequential updates.
#' @param threshold Numeric scalar in the interval \eqn{[0, 1)} giving the ICC
#'   boundary to preserve.
#' @param alignment Numeric scalar or length-`steps` vector of cosine alignments
#'   in the interval \eqn{[-1, 1]}.
#' @param orientation Numeric scalar or length-`steps` vector containing only
#'   `1` or `-1`, passed to [gotham_profile()].
#' @param tolerance Positive finite numeric scalar for numerical decisions.
#'
#' @return A list containing:
#' \describe{
#'   \item{threshold}{Requested ICC threshold.}
#'   \item{steps}{Number of completed sequential updates.}
#'   \item{trajectory}{One row per step summarizing the frontier radius,
#'     requested and achieved alignment, ICC before and after the step, the
#'     zero-radius ICC for that step, and key ANOVA state quantities.}
#'   \item{profiles}{One row per group per step containing the current group
#'     mean, proposed deviation, proposed observation, and updated group mean.}
#'   \item{states}{A list containing the initial state followed by every
#'     updated state.}
#'   \item{initial_state}{The supplied starting state.}
#'   \item{final_state}{The state after the final update.}
#' }
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 3),
#'   y = 1:9
#' )
#'
#' state <- gotham_state(dat, outcome = "y", group = "group")
#' gotham_sequence(
#'   state,
#'   steps = 3,
#'   threshold = 0.80,
#'   alignment = -1
#' )$trajectory

gotham_sequence <- function(
    state,
    steps,
    threshold,
    alignment,
    orientation = 1,
    tolerance = sqrt(.Machine$double.eps)) {
  .validate_gotham_state(state)

  if (!is.numeric(steps) || length(steps) != 1L ||
      is.na(steps) || !is.finite(steps) || steps < 1 ||
      steps != as.integer(steps)) {
    stop("steps must be a positive integer scalar.", call. = FALSE)
  }
  steps <- as.integer(steps)

  if (!is.numeric(threshold) || length(threshold) != 1L ||
      is.na(threshold) || !is.finite(threshold) ||
      threshold < 0 || threshold >= 1) {
    stop(
      "threshold must be a finite numeric scalar in the interval [0, 1).",
      call. = FALSE
    )
  }

  alignment <- .gotham_recycle_sequence_argument(
    alignment,
    steps = steps,
    name = "alignment"
  )

  if (anyNA(alignment) || any(!is.finite(alignment)) ||
      any(alignment < -1) || any(alignment > 1)) {
    stop(
      "alignment must contain finite numeric values in the interval [-1, 1].",
      call. = FALSE
    )
  }

  orientation <- .gotham_recycle_sequence_argument(
    orientation,
    steps = steps,
    name = "orientation"
  )

  if (anyNA(orientation) || any(!is.finite(orientation)) ||
      any(!orientation %in% c(-1, 1))) {
    stop("orientation must contain only -1 or 1.", call. = FALSE)
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be a positive finite numeric scalar.", call. = FALSE)
  }

  current_state <- state
  states <- vector("list", steps + 1L)
  states[[1L]] <- state
  trajectory_rows <- vector("list", steps)
  profile_rows <- vector("list", steps)

  for (step in seq_len(steps)) {
    current_alignment <- alignment[[step]]
    current_orientation <- orientation[[step]]

    frontier <- gotham_frontier(
      state = current_state,
      threshold = threshold,
      alignment = current_alignment,
      tolerance = tolerance
    )

    if (frontier$icc_at_zero < threshold - tolerance) {
      stop(
        paste0(
          "At step ", step,
          ", the zero-radius proposal is already below the requested ICC ",
          "threshold. The connected admissible region containing radius zero ",
          "does not reach the threshold frontier."
        ),
        call. = FALSE
      )
    }

    frontier_row <- frontier$frontier[1L, , drop = FALSE]
    radius <- frontier_row$radius_1[[1L]]

    if (!is.finite(radius)) {
      stop(
        paste0(
          "At step ", step,
          ", no finite first frontier radius is available for alignment ",
          format(current_alignment), "."
        ),
        call. = FALSE
      )
    }

    profile <- gotham_profile(
      state = current_state,
      radius = radius,
      alignment = current_alignment,
      mean_deviation = 0,
      orientation = current_orientation,
      tolerance = tolerance
    )

    updated_state <- profile$updated_state

    icc_tolerance <- tolerance * max(1, abs(threshold))
    if (abs(updated_state$icc - threshold) > 100 * icc_tolerance) {
      stop(
        paste0(
          "At step ", step,
          ", the generated boundary profile did not reproduce the requested ",
          "ICC threshold within numerical tolerance."
        ),
        call. = FALSE
      )
    }

    trajectory_rows[[step]] <- data.frame(
      step = step,
      group_size_before = current_state$group_size,
      group_size_after = updated_state$group_size,
      alignment_requested = current_alignment,
      alignment_achieved = profile$achieved$alignment,
      orientation = current_orientation,
      radius = profile$achieved$radius,
      icc_before = current_state$icc,
      icc_at_zero = frontier$icc_at_zero,
      icc_after = updated_state$icc,
      group_offset_ss_before = current_state$group_offset_ss,
      group_offset_ss_after = updated_state$group_offset_ss,
      ss_between_before = current_state$ss_between,
      ss_between_after = updated_state$ss_between,
      ss_within_before = current_state$ss_within,
      ss_within_after = updated_state$ss_within,
      stringsAsFactors = FALSE
    )

    profile_rows[[step]] <- data.frame(
      step = step,
      alignment_requested = current_alignment,
      alignment_achieved = profile$achieved$alignment,
      radius = profile$achieved$radius,
      group = names(profile$deviations),
      current_group_mean = unname(current_state$group_means),
      deviation = unname(profile$deviations),
      proposed_observation = unname(profile$proposed_observations),
      updated_group_mean = unname(updated_state$group_means),
      achieved_icc = updated_state$icc,
      stringsAsFactors = FALSE
    )

    states[[step + 1L]] <- updated_state
    current_state <- updated_state
  }

  trajectory <- do.call(rbind, trajectory_rows)
  rownames(trajectory) <- NULL

  profiles <- do.call(rbind, profile_rows)
  rownames(profiles) <- NULL

  list(
    threshold = unname(threshold),
    steps = steps,
    trajectory = trajectory,
    profiles = profiles,
    states = states,
    initial_state = state,
    final_state = current_state
  )
}


.gotham_recycle_sequence_argument <- function(x, steps, name) {
  if (!is.numeric(x) || length(x) == 0L) {
    stop(
      paste0(name, " must be a numeric scalar or a numeric vector with one value per step."),
      call. = FALSE
    )
  }

  if (length(x) == 1L) {
    return(rep(x, steps))
  }

  if (length(x) != steps) {
    stop(
      paste0(name, " must have length 1 or length equal to steps."),
      call. = FALSE
    )
  }

  x
}
