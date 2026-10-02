#' Construct a representative Gotham proposal profile
#'
#' Constructs one deterministic vector of proposed deviations having a requested
#' centered magnitude and alignment with the current between-group structure.
#' The resulting deviations can be supplied directly to [gotham_update()].
#'
#' Let \eqn{\mathbf{a}} be the current vector of group offsets and let
#' \eqn{\mathbf{d}_c} be the centered proposed-deviation vector. This function
#' constructs a profile satisfying
#' \deqn{||\mathbf{d}_c|| = r}
#' and, when \eqn{r > 0},
#' \deqn{\frac{\mathbf{a}^T\mathbf{d}_c}
#' {||\mathbf{a}||\,||\mathbf{d}_c||} = c,}
#' where `radius` is \eqn{r} and `alignment` is \eqn{c}.
#'
#' For three or more groups, infinitely many profiles can share the same
#' magnitude and alignment. `gotham_profile()` returns one reproducible
#' representative by constructing a canonical centered direction orthogonal to
#' the current group-offset vector. Changing `orientation` reflects that
#' orthogonal component without changing the requested Gotham macrostate.
#'
#' @param state A state returned by [gotham_state()] or [gotham_update()].
#' @param radius Non-negative finite numeric scalar giving the Euclidean norm
#'   of the centered proposed-deviation vector.
#' @param alignment Finite numeric scalar in the interval \eqn{[-1, 1]} giving
#'   cosine alignment with the current group-offset vector.
#' @param mean_deviation Finite numeric scalar added equally to every proposed
#'   deviation. Defaults to zero. This does not change centered magnitude or
#'   alignment, but it does change within-group variation and therefore can
#'   change the updated ICC.
#' @param orientation Either `1` or `-1`. For alignments strictly between -1
#'   and 1, this chooses between two reflected canonical representatives of the
#'   same centered macrostate.
#' @param tolerance Positive finite numeric scalar for numerical decisions.
#'
#' @return A list containing the proposed deviations, proposed observations,
#'   requested geometry, achieved geometry, and the state returned by
#'   [gotham_update()].
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 3),
#'   y = 1:9
#' )
#'
#' state <- gotham_state(dat, outcome = "y", group = "group")
#' gotham_profile(state, radius = 2, alignment = 0)$deviations

gotham_profile <- function(
    state,
    radius,
    alignment,
    mean_deviation = 0,
    orientation = 1,
    tolerance = sqrt(.Machine$double.eps)) {
  .validate_gotham_state(state)

  if (!is.numeric(radius) || length(radius) != 1L ||
      is.na(radius) || !is.finite(radius) || radius < 0) {
    stop("radius must be a non-negative finite numeric scalar.", call. = FALSE)
  }

  if (!is.numeric(alignment) || length(alignment) != 1L ||
      is.na(alignment) || !is.finite(alignment) ||
      alignment < -1 || alignment > 1) {
    stop(
      "alignment must be a finite numeric scalar in the interval [-1, 1].",
      call. = FALSE
    )
  }

  if (!is.numeric(mean_deviation) || length(mean_deviation) != 1L ||
      is.na(mean_deviation) || !is.finite(mean_deviation)) {
    stop("mean_deviation must be a finite numeric scalar.", call. = FALSE)
  }

  if (!is.numeric(orientation) || length(orientation) != 1L ||
      is.na(orientation) || !is.finite(orientation) ||
      !orientation %in% c(-1, 1)) {
    stop("orientation must be either -1 or 1.", call. = FALSE)
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be a positive finite numeric scalar.", call. = FALSE)
  }

  offsets <- state$group_offsets - mean(state$group_offsets)
  offset_norm <- sqrt(sum(offsets^2))

  if (offset_norm <= tolerance) {
    stop(
      paste0(
        "gotham_profile() requires nonzero between-group offsets so that ",
        "proposal alignment is defined."
      ),
      call. = FALSE
    )
  }

  group_labels <- names(state$group_means)
  u <- offsets / offset_norm

  if (radius <= tolerance) {
    centered <- stats::setNames(rep(0, state$n_groups), group_labels)
  } else if (abs(abs(alignment) - 1) <= tolerance) {
    direction <- if (alignment >= 0) u else -u
    centered <- radius * direction
  } else {
    if (state$n_groups < 3L) {
      stop(
        paste0(
          "With two groups, a nonzero centered proposal can only have ",
          "alignment -1 or 1."
        ),
        call. = FALSE
      )
    }

    v <- .gotham_canonical_orthogonal_direction(
      u = u,
      orientation = orientation,
      tolerance = tolerance
    )

    centered <- radius * (
      alignment * u +
        sqrt(max(0, 1 - alignment^2)) * v
    )
  }

  centered <- as.numeric(centered)
  centered <- centered - mean(centered)
  names(centered) <- group_labels

  # Rescale after the final centering step to eliminate accumulated numerical
  # drift while preserving the requested direction.
  centered_norm <- sqrt(sum(centered^2))
  if (radius > tolerance && centered_norm > 0) {
    centered <- centered * radius / centered_norm
  }

  deviations <- centered + mean_deviation
  names(deviations) <- group_labels

  updated <- gotham_update(state, deviations)
  proposed_observations <- state$group_means + deviations

  list(
    deviations = deviations,
    proposed_observations = proposed_observations,
    requested = list(
      radius = unname(radius),
      alignment = unname(alignment),
      mean_deviation = unname(mean_deviation),
      orientation = unname(orientation)
    ),
    achieved = list(
      radius = sqrt(updated$proposal$centered_deviation_ss),
      alignment = updated$proposal$alignment,
      mean_deviation = updated$proposal$mean_deviation,
      icc = updated$icc
    ),
    updated_state = updated
  )
}


#' Generate representative profiles around an ICC Gotham frontier
#'
#' Generates concrete one-observation-per-group profiles at an ICC threshold
#' and, optionally, at equally spaced ICC targets immediately above and below
#' that threshold. This turns the low-dimensional Gotham frontier into actual
#' candidate observations.
#'
#' With `epsilon > 0`, the three target ICC values are
#' \deqn{\tau + \epsilon, \quad \tau, \quad \tau - \epsilon.}
#' These are reported literally as `threshold_plus_epsilon`, `threshold`, and
#' `threshold_minus_epsilon`; no stronger interpretation is imposed because
#' Gotham rays can be non-monotone in more general settings.
#'
#' This function uses the centered frontier from [gotham_frontier()], so the
#' generated profiles have mean proposed deviation equal to zero.
#'
#' @param state A state returned by [gotham_state()] or [gotham_update()].
#' @param threshold Numeric scalar in the interval \eqn{[0, 1)} giving the ICC
#'   boundary of interest.
#' @param epsilon Non-negative finite numeric scalar. The values
#'   `threshold - epsilon` and `threshold + epsilon` must both remain in the
#'   interval \eqn{[0, 1)}.
#' @param alignment Numeric vector of cosine alignments in the interval
#'   \eqn{[-1, 1]}. Defaults to `c(-1, 0, 1)`.
#' @param root Either `1` or `2`, selecting the first or second finite
#'   non-negative frontier crossing when available.
#' @param orientation Either `1` or `-1`, passed to `gotham_profile()`.
#' @param tolerance Positive finite numeric scalar for numerical decisions.
#'
#' @return A list containing `summary`, with one row per target ICC and
#'   alignment, and `profiles`, with one row per generated group-level proposed
#'   observation. Rows for unavailable frontier roots are retained in `summary`
#'   with status `no_profile` and are omitted from `profiles`.
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 3),
#'   y = 1:9
#' )
#'
#' state <- gotham_state(dat, outcome = "y", group = "group")
#' gotham_profile_set(
#'   state,
#'   threshold = 0.80,
#'   epsilon = 0.02,
#'   alignment = c(-1, 0, 1)
#' )$summary

gotham_profile_set <- function(
    state,
    threshold,
    epsilon = 0,
    alignment = c(-1, 0, 1),
    root = 1L,
    orientation = 1,
    tolerance = sqrt(.Machine$double.eps)) {
  .validate_gotham_state(state)

  if (!is.numeric(threshold) || length(threshold) != 1L ||
      is.na(threshold) || !is.finite(threshold) ||
      threshold < 0 || threshold >= 1) {
    stop(
      "threshold must be a finite numeric scalar in the interval [0, 1).",
      call. = FALSE
    )
  }

  if (!is.numeric(epsilon) || length(epsilon) != 1L ||
      is.na(epsilon) || !is.finite(epsilon) || epsilon < 0) {
    stop("epsilon must be a non-negative finite numeric scalar.", call. = FALSE)
  }

  if (threshold - epsilon < 0 || threshold + epsilon >= 1) {
    stop(
      paste0(
        "threshold - epsilon and threshold + epsilon must both lie in ",
        "the interval [0, 1)."
      ),
      call. = FALSE
    )
  }

  if (!is.numeric(alignment) || length(alignment) == 0L ||
      anyNA(alignment) || any(!is.finite(alignment)) ||
      any(alignment < -1) || any(alignment > 1)) {
    stop(
      "alignment must contain finite numeric values in the interval [-1, 1].",
      call. = FALSE
    )
  }

  if (!is.numeric(root) || length(root) != 1L || is.na(root) ||
      !is.finite(root) || !root %in% c(1, 2)) {
    stop("root must be either 1 or 2.", call. = FALSE)
  }
  root <- as.integer(root)

  if (!is.numeric(orientation) || length(orientation) != 1L ||
      is.na(orientation) || !is.finite(orientation) ||
      !orientation %in% c(-1, 1)) {
    stop("orientation must be either -1 or 1.", call. = FALSE)
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be a positive finite numeric scalar.", call. = FALSE)
  }

  if (epsilon == 0) {
    targets <- data.frame(
      target = "threshold",
      relation = "equal",
      target_icc = threshold,
      stringsAsFactors = FALSE
    )
  } else {
    targets <- data.frame(
      target = c(
        "threshold_plus_epsilon",
        "threshold",
        "threshold_minus_epsilon"
      ),
      relation = c("above", "equal", "below"),
      target_icc = c(
        threshold + epsilon,
        threshold,
        threshold - epsilon
      ),
      stringsAsFactors = FALSE
    )
  }

  summary_rows <- list()
  profile_rows <- list()
  summary_index <- 0L
  profile_index <- 0L

  for (target_index in seq_len(nrow(targets))) {
    current_target <- targets[target_index, , drop = FALSE]

    frontier <- gotham_frontier(
      state = state,
      threshold = current_target$target_icc,
      alignment = alignment,
      tolerance = tolerance
    )

    radius_name <- paste0("radius_", root)

    for (alignment_index in seq_len(nrow(frontier$frontier))) {
      frontier_row <- frontier$frontier[alignment_index, , drop = FALSE]
      radius <- frontier_row[[radius_name]]
      available <- is.finite(radius)

      summary_index <- summary_index + 1L

      if (!available) {
        summary_rows[[summary_index]] <- data.frame(
          target = current_target$target,
          relation = current_target$relation,
          target_icc = current_target$target_icc,
          alignment = frontier_row$alignment,
          root = root,
          radius = NA_real_,
          achieved_alignment = NA_real_,
          achieved_icc = NA_real_,
          status = "no_profile",
          stringsAsFactors = FALSE
        )
        next
      }

      profile <- gotham_profile(
        state = state,
        radius = radius,
        alignment = frontier_row$alignment,
        mean_deviation = 0,
        orientation = orientation,
        tolerance = tolerance
      )

      summary_rows[[summary_index]] <- data.frame(
        target = current_target$target,
        relation = current_target$relation,
        target_icc = current_target$target_icc,
        alignment = frontier_row$alignment,
        root = root,
        radius = radius,
        achieved_alignment = profile$achieved$alignment,
        achieved_icc = profile$achieved$icc,
        status = "generated",
        stringsAsFactors = FALSE
      )

      profile_index <- profile_index + 1L
      profile_rows[[profile_index]] <- data.frame(
        target = current_target$target,
        relation = current_target$relation,
        target_icc = current_target$target_icc,
        alignment = frontier_row$alignment,
        root = root,
        radius = radius,
        group = names(profile$deviations),
        current_group_mean = unname(state$group_means),
        deviation = unname(profile$deviations),
        proposed_observation = unname(profile$proposed_observations),
        achieved_icc = profile$achieved$icc,
        stringsAsFactors = FALSE
      )
    }
  }

  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL

  profiles <- if (length(profile_rows) == 0L) {
    data.frame(
      target = character(),
      relation = character(),
      target_icc = numeric(),
      alignment = numeric(),
      root = integer(),
      radius = numeric(),
      group = character(),
      current_group_mean = numeric(),
      deviation = numeric(),
      proposed_observation = numeric(),
      achieved_icc = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, profile_rows)
    rownames(out) <- NULL
    out
  }

  list(
    threshold = unname(threshold),
    epsilon = unname(epsilon),
    root = root,
    orientation = unname(orientation),
    summary = summary,
    profiles = profiles
  )
}


.gotham_canonical_orthogonal_direction <- function(
    u,
    orientation,
    tolerance) {
  g <- length(u)
  candidates <- vector("list", g)
  norms <- numeric(g)

  for (j in seq_len(g)) {
    candidate <- rep(0, g)
    candidate[j] <- 1
    candidate <- candidate - mean(candidate)
    candidate <- candidate - sum(candidate * u) * u
    candidate <- candidate - mean(candidate)

    candidate_norm <- sqrt(sum(candidate^2))
    candidates[[j]] <- candidate
    norms[j] <- candidate_norm
  }

  best <- which.max(norms)

  if (length(best) == 0L || norms[best] <= tolerance) {
    stop(
      "Could not construct a centered direction orthogonal to group offsets.",
      call. = FALSE
    )
  }

  orientation * candidates[[best]] / norms[best]
}
