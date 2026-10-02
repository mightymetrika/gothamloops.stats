#' Quantify an ICC Gotham threshold band in observation space
#'
#' Converts an ICC-scale epsilon band around a Gotham frontier into radial
#' distances in the centered proposal space. This makes it possible to compare
#' how much the proposed observations must move, along a fixed alignment ray,
#' to go from just above a clustering threshold to the threshold and then just
#' below it.
#'
#' For a threshold \eqn{\tau} and positive `epsilon`, the function evaluates
#' the centered Gotham frontier at
#' \deqn{\tau + \epsilon, \quad \tau, \quad \tau - \epsilon.}
#' For each requested alignment it reports the selected frontier radius at all
#' three ICC values. The Euclidean distance between two representative
#' centered profiles on the same ray is the absolute difference between their
#' radii, so the returned radial gaps are also observation-space distances for
#' such paired profiles.
#'
#' The function does not assume that every frontier ray is monotone or that a
#' selected root exists at every target ICC. If any selected radius is
#' unavailable, the corresponding distances are returned as `NA` and the row
#' status is `incomplete`.
#'
#' @param state A state returned by [gotham_state()] or [gotham_update()].
#' @param threshold Numeric scalar in the interval \eqn{(0, 1)} giving the ICC
#'   boundary of interest.
#' @param epsilon Positive finite numeric scalar. `threshold - epsilon` and
#'   `threshold + epsilon` must both lie in the interval \eqn{[0, 1)}.
#' @param alignment Numeric vector of cosine alignments in the interval
#'   \eqn{[-1, 1]}. Defaults to 201 equally spaced values from -1 to 1.
#' @param root Either `1` or `2`, selecting the first or second finite
#'   non-negative frontier crossing when available.
#' @param tolerance Positive finite numeric scalar for numerical decisions.
#'
#' @return A data frame with one row per requested alignment and columns for
#'   the three ICC targets, their frontier radii, the Euclidean radial distance
#'   from `threshold + epsilon` to the threshold, the distance from the
#'   threshold to `threshold - epsilon`, the full band width, and a status.
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 3),
#'   y = 1:9
#' )
#'
#' state <- gotham_state(dat, outcome = "y", group = "group")
#' gotham_band(
#'   state,
#'   threshold = 0.80,
#'   epsilon = 0.02,
#'   alignment = c(-1, 0, 1)
#' )
gotham_band <- function(
    state,
    threshold,
    epsilon,
    alignment = seq(-1, 1, length.out = 201L),
    root = 1L,
    tolerance = sqrt(.Machine$double.eps)) {
  .validate_gotham_state(state)

  if (!is.numeric(threshold) || length(threshold) != 1L ||
      is.na(threshold) || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1) {
    stop(
      "threshold must be a finite numeric scalar in the interval (0, 1).",
      call. = FALSE
    )
  }

  if (!is.numeric(epsilon) || length(epsilon) != 1L ||
      is.na(epsilon) || !is.finite(epsilon) || epsilon <= 0) {
    stop("epsilon must be a positive finite numeric scalar.", call. = FALSE)
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

  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be a positive finite numeric scalar.", call. = FALSE)
  }

  targets <- c(
    above = threshold + epsilon,
    threshold = threshold,
    below = threshold - epsilon
  )

  frontiers <- lapply(
    targets,
    function(current_threshold) {
      gotham_frontier(
        state = state,
        threshold = current_threshold,
        alignment = alignment,
        tolerance = tolerance
      )$frontier
    }
  )

  radius_name <- paste0("radius_", root)

  radius_above <- frontiers$above[[radius_name]]
  radius_threshold <- frontiers$threshold[[radius_name]]
  radius_below <- frontiers$below[[radius_name]]

  complete <- is.finite(radius_above) &
    is.finite(radius_threshold) &
    is.finite(radius_below)

  distance_above_to_threshold <- ifelse(
    complete,
    abs(radius_threshold - radius_above),
    NA_real_
  )

  distance_threshold_to_below <- ifelse(
    complete,
    abs(radius_below - radius_threshold),
    NA_real_
  )

  band_width <- ifelse(
    complete,
    abs(radius_below - radius_above),
    NA_real_
  )

  data.frame(
    alignment = alignment,
    root = root,
    icc_above = unname(targets[["above"]]),
    radius_above = radius_above,
    icc_threshold = unname(targets[["threshold"]]),
    radius_threshold = radius_threshold,
    icc_below = unname(targets[["below"]]),
    radius_below = radius_below,
    distance_above_to_threshold = distance_above_to_threshold,
    distance_threshold_to_below = distance_threshold_to_below,
    band_width = band_width,
    status = ifelse(complete, "complete", "incomplete"),
    stringsAsFactors = FALSE
  )
}
