#' Update a Gotham ANOVA state with one new observation per group
#'
#' Analytically updates a balanced one-way ANOVA Gotham state after adding one
#' new observation to every group. No raw data are required once the starting
#' state has been constructed with [gotham_state()].
#'
#' For group \eqn{g}, let the proposed new observation be
#' \deqn{x_g = \bar{y}_g + d_g,}
#' where \eqn{\bar{y}_g} is the current group mean and \eqn{d_g} is supplied in
#' `deviations`. If each group currently contains \eqn{n} observations, then
#' the within-group sum of squares updates as
#' \deqn{SS_W' = SS_W + \frac{n}{n+1}\sum_g d_g^2.}
#'
#' Let \eqn{a_g = \bar{y}_g - \bar{y}}, \eqn{\bar{d}} be the mean of the
#' proposed deviations, and \eqn{d_{cg} = d_g - \bar{d}}. The updated
#' between-group sum of squares is
#' \deqn{SS_B' = (n+1)\sum_g\left(a_g + \frac{d_{cg}}{n+1}\right)^2.}
#'
#' The returned object has the same core state fields as [gotham_state()], so
#' it can be supplied to a later call to `gotham_update()` for sequential
#' updates. The `proposal` component records the geometry of the just-added
#' deviations.
#'
#' @param state A state returned by [gotham_state()] or a previous call to
#'   `gotham_update()`.
#' @param deviations Numeric vector with one proposed deviation per group. A
#'   deviation is the proposed new observation minus that group's current mean.
#'   If named, the names must match the state's group labels; named deviations
#'   are reordered to the state order before calculation.
#'
#' @return A named list containing the updated Gotham state plus a `proposal`
#'   component with:
#' \describe{
#'   \item{deviations}{Named vector of proposed deviations in state order.}
#'   \item{mean_deviation}{Mean proposed deviation across groups.}
#'   \item{centered_deviations}{Deviations after subtracting their mean.}
#'   \item{deviation_ss}{Sum of squared raw deviations.}
#'   \item{centered_deviation_ss}{Sum of squared centered deviations.}
#'   \item{alignment}{Cosine alignment between the current group-offset vector
#'     and the centered deviation vector. `NA` when either vector has zero
#'     length.}
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
#' gotham_update(state, deviations = c(A = -1, B = 0, C = 1))
gotham_update <- function(state, deviations) {
  .validate_gotham_state(state)

  group_labels <- names(state$group_means)

  if (!is.numeric(deviations) || length(deviations) != state$n_groups) {
    stop(
      "deviations must be a numeric vector with one value per group.",
      call. = FALSE
    )
  }

  if (anyNA(deviations) || any(!is.finite(deviations))) {
    stop(
      "deviations must contain only finite, non-missing values.",
      call. = FALSE
    )
  }

  if (!is.null(names(deviations))) {
    if (anyNA(names(deviations)) || any(!nzchar(names(deviations))) ||
        anyDuplicated(names(deviations))) {
      stop(
        "Named deviations must have unique, non-empty group names.",
        call. = FALSE
      )
    }

    if (!setequal(names(deviations), group_labels)) {
      stop(
        "Named deviations must match the state's group labels exactly.",
        call. = FALSE
      )
    }

    deviations <- deviations[group_labels]
  } else {
    names(deviations) <- group_labels
  }

  n <- state$group_size
  n_new <- n + 1L

  mean_deviation <- mean(deviations)
  centered_deviations <- deviations - mean_deviation

  new_grand_mean <- state$grand_mean + mean_deviation / n_new
  new_group_means <- state$group_means + deviations / n_new
  new_group_offsets <- state$group_offsets + centered_deviations / n_new

  new_ss_within <- state$ss_within +
    (n / n_new) * sum(deviations^2)
  new_group_offset_ss <- sum(new_group_offsets^2)
  new_ss_between <- n_new * new_group_offset_ss

  df_between <- state$n_groups - 1L
  df_within <- state$n_groups * (n_new - 1L)
  new_ms_between <- new_ss_between / df_between
  new_ms_within <- new_ss_within / df_within

  icc_components <- list(
    group_size = n_new,
    ms_between = new_ms_between,
    ms_within = new_ms_within
  )

  offset_norm <- sqrt(state$group_offset_ss)
  centered_deviation_norm <- sqrt(sum(centered_deviations^2))

  alignment <- if (offset_norm > 0 && centered_deviation_norm > 0) {
    unname(
      sum(state$group_offsets * centered_deviations) /
        (offset_norm * centered_deviation_norm)
    )
  } else {
    NA_real_
  }

  list(
    n_groups = state$n_groups,
    group_size = n_new,
    grand_mean = unname(new_grand_mean),
    group_means = new_group_means,
    group_offsets = new_group_offsets,
    group_offset_ss = unname(new_group_offset_ss),
    ss_between = unname(new_ss_between),
    ss_within = unname(new_ss_within),
    ms_between = unname(new_ms_between),
    ms_within = unname(new_ms_within),
    icc = .anova_icc_from_components(icc_components),
    proposal = list(
      deviations = deviations,
      mean_deviation = unname(mean_deviation),
      centered_deviations = centered_deviations,
      deviation_ss = unname(sum(deviations^2)),
      centered_deviation_ss = unname(sum(centered_deviations^2)),
      alignment = alignment
    )
  )
}


.validate_gotham_state <- function(state) {
  required <- c(
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

  if (!is.list(state) || !all(required %in% names(state))) {
    stop(
      "state must be a Gotham state returned by gotham_state() or gotham_update().",
      call. = FALSE
    )
  }

  if (!is.numeric(state$n_groups) || length(state$n_groups) != 1L ||
      !is.finite(state$n_groups) || state$n_groups < 2 ||
      state$n_groups != as.integer(state$n_groups)) {
    stop("state has an invalid n_groups value.", call. = FALSE)
  }

  if (!is.numeric(state$group_size) || length(state$group_size) != 1L ||
      !is.finite(state$group_size) || state$group_size < 2 ||
      state$group_size != as.integer(state$group_size)) {
    stop("state has an invalid group_size value.", call. = FALSE)
  }

  if (!is.numeric(state$group_means) ||
      length(state$group_means) != state$n_groups ||
      is.null(names(state$group_means)) ||
      anyNA(state$group_means) || any(!is.finite(state$group_means))) {
    stop("state has invalid group_means.", call. = FALSE)
  }

  if (!is.numeric(state$group_offsets) ||
      length(state$group_offsets) != state$n_groups ||
      !identical(names(state$group_offsets), names(state$group_means)) ||
      anyNA(state$group_offsets) || any(!is.finite(state$group_offsets))) {
    stop("state has invalid group_offsets.", call. = FALSE)
  }

  scalar_fields <- c(
    "grand_mean",
    "group_offset_ss",
    "ss_between",
    "ss_within",
    "ms_between",
    "ms_within",
    "icc"
  )

  scalar_ok <- vapply(
    scalar_fields,
    function(field) {
      value <- state[[field]]
      is.numeric(value) && length(value) == 1L &&
        !is.na(value) && is.finite(value)
    },
    logical(1)
  )

  if (!all(scalar_ok)) {
    stop("state contains invalid numeric summary fields.", call. = FALSE)
  }

  invisible(TRUE)
}
