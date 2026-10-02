#' Derive the centered ANOVA ICC Gotham frontier
#'
#' Derives the boundary between centered one-observation-per-group proposals
#' whose updated balanced one-way ANOVA ICC is above versus below a specified
#' threshold. The frontier is parameterized by proposal magnitude and by its
#' alignment with the current vector of group offsets.
#'
#' This first frontier fixes the mean proposed deviation at zero. Let
#' \eqn{\mathbf{a}} be the current vector of group means minus the grand mean,
#' let \eqn{r = ||\mathbf{d}||} be the Euclidean norm of a centered proposal,
#' and let \eqn{c} be its cosine alignment with \eqn{\mathbf{a}}. If each group
#' currently contains \eqn{n} observations, then after adding one observation
#' per group,
#' \deqn{SS_W' = SS_W + \frac{n}{n+1}r^2}
#' and
#' \deqn{SS_B' = (n+1)||\mathbf{a}||^2 + 2||\mathbf{a}||rc +
#'       \frac{r^2}{n+1}.}
#'
#' Setting the updated ICC equal to `threshold` gives a quadratic equation in
#' \eqn{r}. Depending on the starting state, threshold, and alignment, a ray can
#' have zero, one, or two finite non-negative crossings. Two crossings are
#' possible because the ANOVA ICC along a ray need not be monotone in \eqn{r}.
#'
#' @param state A state returned by [gotham_state()] or [gotham_update()].
#' @param threshold Numeric scalar in the interval \eqn{[0, 1)} giving the ICC
#'   boundary of interest.
#' @param alignment Numeric vector of cosine alignments in the interval
#'   \eqn{[-1, 1]}. Defaults to 201 equally spaced values from -1 to 1.
#' @param tolerance Positive numeric scalar used only for numerical decisions
#'   about effectively zero coefficients, discriminants, and duplicate roots.
#'
#' @return A list with:
#' \describe{
#'   \item{threshold}{Requested ICC threshold.}
#'   \item{icc_at_zero}{Updated ICC when the added observation in every group
#'     equals that group's current mean, so \eqn{r=0}.}
#'   \item{icc_at_infinity}{Limiting ICC as \eqn{r} tends to infinity for a
#'     centered proposal ray.}
#'   \item{boundary_ss_ratio}{Required \eqn{SS_B'/SS_W'} ratio at the ICC
#'     threshold.}
#'   \item{frontier}{Data frame containing the alignment, quadratic
#'     coefficients, discriminant, number of finite non-negative roots, first
#'     and second frontier radii, and a status label.}
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
#' gotham_frontier(
#'   state,
#'   threshold = 0.80,
#'   alignment = c(-1, 0, 1)
#' )$frontier
#gotham_frontier <- function(
#    state,
#    threshold,
#    alignment = seq(-1, 1, length.out = 201L),
#    tolerance = sqrt(.Machine$double.eps)) {
#  NULL
#}
gotham_frontier <- function(
    state,
    threshold,
    alignment = seq(-1, 1, length.out = 201L),
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

  if (!is.numeric(alignment) || length(alignment) == 0L ||
      anyNA(alignment) || any(!is.finite(alignment)) ||
      any(alignment < -1) || any(alignment > 1)) {
    stop(
      "alignment must contain finite numeric values in the interval [-1, 1].",
      call. = FALSE
    )
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be a positive finite numeric scalar.", call. = FALSE)
  }

  if (state$group_offset_ss <= tolerance) {
    stop(
      paste0(
        "gotham_frontier() requires nonzero between-group offsets so that ",
        "proposal alignment is defined."
      ),
      call. = FALSE
    )
  }

  g <- state$n_groups
  n <- state$group_size
  n_new <- n + 1L
  offset_norm <- sqrt(state$group_offset_ss)

  # At the threshold, MS_B = k * MS_W. Converting mean squares to
  # sums of squares gives SS_B = boundary_ss_ratio * SS_W.
  k <- (1 + n * threshold) / (1 - threshold)
  boundary_ss_ratio <- (g - 1) * k / (g * n)

  quadratic <- (1 - boundary_ss_ratio * n) / n_new
  constant <- n_new * state$group_offset_ss -
    boundary_ss_ratio * state$ss_within

  frontier_rows <- lapply(
    alignment,
    function(current_alignment) {
      linear <- 2 * offset_norm * current_alignment
      roots <- .gotham_frontier_roots(
        quadratic = quadratic,
        linear = linear,
        constant = constant,
        tolerance = tolerance
      )

      data.frame(
        alignment = current_alignment,
        quadratic = quadratic,
        linear = linear,
        constant = constant,
        discriminant = roots$discriminant,
        n_roots = roots$n_roots,
        radius_1 = roots$radius_1,
        radius_2 = roots$radius_2,
        status = roots$status,
        stringsAsFactors = FALSE
      )
    }
  )

  frontier <- do.call(rbind, frontier_rows)
  rownames(frontier) <- NULL

  zero_ms_between <-
    n_new * state$group_offset_ss / (g - 1L)
  zero_ms_within <- state$ss_within / (g * n)

  icc_at_zero <- .anova_icc_from_components(
    list(
      group_size = n_new,
      ms_between = zero_ms_between,
      ms_within = zero_ms_within
    )
  )

  icc_at_infinity <- 1 / (g + n * (g - 1))

  list(
    threshold = unname(threshold),
    icc_at_zero = unname(icc_at_zero),
    icc_at_infinity = unname(icc_at_infinity),
    boundary_ss_ratio = unname(boundary_ss_ratio),
    frontier = frontier
  )
}


.gotham_frontier_roots <- function(
    quadratic,
    linear,
    constant,
    tolerance) {
  coefficient_scale <- max(
    1,
    abs(quadratic),
    abs(linear),
    abs(constant)
  )
  coefficient_tolerance <- tolerance * coefficient_scale

  if (abs(quadratic) <= coefficient_tolerance) {
    if (abs(linear) <= coefficient_tolerance) {
      if (abs(constant) <= coefficient_tolerance) {
        return(list(
          discriminant = NA_real_,
          n_roots = Inf,
          radius_1 = NA_real_,
          radius_2 = NA_real_,
          status = "all_radii"
        ))
      }

      return(list(
        discriminant = NA_real_,
        n_roots = 0,
        radius_1 = NA_real_,
        radius_2 = NA_real_,
        status = "no_root"
      ))
    }

    root <- -constant / linear
    roots <- .gotham_keep_nonnegative_roots(root, tolerance)

    return(.gotham_pack_frontier_roots(
      roots = roots,
      discriminant = NA_real_
    ))
  }

  discriminant <- linear^2 - 4 * quadratic * constant
  discriminant_scale <- max(
    1,
    linear^2,
    abs(4 * quadratic * constant)
  )
  discriminant_tolerance <- tolerance * discriminant_scale

  if (discriminant < -discriminant_tolerance) {
    return(list(
      discriminant = unname(discriminant),
      n_roots = 0,
      radius_1 = NA_real_,
      radius_2 = NA_real_,
      status = "no_root"
    ))
  }

  if (abs(discriminant) <= discriminant_tolerance) {
    discriminant <- 0
  }

  sqrt_discriminant <- sqrt(discriminant)
  candidate_roots <- c(
    (-linear - sqrt_discriminant) / (2 * quadratic),
    (-linear + sqrt_discriminant) / (2 * quadratic)
  )

  roots <- .gotham_keep_nonnegative_roots(
    candidate_roots,
    tolerance
  )

  .gotham_pack_frontier_roots(
    roots = roots,
    discriminant = discriminant
  )
}


.gotham_keep_nonnegative_roots <- function(roots, tolerance) {
  roots <- as.numeric(roots)
  roots <- roots[is.finite(roots) & roots >= -tolerance]

  if (length(roots) == 0L) {
    return(numeric())
  }

  roots[abs(roots) <= tolerance] <- 0
  roots <- sort(roots)

  if (length(roots) > 1L) {
    keep <- c(
      TRUE,
      diff(roots) > tolerance * pmax(1, abs(roots[-length(roots)]))
    )
    roots <- roots[keep]
  }

  roots
}


.gotham_pack_frontier_roots <- function(roots, discriminant) {
  n_roots <- length(roots)

  list(
    discriminant = unname(discriminant),
    n_roots = n_roots,
    radius_1 = if (n_roots >= 1L) roots[[1L]] else NA_real_,
    radius_2 = if (n_roots >= 2L) roots[[2L]] else NA_real_,
    status = switch(
      as.character(n_roots),
      "0" = "no_root",
      "1" = "one_root",
      "2" = "two_roots",
      "multiple_roots"
    )
  )
}
