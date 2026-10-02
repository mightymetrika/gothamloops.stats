#' Balanced one-way ANOVA intraclass correlation
#'
#' Computes the single-measure one-way random-effects ANOVA intraclass
#' correlation coefficient (ICC) for balanced grouped data.
#'
#' The estimator is
#' \deqn{ICC = (MS_B - MS_W) / (MS_B + (n - 1) MS_W),}
#' where \eqn{MS_B} is the between-group mean square, \eqn{MS_W} is the
#' within-group mean square, and \eqn{n} is the common group size.
#'
#' The estimate is not truncated to the interval \eqn{[0, 1]}. Negative estimates
#' are retained because they can arise from the ANOVA estimator when the
#' observed between-group variation is smaller than the within-group variation.
#'
#' @param data A data frame containing the outcome and grouping variable.
#' @param outcome Character scalar naming the numeric outcome column.
#' @param group Character scalar naming the grouping column.
#'
#' @return A numeric scalar containing the ANOVA ICC estimate.
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 3),
#'   y = 1:9
#' )
#'
#' anova_icc(dat, outcome = "y", group = "group")
anova_icc <- function(data, outcome, group) {
  components <- .anova_icc_components(
    data = data,
    outcome = outcome,
    group = group
  )

  .anova_icc_from_components(components)
}


.anova_icc_from_components <- function(components) {
  denominator <- components$ms_between +
    (components$group_size - 1) * components$ms_within

  if (!is.finite(denominator) || denominator <= 0) {
    stop(
      "The ANOVA ICC is undefined because its denominator is not positive.",
      call. = FALSE
    )
  }

  unname(
    (components$ms_between - components$ms_within) /
      denominator
  )
}


.anova_icc_components <- function(data, outcome, group) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }

  if (!is.character(outcome) || length(outcome) != 1L ||
      is.na(outcome) || !nzchar(outcome)) {
    stop("outcome must be a single non-empty column name.", call. = FALSE)
  }

  if (!is.character(group) || length(group) != 1L ||
      is.na(group) || !nzchar(group)) {
    stop("group must be a single non-empty column name.", call. = FALSE)
  }

  if (!outcome %in% names(data)) {
    stop(
      paste0("Outcome column not found: ", outcome, "."),
      call. = FALSE
    )
  }

  if (!group %in% names(data)) {
    stop(
      paste0("Grouping column not found: ", group, "."),
      call. = FALSE
    )
  }

  y <- data[[outcome]]
  g <- data[[group]]

  if (!is.numeric(y)) {
    stop("The outcome column must be numeric.", call. = FALSE)
  }

  if (length(y) == 0L) {
    stop("data must contain at least one row.", call. = FALSE)
  }

  if (anyNA(y) || any(!is.finite(y))) {
    stop("The outcome column must contain only finite, non-missing values.", call. = FALSE)
  }

  if (anyNA(g)) {
    stop("The grouping column must not contain missing values.", call. = FALSE)
  }

  g <- factor(g)
  group_sizes <- table(g)
  n_groups <- length(group_sizes)

  if (n_groups < 2L) {
    stop("At least two groups are required.", call. = FALSE)
  }

  if (length(unique(as.integer(group_sizes))) != 1L) {
    stop(
      "anova_icc() currently requires a balanced design with equal group sizes.",
      call. = FALSE
    )
  }

  group_size <- unname(as.integer(group_sizes[[1L]]))

  if (group_size < 2L) {
    stop("Each group must contain at least two observations.", call. = FALSE)
  }

  group_means <- tapply(y, g, mean)
  grand_mean <- mean(y)

  ss_between <- group_size * sum((group_means - grand_mean)^2)
  fitted_group_mean <- group_means[as.integer(g)]
  ss_within <- sum((y - fitted_group_mean)^2)

  df_between <- n_groups - 1L
  df_within <- n_groups * (group_size - 1L)
  ms_between <- ss_between / df_between
  ms_within <- ss_within / df_within

  list(
    n_groups = n_groups,
    group_size = group_size,
    group_labels = names(group_means),
    grand_mean = grand_mean,
    group_means = as.numeric(group_means),
    ss_between = unname(ss_between),
    ss_within = unname(ss_within),
    df_between = df_between,
    df_within = df_within,
    ms_between = unname(ms_between),
    ms_within = unname(ms_within)
  )
}
