#' Summarize the current balanced ANOVA grouping state
#'
#' Constructs the core quantities used by the statistical Gotham Loops
#' framework for a balanced one-way grouped outcome. The returned state is a
#' descriptive snapshot of the observed grouping structure; it does not add or
#' simulate any observations.
#'
#' Let \eqn{\bar{y}_g} denote the mean for group \eqn{g} and \eqn{\bar{y}}
#' the grand mean. The group-offset vector is
#' \deqn{a_g = \bar{y}_g - \bar{y},}
#' and `group_offset_ss` is \eqn{\sum_g a_g^2}. For a balanced design with
#' common group size \eqn{n},
#' \deqn{SS_B = n \sum_g a_g^2.}
#'
#' The ANOVA ICC is the same single-measure one-way random-effects estimator
#' returned by [anova_icc()].
#'
#' @param data A data frame containing the outcome and grouping variable.
#' @param outcome Character scalar naming the numeric outcome column.
#' @param group Character scalar naming the grouping column.
#'
#' @return A named list with the current grouping state:
#' \describe{
#'   \item{n_groups}{Number of groups.}
#'   \item{group_size}{Common number of observations per group.}
#'   \item{grand_mean}{Overall outcome mean.}
#'   \item{group_means}{Named vector of group means.}
#'   \item{group_offsets}{Named vector of group means minus the grand mean.}
#'   \item{group_offset_ss}{Sum of squared group offsets.}
#'   \item{ss_between}{Between-group sum of squares.}
#'   \item{ss_within}{Within-group sum of squares.}
#'   \item{ms_between}{Between-group mean square.}
#'   \item{ms_within}{Within-group mean square.}
#'   \item{icc}{Balanced one-way ANOVA ICC.}
#' }
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 3),
#'   y = 1:9
#' )
#'
#' gotham_state(dat, outcome = "y", group = "group")
gotham_state <- function(data, outcome, group) {
  components <- .anova_icc_components(
    data = data,
    outcome = outcome,
    group = group
  )

  group_means <- stats::setNames(
    components$group_means,
    components$group_labels
  )
  group_offsets <- group_means - components$grand_mean

  list(
    n_groups = components$n_groups,
    group_size = components$group_size,
    grand_mean = components$grand_mean,
    group_means = group_means,
    group_offsets = group_offsets,
    group_offset_ss = unname(sum(group_offsets^2)),
    ss_between = components$ss_between,
    ss_within = components$ss_within,
    ms_between = components$ms_between,
    ms_within = components$ms_within,
    icc = .anova_icc_from_components(components)
  )
}
