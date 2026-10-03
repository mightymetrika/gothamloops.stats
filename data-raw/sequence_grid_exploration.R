# Gotham Loops sequential-grid exploration harness
#
# Analysis-only layer. This file is intentionally kept under data-raw/ and is
# excluded from the built package. It does not modify the package's production
# sequence/frontier functions.
#
# Typical development use:
#   devtools::load_all()
#   source("data-raw/sequence_grid_exploration.R")
#   review <- run_gotham_sequence_grid_review(state, steps = 50)
#
# The harness advances each fixed alignment one step at a time by calling
# gotham_sequence(..., steps = 1). This preserves completed history if a later
# step loses the first frontier root or otherwise becomes inadmissible.

run_gotham_sequence_grid_review <- function(
    state,
    steps = 50L,
    threshold = 0.80,
    alignment = seq(-1, 1, by = 0.25),
    orientation = 1,
    tolerance = sqrt(.Machine$double.eps)) {
  required_functions <- c(
    "gotham_frontier",
    "gotham_sequence"
  )

  missing_functions <- required_functions[
    !vapply(
      required_functions,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE
    )
  ]

  if (length(missing_functions) > 0L) {
    stop(
      paste0(
        "Load gothamloops.stats before running this harness. Missing: ",
        paste(missing_functions, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (!is.list(state) || is.null(state$icc) || is.null(state$group_size) ||
      is.null(state$group_offset_ss) || is.null(state$ss_within)) {
    stop("state must be a Gotham state object.", call. = FALSE)
  }

  if (!is.numeric(steps) || length(steps) != 1L || is.na(steps) ||
      !is.finite(steps) || steps < 1 || steps != as.integer(steps)) {
    stop("steps must be a positive integer scalar.", call. = FALSE)
  }
  steps <- as.integer(steps)

  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
      !is.finite(threshold) || threshold < 0 || threshold >= 1) {
    stop(
      "threshold must be a finite numeric scalar in the interval [0, 1).",
      call. = FALSE
    )
  }

  if (!is.numeric(alignment) || length(alignment) == 0L || anyNA(alignment) ||
      any(!is.finite(alignment)) || any(alignment < -1) || any(alignment > 1)) {
    stop(
      "alignment must contain finite numeric values in the interval [-1, 1].",
      call. = FALSE
    )
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

  alignment <- unique(as.numeric(alignment))

  alignment_results <- lapply(
    alignment,
    function(current_alignment) {
      .run_single_alignment_review(
        state = state,
        steps = steps,
        threshold = threshold,
        alignment = current_alignment,
        orientation = orientation,
        tolerance = tolerance
      )
    }
  )

  trajectory <- .bind_review_rows(
    lapply(alignment_results, `[[`, "trajectory")
  )
  frontier_diagnostics <- .bind_review_rows(
    lapply(alignment_results, `[[`, "frontier_diagnostics")
  )
  summary <- do.call(
    rbind,
    lapply(alignment_results, `[[`, "summary")
  )
  rownames(summary) <- NULL

  list(
    settings = list(
      requested_steps = steps,
      threshold = unname(threshold),
      alignment = alignment,
      orientation = orientation,
      tolerance = tolerance
    ),
    summary = summary,
    trajectory = trajectory,
    frontier_diagnostics = frontier_diagnostics,
    initial_state = state
  )
}


.run_single_alignment_review <- function(
    state,
    steps,
    threshold,
    alignment,
    orientation,
    tolerance) {
  current_state <- state
  trajectory_rows <- list()
  diagnostic_rows <- list()
  status <- "complete"
  termination_step <- NA_integer_
  termination_reason <- NA_character_

  for (step in seq_len(steps)) {
    frontier <- tryCatch(
      gotham_frontier(
        state = current_state,
        threshold = threshold,
        alignment = alignment,
        tolerance = tolerance
      ),
      error = identity
    )

    if (inherits(frontier, "error")) {
      status <- "frontier_error"
      termination_step <- step
      termination_reason <- conditionMessage(frontier)
      break
    }

    frontier_row <- frontier$frontier[1L, , drop = FALSE]

    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
      alignment = alignment,
      step = step,
      group_size = current_state$group_size,
      icc_before = current_state$icc,
      icc_at_zero = frontier$icc_at_zero,
      icc_at_infinity = frontier$icc_at_infinity,
      n_roots = frontier_row$n_roots[[1L]],
      radius_1 = frontier_row$radius_1[[1L]],
      radius_2 = frontier_row$radius_2[[1L]],
      frontier_status = frontier_row$status[[1L]],
      stringsAsFactors = FALSE
    )

    if (frontier$icc_at_zero < threshold - tolerance) {
      status <- "zero_below_threshold"
      termination_step <- step
      termination_reason <- paste0(
        "Zero-radius ICC (",
        format(frontier$icc_at_zero, digits = 10),
        ") is below threshold (",
        format(threshold, digits = 10),
        ")."
      )
      break
    }

    if (!is.finite(frontier_row$radius_1[[1L]])) {
      status <- "no_first_root"
      termination_step <- step
      termination_reason <- paste0(
        "No finite first frontier radius is available; frontier status = ",
        frontier_row$status[[1L]],
        "."
      )
      break
    }

    one_step <- tryCatch(
      gotham_sequence(
        state = current_state,
        steps = 1L,
        threshold = threshold,
        alignment = alignment,
        orientation = orientation,
        tolerance = tolerance
      ),
      error = identity
    )

    if (inherits(one_step, "error")) {
      status <- "sequence_error"
      termination_step <- step
      termination_reason <- conditionMessage(one_step)
      break
    }

    row <- one_step$trajectory
    # gotham_sequence(..., steps = 1) restarts its local step counter at 1.
    # Preserve the outer review step so trajectories span 1, 2, ..., steps.
    row$step <- step
    row$alignment_fixed <- alignment
    row$frontier_n_roots <- frontier_row$n_roots[[1L]]
    row$frontier_radius_2 <- frontier_row$radius_2[[1L]]
    row$frontier_status <- frontier_row$status[[1L]]
    trajectory_rows[[length(trajectory_rows) + 1L]] <- row

    current_state <- one_step$final_state
  }

  trajectory <- .bind_review_rows(trajectory_rows)
  frontier_diagnostics <- .bind_review_rows(diagnostic_rows)
  completed_steps <- nrow(trajectory)

  if (completed_steps > 0L) {
    first_radius <- trajectory$radius[[1L]]
    last_radius <- trajectory$radius[[completed_steps]]
    min_radius <- min(trajectory$radius)
    max_radius <- max(trajectory$radius)
    min_icc_at_zero <- min(trajectory$icc_at_zero)
    final_icc <- trajectory$icc_after[[completed_steps]]
  } else {
    first_radius <- NA_real_
    last_radius <- NA_real_
    min_radius <- NA_real_
    max_radius <- NA_real_
    min_icc_at_zero <- NA_real_
    final_icc <- state$icc
  }

  any_two_roots <- if (nrow(frontier_diagnostics) > 0L) {
    any(frontier_diagnostics$n_roots >= 2, na.rm = TRUE)
  } else {
    FALSE
  }

  first_two_root_step <- if (any_two_roots) {
    min(frontier_diagnostics$step[
      frontier_diagnostics$n_roots >= 2
    ])
  } else {
    NA_integer_
  }

  summary <- data.frame(
    alignment = alignment,
    requested_steps = steps,
    completed_steps = completed_steps,
    status = status,
    termination_step = termination_step,
    termination_reason = termination_reason,
    initial_icc = state$icc,
    final_icc = final_icc,
    first_radius = first_radius,
    last_radius = last_radius,
    radius_ratio = if (is.finite(first_radius) && first_radius != 0) {
      last_radius / first_radius
    } else {
      NA_real_
    },
    min_radius = min_radius,
    max_radius = max_radius,
    initial_group_offset_ss = state$group_offset_ss,
    final_group_offset_ss = current_state$group_offset_ss,
    group_offset_ss_ratio = current_state$group_offset_ss /
      state$group_offset_ss,
    initial_ss_within = state$ss_within,
    final_ss_within = current_state$ss_within,
    ss_within_ratio = if (state$ss_within != 0) {
      current_state$ss_within / state$ss_within
    } else {
      NA_real_
    },
    min_icc_at_zero = min_icc_at_zero,
    any_two_roots = any_two_roots,
    first_two_root_step = first_two_root_step,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    trajectory = trajectory,
    frontier_diagnostics = frontier_diagnostics,
    final_state = current_state
  )
}


.bind_review_rows <- function(rows) {
  rows <- rows[
    vapply(rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))
  ]

  if (length(rows) == 0L) {
    return(data.frame())
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


plot_gotham_sequence_grid <- function(
    review,
    metric = c("radius", "group_offset_ss", "ss_within", "icc_at_zero"),
    log_y = FALSE,
    relative = FALSE) {
  metric <- match.arg(metric)

  if (!is.list(review) || is.null(review$trajectory) ||
      !is.data.frame(review$trajectory) || nrow(review$trajectory) == 0L) {
    stop("review must contain a non-empty trajectory table.", call. = FALSE)
  }

  trajectory <- review$trajectory
  value_column <- switch(
    metric,
    radius = "radius",
    group_offset_ss = "group_offset_ss_after",
    ss_within = "ss_within_after",
    icc_at_zero = "icc_at_zero"
  )

  y_label <- switch(
    metric,
    radius = "Frontier radius",
    group_offset_ss = "Group offset SS",
    ss_within = "Within-group SS",
    icc_at_zero = "ICC at zero-radius proposal"
  )

  plot_values <- trajectory[[value_column]]

  if (isTRUE(relative)) {
    if (metric == "radius") {
      first_radius_by_alignment <- tapply(
        trajectory$radius,
        trajectory$alignment_fixed,
        function(x) x[[1L]]
      )
      denominators <- unname(
        first_radius_by_alignment[
          as.character(trajectory$alignment_fixed)
        ]
      )
      plot_values <- trajectory$radius / denominators
      y_label <- "Frontier radius / step 1 radius"
    } else if (metric == "group_offset_ss") {
      denominator <- review$initial_state$group_offset_ss
      if (!is.finite(denominator) || denominator <= 0) {
        stop(
          "Relative group_offset_ss plotting requires a positive initial value.",
          call. = FALSE
        )
      }
      plot_values <- trajectory$group_offset_ss_after / denominator
      y_label <- "Group offset SS / initial group offset SS"
    } else if (metric == "ss_within") {
      denominator <- review$initial_state$ss_within
      if (!is.finite(denominator) || denominator <= 0) {
        stop(
          "Relative ss_within plotting requires a positive initial value.",
          call. = FALSE
        )
      }
      plot_values <- trajectory$ss_within_after / denominator
      y_label <- "Within-group SS / initial within-group SS"
    } else {
      stop(
        "relative = TRUE is not defined for metric = 'icc_at_zero'.",
        call. = FALSE
      )
    }
  }

  trajectory$.plot_value <- plot_values

  alignment_values <- sort(unique(trajectory$alignment_fixed))
  pieces <- lapply(
    alignment_values,
    function(current_alignment) {
      trajectory[
        trajectory$alignment_fixed == current_alignment,
        c("step", ".plot_value"),
        drop = FALSE
      ]
    }
  )

  max_step <- max(trajectory$step)
  y_values <- trajectory$.plot_value

  if (isTRUE(log_y) && any(y_values <= 0, na.rm = TRUE)) {
    stop("log_y = TRUE requires strictly positive plotted values.", call. = FALSE)
  }

  line_colours <- grDevices::hcl.colors(
    length(alignment_values),
    palette = "Dark 3"
  )
  line_types <- ((seq_along(alignment_values) - 1L) %% 6L) + 1L

  graphics::plot(
    NA,
    xlim = c(1, max_step),
    ylim = range(y_values, finite = TRUE),
    log = if (isTRUE(log_y)) "y" else "",
    xlab = "Step",
    ylab = y_label,
    main = paste("Gotham sequence grid:", y_label)
  )

  if (isTRUE(relative)) {
    graphics::abline(h = 1, lty = 3, col = "grey50")
  }

  for (index in seq_along(pieces)) {
    graphics::lines(
      pieces[[index]]$step,
      pieces[[index]]$.plot_value,
      lty = line_types[[index]],
      lwd = 2,
      col = line_colours[[index]]
    )
  }

  graphics::legend(
    "topleft",
    legend = paste0("c = ", format(alignment_values, trim = TRUE)),
    lty = line_types,
    lwd = 2,
    col = line_colours,
    bty = "n",
    cex = 0.8,
    ncol = 2
  )

  invisible(review)
}


plot_gotham_sequence_grid_endpoint <- function(
    review,
    metric = c("radius_ratio", "group_offset_ss_ratio", "ss_within_ratio"),
    log_y = TRUE) {
  metric <- match.arg(metric)

  if (!is.list(review) || is.null(review$summary) ||
      !is.data.frame(review$summary) || nrow(review$summary) == 0L) {
    stop("review must contain a non-empty summary table.", call. = FALSE)
  }

  summary <- review$summary
  y_values <- summary[[metric]]

  if (isTRUE(log_y) && any(y_values <= 0, na.rm = TRUE)) {
    stop("log_y = TRUE requires strictly positive plotted values.", call. = FALSE)
  }

  y_label <- switch(
    metric,
    radius_ratio = "Last radius / first radius",
    group_offset_ss_ratio = "Final / initial group offset SS",
    ss_within_ratio = "Final / initial within-group SS"
  )

  point_colours <- grDevices::hcl.colors(
    nrow(summary),
    palette = "Dark 3"
  )

  graphics::plot(
    summary$alignment,
    y_values,
    type = "b",
    pch = 19,
    lwd = 2,
    col = "grey35",
    log = if (isTRUE(log_y)) "y" else "",
    xlab = "Fixed alignment",
    ylab = y_label,
    main = paste("50-step Gotham endpoint:", y_label)
  )

  graphics::points(
    summary$alignment,
    y_values,
    pch = 19,
    col = point_colours
  )

  graphics::abline(h = 1, lty = 3, col = "grey50")

  graphics::text(
    summary$alignment,
    y_values,
    labels = format(summary$alignment, trim = TRUE),
    pos = 3,
    cex = 0.75
  )

  invisible(review)
}
