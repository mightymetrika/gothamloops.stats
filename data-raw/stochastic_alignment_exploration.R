# Gotham Loops stochastic-alignment exploration harness
#
# Analysis-only layer. This file is intentionally kept under data-raw/ and is
# excluded from the built package.
#
# It studies random alignment sequences using the exact normalized recurrence
# derived in the Gotham hypothesis notebook.

.gotham_z_infinity <- function(
    g,
    tau,
    alignment) {
  a_infinity <- tau * (g - 1) /
    (
      g * (1 - tau)
    )

  (
    alignment +
      sqrt(
        alignment^2 +
          a_infinity
      )
  ) / a_infinity
}

gotham_asymptotic_drift <- function(
    g,
    tau,
    alignment) {
  2 * alignment *
    .gotham_z_infinity(
      g = g,
      tau = tau,
      alignment = alignment
    )
}

.gotham_normalized_radius <- function(
    n,
    g,
    tau,
    alignment) {
  d_n <- g + n * (g - 1)

  a_n <- (
    tau * d_n - 1
  ) / (
    g * (n + 1) * (1 - tau)
  )

  gamma_n <- (
    tau * (n - 1) + 2
  ) / (
    tau * (n - 1) + 1
  )

  if (any(a_n <= 0)) {
    stop(
      "The normalized frontier formula requires the H2 regime with a_n > 0.",
      call. = FALSE
    )
  }

  (
    alignment +
      sqrt(
        alignment^2 +
          a_n * gamma_n
      )
  ) / a_n
}

gotham_normalized_multiplier <- function(
    n,
    g,
    tau,
    alignment) {
  z <- .gotham_normalized_radius(
    n = n,
    g = g,
    tau = tau,
    alignment = alignment
  )

  1 +
    2 * alignment * z / (n + 1) +
    z^2 / (n + 1)^2
}

.validate_alignment_distribution <- function(
    distribution,
    name = "distribution") {
  if (
    !is.list(distribution) ||
      is.null(distribution$values) ||
      is.null(distribution$prob)
  ) {
    stop(
      name,
      " must be a list with numeric 'values' and 'prob' components.",
      call. = FALSE
    )
  }

  values <- distribution$values
  prob <- distribution$prob

  if (
    !is.numeric(values) ||
      length(values) == 0L ||
      anyNA(values) ||
      any(!is.finite(values)) ||
      any(values < -1) ||
      any(values > 1)
  ) {
    stop(
      name,
      "$values must be finite numeric alignments in [-1, 1].",
      call. = FALSE
    )
  }

  if (
    !is.numeric(prob) ||
      length(prob) != length(values) ||
      anyNA(prob) ||
      any(!is.finite(prob)) ||
      any(prob < 0) ||
      sum(prob) <= 0
  ) {
    stop(
      name,
      "$prob must be non-negative finite weights matching values.",
      call. = FALSE
    )
  }

  prob <- prob / sum(prob)

  list(
    values = as.numeric(values),
    prob = as.numeric(prob)
  )
}

run_gotham_iid_alignment_review <- function(
    distributions,
    initial_n,
    steps = 5000L,
    reps = 400L,
    g,
    tau,
    seed = 20261005L) {
  if (
    !is.list(distributions) ||
      length(distributions) == 0L ||
      is.null(names(distributions)) ||
      any(names(distributions) == "")
  ) {
    stop(
      "distributions must be a named, non-empty list.",
      call. = FALSE
    )
  }

  steps <- as.integer(steps)
  reps <- as.integer(reps)

  if (
    steps < 1L ||
      reps < 1L
  ) {
    stop(
      "steps and reps must be positive integers.",
      call. = FALSE
    )
  }

  n_values <- initial_n + seq_len(steps) - 1L
  log_scale <- log(
    (initial_n + steps) /
      initial_n
  )

  set.seed(seed)

  summaries <- lapply(
    names(distributions),
    function(scenario_name) {
      distribution <- .validate_alignment_distribution(
        distributions[[scenario_name]],
        scenario_name
      )

      values <- distribution$values
      prob <- distribution$prob

      mean_alignment <- sum(
        values * prob
      )

      asymptotic_drift <- sum(
        gotham_asymptotic_drift(
          g = g,
          tau = tau,
          alignment = values
        ) * prob
      )

      finite_expected_log_growth <- sum(
        vapply(
          n_values,
          function(n_value) {
            sum(
              prob *
                log(
                  gotham_normalized_multiplier(
                    n = n_value,
                    g = g,
                    tau = tau,
                    alignment = values
                  )
                )
            )
          },
          numeric(1)
        )
      )

      finite_expected_exponent <-
        finite_expected_log_growth /
        log_scale

      log_growth <- numeric(reps)

      for (n_value in n_values) {
        draws <- sample(
          values,
          size = reps,
          replace = TRUE,
          prob = prob
        )

        log_growth <- log_growth +
          log(
            gotham_normalized_multiplier(
              n = n_value,
              g = g,
              tau = tau,
              alignment = draws
            )
          )
      }

      empirical_exponent <-
        log_growth /
        log_scale

      data.frame(
        scenario = scenario_name,
        mean_alignment = mean_alignment,
        asymptotic_drift = asymptotic_drift,
        finite_expected_exponent =
          finite_expected_exponent,
        empirical_mean_exponent =
          mean(empirical_exponent),
        empirical_sd_exponent =
          stats::sd(empirical_exponent),
        q025_exponent =
          unname(
            stats::quantile(
              empirical_exponent,
              0.025
            )
          ),
        q975_exponent =
          unname(
            stats::quantile(
              empirical_exponent,
              0.975
            )
          ),
        median_final_A_ratio =
          stats::median(
            exp(log_growth)
          ),
        probability_final_expansion =
          mean(log_growth > 0),
        stringsAsFactors = FALSE
      )
    }
  )

  out <- do.call(
    rbind,
    summaries
  )

  rownames(out) <- NULL
  out
}

run_gotham_markov_alignment_review <- function(
    negative_alignment = -0.5,
    positive_alignment = 0.5,
    stationary_positive_share = 0.4,
    rho = c(-0.5, 0, 0.8),
    initial_n,
    steps = 5000L,
    reps = 400L,
    g,
    tau,
    seed = 20261006L) {
  if (
    negative_alignment >= 0 ||
      positive_alignment <= 0
  ) {
    stop(
      "Require negative_alignment < 0 < positive_alignment.",
      call. = FALSE
    )
  }

  p <- stationary_positive_share

  if (
    !is.numeric(p) ||
      length(p) != 1L ||
      !is.finite(p) ||
      p <= 0 ||
      p >= 1
  ) {
    stop(
      "stationary_positive_share must lie strictly between 0 and 1.",
      call. = FALSE
    )
  }

  rho <- as.numeric(rho)

  transition_ok <- vapply(
    rho,
    function(rho_value) {
      p11 <- p + rho_value * (1 - p)
      p01 <- p * (1 - rho_value)

      p11 >= 0 &&
        p11 <= 1 &&
        p01 >= 0 &&
        p01 <= 1
    },
    logical(1)
  )

  if (!all(transition_ok)) {
    stop(
      "At least one rho value implies invalid two-state transition probabilities.",
      call. = FALSE
    )
  }

  steps <- as.integer(steps)
  reps <- as.integer(reps)

  n_values <- initial_n + seq_len(steps) - 1L
  log_scale <- log(
    (initial_n + steps) /
      initial_n
  )

  values <- c(
    negative_alignment,
    positive_alignment
  )

  prob <- c(
    1 - p,
    p
  )

  mean_alignment <- sum(
    values * prob
  )

  asymptotic_drift <- sum(
    gotham_asymptotic_drift(
      g = g,
      tau = tau,
      alignment = values
    ) * prob
  )

  finite_expected_log_growth <- sum(
    vapply(
      n_values,
      function(n_value) {
        sum(
          prob *
            log(
              gotham_normalized_multiplier(
                n = n_value,
                g = g,
                tau = tau,
                alignment = values
              )
            )
        )
      },
      numeric(1)
    )
  )

  finite_expected_exponent <-
    finite_expected_log_growth /
    log_scale

  set.seed(seed)

  summaries <- lapply(
    rho,
    function(rho_value) {
      p11 <- p + rho_value * (1 - p)
      p01 <- p * (1 - rho_value)

      state_positive <-
        stats::runif(reps) < p

      log_growth <- numeric(reps)

      for (n_value in n_values) {
        alignments <- ifelse(
          state_positive,
          positive_alignment,
          negative_alignment
        )

        log_growth <- log_growth +
          log(
            gotham_normalized_multiplier(
              n = n_value,
              g = g,
              tau = tau,
              alignment = alignments
            )
          )

        u <- stats::runif(reps)

        state_positive <- ifelse(
          state_positive,
          u < p11,
          u < p01
        )
      }

      empirical_exponent <-
        log_growth /
        log_scale

      data.frame(
        rho = rho_value,
        stationary_positive_share = p,
        mean_alignment = mean_alignment,
        asymptotic_drift = asymptotic_drift,
        finite_expected_exponent =
          finite_expected_exponent,
        empirical_mean_exponent =
          mean(empirical_exponent),
        empirical_sd_exponent =
          stats::sd(empirical_exponent),
        sd_final_log_A_ratio =
          stats::sd(log_growth),
        q025_exponent =
          unname(
            stats::quantile(
              empirical_exponent,
              0.025
            )
          ),
        q975_exponent =
          unname(
            stats::quantile(
              empirical_exponent,
              0.975
            )
          ),
        stringsAsFactors = FALSE
      )
    }
  )

  out <- do.call(
    rbind,
    summaries
  )

  rownames(out) <- NULL
  out
}
