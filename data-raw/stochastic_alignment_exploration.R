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


# Neutral-boundary two-state Markov review ---------------------------------

.gotham_two_state_markov_log_moments <- function(
    negative_alignment,
    positive_alignment,
    stationary_positive_share,
    rho,
    initial_n,
    steps,
    g,
    tau) {
  p <- stationary_positive_share
  n_values <- initial_n + seq_len(steps) - 1L

  log_negative <- log(
    gotham_normalized_multiplier(
      n = n_values,
      g = g,
      tau = tau,
      alignment = negative_alignment
    )
  )

  log_positive <- log(
    gotham_normalized_multiplier(
      n = n_values,
      g = g,
      tau = tau,
      alignment = positive_alignment
    )
  )

  expected_increment <-
    (1 - p) * log_negative +
    p * log_positive

  delta <-
    log_positive -
    log_negative

  # For a stationary two-state chain, if I_n is the indicator of the
  # positive state then
  #
  #   Cov(I_n, I_{n+k}) = p(1-p) rho^k.
  #
  # Therefore the exact finite-horizon variance of the random part of
  # sum log m_n(C_n) can be accumulated in O(steps) time.
  lagged_delta <- 0
  quadratic_form <- 0

  for (delta_n in delta) {
    quadratic_form <-
      quadratic_form +
      delta_n^2 +
      2 * delta_n * lagged_delta

    lagged_delta <-
      rho * (
        lagged_delta +
        delta_n
      )
  }

  list(
    mean_log_ratio =
      sum(expected_increment),
    variance_log_ratio =
      p * (1 - p) *
      quadratic_form,
    sd_log_ratio =
      sqrt(
        p * (1 - p) *
          quadratic_form
      ),
    delta = delta,
    expected_increment =
      expected_increment
  )
}

run_gotham_neutral_markov_boundary_review <- function(
    negative_alignment = -0.5,
    positive_alignment = 0.5,
    rho = c(-0.5, 0, 0.5, 0.8),
    initial_n,
    steps = 5000L,
    reps = 800L,
    g,
    tau,
    seed = 20261009L) {
  if (
    negative_alignment >= 0 ||
      positive_alignment <= 0
  ) {
    stop(
      "Require negative_alignment < 0 < positive_alignment.",
      call. = FALSE
    )
  }

  h_negative <- gotham_asymptotic_drift(
    g = g,
    tau = tau,
    alignment = negative_alignment
  )

  h_positive <- gotham_asymptotic_drift(
    g = g,
    tau = tau,
    alignment = positive_alignment
  )

  p <- -h_negative /
    (
      h_positive -
        h_negative
    )

  lower_rho <- -min(
    p / (1 - p),
    (1 - p) / p
  )

  if (
    any(rho <= lower_rho) ||
      any(rho >= 1)
  ) {
    stop(
      "rho values must define an irreducible stationary two-state chain. ",
      "For this neutral marginal, require rho > ",
      format(lower_rho, digits = 8),
      " and rho < 1.",
      call. = FALSE
    )
  }

  steps <- as.integer(steps)
  reps <- as.integer(reps)

  n_values <- initial_n + seq_len(steps) - 1L

  log_negative <- log(
    gotham_normalized_multiplier(
      n = n_values,
      g = g,
      tau = tau,
      alignment = negative_alignment
    )
  )

  log_positive <- log(
    gotham_normalized_multiplier(
      n = n_values,
      g = g,
      tau = tau,
      alignment = positive_alignment
    )
  )

  set.seed(seed)

  summaries <- lapply(
    rho,
    function(rho_value) {
      p11 <-
        p +
        rho_value * (1 - p)

      p01 <-
        p * (1 - rho_value)

      theory <- .gotham_two_state_markov_log_moments(
        negative_alignment = negative_alignment,
        positive_alignment = positive_alignment,
        stationary_positive_share = p,
        rho = rho_value,
        initial_n = initial_n,
        steps = steps,
        g = g,
        tau = tau
      )

      state_positive <-
        stats::runif(reps) < p

      log_ratio <- numeric(reps)

      for (step_index in seq_len(steps)) {
        log_ratio <- log_ratio +
          ifelse(
            state_positive,
            log_positive[[step_index]],
            log_negative[[step_index]]
          )

        u <- stats::runif(reps)

        state_positive <- ifelse(
          state_positive,
          u < p11,
          u < p01
        )
      }

      data.frame(
        rho = rho_value,
        stationary_positive_share = p,
        mean_alignment =
          (1 - p) * negative_alignment +
          p * positive_alignment,
        asymptotic_drift =
          (1 - p) * h_negative +
          p * h_positive,
        theoretical_mean_log_ratio =
          theory$mean_log_ratio,
        empirical_mean_log_ratio =
          mean(log_ratio),
        theoretical_sd_log_ratio =
          theory$sd_log_ratio,
        empirical_sd_log_ratio =
          stats::sd(log_ratio),
        median_final_A_ratio =
          stats::median(
            exp(log_ratio)
          ),
        q025_final_A_ratio =
          unname(
            stats::quantile(
              exp(log_ratio),
              0.025
            )
          ),
        q975_final_A_ratio =
          unname(
            stats::quantile(
              exp(log_ratio),
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
