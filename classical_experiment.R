# Classical fixed-design experiments for the results section.
#
# The script uses independent calibration and inference halves for every
# method. Credible and Oracle use a direct 2.5% failure level. The OLS split
# interval divides that total budget equally between calibration failure
# (delta) and inference noise (alpha_inf), so
#
#   delta + alpha_inf = 0.025.

library(tidyverse)

# set.seed(20260925)
set.seed(1000)

standardize_columns <- function(X) {
  # Center each feature and rescale it to have Euclidean norm sqrt(n).
  X <- as.matrix(X)
  X <- sweep(X, 2, colMeans(X), "-")
  norms <- sqrt(colSums(X * X))
  if (any(!is.finite(norms)) || any(norms == 0)) {
    stop("Design contains a constant or non-finite column")
  }
  sweep(X, 2, norms, "/") * sqrt(nrow(X))
}

# Generate covariance matrices
augmentation_penalty <- function(X, kind, strength) {
  n <- nrow(X)
  p <- ncol(X)
  if (kind == "Dropout") {
    beta <- strength
    return(beta / (1 - beta) * diag(diag(crossprod(X)), p, p))
  }
  if (kind == "Gaussian") {
    tau <- strength
    return(n * tau^2 * diag(p))
  }
}

load_designs <- function() {
  wine_raw <- read.csv("data/wine.data", header = FALSE)
  wine <- standardize_columns(wine_raw[, 2:14])

  cancer_raw <- read.csv("data/wdbc.data", header = FALSE)
  cancer_standardized <- standardize_columns(cancer_raw[, 3:32])
  cancer_svd <- svd(cancer_standardized, nu = 0)
  keep <- which((cancer_svd$d[1] / cancer_svd$d)^2 <= 1e3)
  cancer <- cancer_standardized %*% cancer_svd$v[, keep, drop = FALSE]

  auto_raw <- read.table(
    "data/auto-mpg.data",
    header = FALSE,
    na.strings = "?",
    quote = '"',
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  auto_raw <- auto_raw[complete.cases(auto_raw[, 1:8]), ]
  auto <- standardize_columns(auto_raw[, 2:8])

  gaussian <- standardize_columns(matrix(rnorm(420 * 20), nrow = 420, ncol = 20))

  list(
    "Gaussian design" = gaussian,
    "UCI Wine" = wine,
    "UCI Breast Cancer" = cancer,
    "UCI Auto MPG" = auto
  )
}

# Calculate target directions
target_directions <- function(X, theta) {
  eig <- eigen(crossprod(X), symmetric = TRUE) # get eigenvectors of X^T * X
  list(
    signal = theta / sqrt(sum(theta * theta)), # signal direction
    leading = eig$vectors[, 1], # leading eigenvector
    trailing = eig$vectors[, ncol(eig$vectors)] # trailing eigenvector
  )
}


# Create list of target direction, gamma (target functional), g (bias direction),
# sampling and posterior standard deviation, bias and matrix product of A^(-1) * X^T
prepare_direction <- function(X, theta, sigma, v, lambda_aug) {
  p <- ncol(X)
  XtX <- crossprod(X)
  A <- XtX + lambda_aug
  A_inv <- solve(A)
  S <- A_inv %*% XtX
  g <- as.vector(crossprod(S - diag(p), v))
  list(
    v = v,
    gamma = as.numeric(crossprod(v, theta)),
    g = g,
    sampling_sd = sigma * sqrt(as.numeric(crossprod(v, A_inv %*% XtX %*% A_inv %*% v))),
    posterior_sd = sigma * sqrt(as.numeric(crossprod(v, A_inv %*% v))),
    exact_bias = abs(as.numeric(crossprod(g, theta))),
    A_inv_Xt = A_inv %*% t(X)
  )
}


# Run experiments
alpha_total <- 0.025
alpha_baseline <- alpha_total
delta <- alpha_total / 2
alpha_inf <- alpha_total / 2
z_baseline <- qnorm(1 - alpha_baseline / 2)
z_alpha_inf <- qnorm(1 - alpha_inf / 2)
reps <- 1000

designs <- load_designs()
settings <- tibble(
  augmentation = c("Dropout", "Dropout", "Dropout", "Gaussian", "Gaussian", "Gaussian"),
  strength = c(0.10, 0.25, 0.40, 0.05, 0.10, 0.20),
  stringsAsFactors = FALSE
)

rows <- list()
direction_rows <- list()
design_rows <- list()

# Iterate over the four designs (Gaussian and 3 UCI benchmarks)
for (design_name in names(designs)) {
  X_full <- designs[[design_name]]
  n <- nrow(X_full)
  p <- ncol(X_full)

  # Create calibration and inference data
  calibration_rows <- sort(sample.int(n, floor(n / 2), replace = FALSE))
  X_cal <- X_full[calibration_rows, , drop = FALSE]
  X_inf <- X_full[-calibration_rows, , drop = FALSE]

  theta <- sample(c(-1, 1), p, replace = TRUE) # simulate theta
  theta <- theta * sqrt(p) / sqrt(sum(theta * theta))
  sigma <- sd(as.vector(X_full %*% theta)) / 2 # calculate sigma
  directions <- target_directions(X_inf, theta) # get target directions
  M_cal <- crossprod(X_cal) # calculate matrix M of localization set
  M_cal_inv <- solve(M_cal)
  ellipsoid_radius <- sigma^2 * qchisq(1 - delta, df = p)

  design_rows[[length(design_rows) + 1]] <- data.frame(
    design = design_name,
    n = n,
    p = p,
    n_cal = nrow(X_cal),
    n_inf = nrow(X_inf),
    sigma = sigma,
    stringsAsFactors = FALSE
  )

  # Simulate calibration and inference noise
  errors_cal <- matrix(rnorm(nrow(X_cal) * reps, sd = sigma), nrow = nrow(X_cal))
  errors_inf <- matrix(rnorm(nrow(X_inf) * reps, sd = sigma), nrow = nrow(X_inf))

  # Calculate calibration and inference response (1000 samples)
  y_cal <- sweep(errors_cal, 1, as.vector(X_cal %*% theta), "+")
  y_inf <- sweep(errors_inf, 1, as.vector(X_inf %*% theta), "+")
  # Get 1000 centers as OLS estimators from calibration data
  centers <- solve(M_cal, crossprod(X_cal, y_cal)) #

  # Iterate over the six augmentation settings
  for (setting_id in seq_len(nrow(settings))) {
    kind <- settings$augmentation[setting_id]
    strength <- settings$strength[setting_id]

    # Get covariance matrix for setting and strength
    lambda_aug <- augmentation_penalty(X_inf, kind, strength)
    # Calculate matrices, standard deviation and bias
    temps <- lapply(directions, function(v) prepare_direction(X_inf, theta, sigma, v, lambda_aug))

    # Allocate space for interval coverage and widths
    covered <- list(Credible = numeric(), Oracle = numeric(), Split = numeric())
    widths <- list(Credible = numeric(), Oracle = numeric(), Split = numeric())

    # Iterate over the 3 target directions (signal, leading, trailing)
    for (direction_name in names(temps)) {
      temp <- temps[[direction_name]]
      gamma_hat <- as.vector(crossprod(temp$v, temp$A_inv_Xt %*% y_inf))

      # Calculate credible and oracle intervals
      credible_half <- z_baseline * temp$posterior_sd
      oracle_half <- temp$exact_bias + z_baseline * temp$sampling_sd

      # Calculate split interval
      center_bias <- as.vector(crossprod(temp$g, centers))
      support_term <- sqrt(ellipsoid_radius * as.numeric(crossprod(temp$g, M_cal_inv %*% temp$g)))
      split_half <- abs(center_bias) + support_term + z_alpha_inf * temp$sampling_sd

      stats <- list(
        Credible = list(half = rep(credible_half, reps)),
        Oracle = list(half = rep(oracle_half, reps)),
        Split = list(half = split_half)
      )
      oracle_width <- 2 * oracle_half # calculate oracle half as reference

      # Iterate over 3 methods (Credible, Oracle, Split)
      for (method in names(stats)) {
        half <- stats[[method]]$half

        # Check if gamma in interval (for 1000 replicates of gamma_hat)
        direction_covered <- as.numeric(abs(gamma_hat - temp$gamma) <= half)
        # Calculate interval width relative to oracle
        direction_width <- 2 * half / oracle_width

        # Add results to result table
        covered[[method]] <- c(covered[[method]], direction_covered)
        widths[[method]] <- c(widths[[method]], direction_width)
        direction_rows[[length(direction_rows) + 1]] <- data.frame(
          design = design_name,
          augmentation = kind,
          strength = strength,
          direction = direction_name,
          method = method,
          coverage = mean(direction_covered),
          coverage_se = sd(direction_covered) / sqrt(reps),
          rel_width = mean(direction_width),
          rel_width_se = sd(direction_width) / sqrt(reps),
          reps = reps,
          stringsAsFactors = FALSE
        )
      }
    }
    # Average the results over the three target directions
    for (method in names(covered)) {
      covered_by_rep <- rowMeans(matrix(covered[[method]], nrow = reps, ncol = length(directions)))
      widths_by_rep <- rowMeans(matrix(widths[[method]], nrow = reps, ncol = length(directions)))
      rows[[length(rows) + 1]] <- data.frame(
        design = design_name,
        n = n,
        p = p,
        n_cal = nrow(X_cal),
        n_inf = nrow(X_inf),
        augmentation = kind,
        strength = strength,
        method = method,
        coverage = mean(covered_by_rep),
        coverage_se = sd(covered_by_rep) / sqrt(reps),
        rel_width = mean(widths_by_rep),
        rel_width_se = sd(widths_by_rep) / sqrt(reps),
        reps = reps,
        alpha_total = alpha_total,
        alpha_baseline = alpha_baseline,
        alpha_inf = if (method == "Split") alpha_inf else alpha_baseline,
        eta = 0,
        delta = if (method == "Split") delta else 0,
        stringsAsFactors = FALSE
      )
    }
  }
  message(sprintf("finished %s", design_name))
}

# Create results table (averaged over all three directions)
results <- do.call(rbind, rows)
# Create results table with results for each direction
direction_results <- do.call(rbind, direction_rows)

# Summary table (Appendix D.1.3)
design_summary <- do.call(rbind, design_rows)
output_dir <- file.path("results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(output_dir, "classical_results.csv"), row.names = FALSE)
write.csv(direction_results, file.path(output_dir, "classical_results_by_direction.csv"), row.names = FALSE)
write.csv(design_summary, file.path(output_dir, "classical_design_summary.csv"), row.names = FALSE)

# Create statistics for Table 1 in Section 6.1
for (method in unique(results$method)) {
  sub <- results[results$method == method, ]
  message(sprintf(
    "%s min_cov=%.3f med_cov=%.3f med_width=%.3f",
    method,
    min(sub$coverage),
    median(sub$coverage),
    median(sub$rel_width)
  ))
}
