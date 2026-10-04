# High-dimensional sparse experiments for the results section.
#
# The experiment follows the notation used in the paper:
#   n_cal: sample size used to fit the split Lasso center.
#   n_inf: sample size used to form the data-augmentation estimator.
#   p:     ambient dimension, with p > n_inf in this benchmark.
#   s:     sparsity, i.e. the number of nonzero coefficients in theta.
#   alpha: target miscoverage level for the final confidence interval.
#   eta:   auxiliary probability level used in the certified bias bound.
#   delta: auxiliary probability level used in the Lasso radius.
#
#
set.seed(1000)


# Center each feature and rescale it to have Euclidean norm sqrt(n)
standardize_columns <- function(X) {
  X <- sweep(X, 2, colMeans(X), "-")
  norms <- sqrt(colSums(X * X))
  norms[norms == 0] <- 1
  sweep(X, 2, norms, "/") * sqrt(nrow(X))
}


# Generate the fixed Gaussian design matrix
toeplitz_design <- function(n, p, rho) {
  z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  if (rho == 0) {
    return(standardize_columns(z))
  }
  X <- matrix(0, nrow = n, ncol = p)
  X[, 1] <- z[, 1]
  scale <- sqrt(1 - rho^2)
  for (j in 2:p) {
    X[, j] <- rho * X[, j - 1] + scale * z[, j]
  }
  standardize_columns(X)
}

# Construct a fixed high-dimensional design with independent signal and nuisance blocks
compatible_toeplitz_design <- function(n, p, sparsity, rho) {
  # Signal block and nuisance block with Toeplitz correlation
  X_signal <- toeplitz_design(n, sparsity, rho)
  X_nuisance <- toeplitz_design(n, p - sparsity, rho)

  signal_basis <- qr.Q(qr(cbind(rep(1, n), X_signal)))

  # project nuisance block onto the orthogonal complement of span(1, X_signal)
  X_nuisance <- X_nuisance - signal_basis %*% crossprod(signal_basis, X_nuisance)
  # normalize columns
  X_nuisance <- standardize_columns(X_nuisance)
  X <- cbind(X_signal, X_nuisance)

  cross_block <- max(abs(crossprod(X_signal, X_nuisance))) / n
  if (cross_block > 1e-10) {
    stop(sprintf("Block orthogonality check failed: %.3e", cross_block))
  }

  # calculate certified compatibility constant
  phi0_sq <- min(eigen(crossprod(X_signal) / n, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(phi0_sq) || phi0_sq <= 0) {
    stop("The signal block does not have a positive compatibility constant")
  }

  list(X = X, phi0_sq = phi0_sq)
}


# Coordinate-wise soft thresholding operator for the Lasso update
soft_threshold <- function(X, lambda) {
  sign(X) * pmax(abs(X) - lambda, 0)
}

# Fit the Lasso center used by the split interval
lasso_cd <- function(X, y, lambda, max_iter = 5000, tol = 1e-10) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  residual <- y
  col_norm <- colSums(X * X) / n
  col_norm[col_norm == 0] <- 1

  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    old <- beta
    for (j in seq_len(p)) {
      residual <- residual + X[, j] * beta[j]
      raw <- sum(X[, j] * residual) / n
      beta[j] <- soft_threshold(raw, lambda) / col_norm[j]
      residual <- residual - X[, j] * beta[j]
    }
    if (max(abs(beta - old)) <= tol * max(1, max(abs(old)))) {
      converged <- TRUE
      break
    }
  }
  gradient <- -as.vector(crossprod(X, residual)) / n
  active <- abs(beta) > 1e-9
  kkt_error <- 0
  if (any(active)) {
    kkt_error <- max(kkt_error, max(abs(gradient[active] + lambda * sign(beta[active]))))
  }
  if (any(!active)) {
    kkt_error <- max(kkt_error, max(pmax(abs(gradient[!active]) - lambda, 0)))
  }
  if (!converged || kkt_error > 1e-7) {
    stop(sprintf("Lasso coordinate descent failed to converge (KKT error %.3e)", kkt_error))
  }
  beta
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

# Choose a certificate vector a for the certified split-Lasso bias bound
# The chosen a minimizes, over the deterministic ridge grid and a = 0, the
# upper bound
choose_certificate <- function(X, g, radius, sigma, z_eta) {
  XX_t <- tcrossprod(X)
  best_a <- rep(0, nrow(X))
  best_value <- radius * max(abs(g))
  scale <- sum(diag(XX_t)) / nrow(XX_t)

  for (mult in 10^seq(-5, 3, length.out = 25)) {
    kappa <- mult * scale
    a <- solve(XX_t + kappa * diag(nrow(XX_t)), X %*% g)
    value <- z_eta * sigma * sqrt(sum(a * a)) + radius * max(abs(g - as.vector(crossprod(X, a))))
    if (value < best_value) {
      best_value <- value
      best_a <- as.vector(a)
    }
  }
  best_a
}

# Get prediction directions
target_directions <- function(X, theta) {
  sv <- svd(X, nu = 0, nv = min(dim(X)))
  rank_tol <- max(dim(X)) * .Machine$double.eps * max(sv$d)
  nonzero_rank <- sum(sv$d > rank_tol)
  list(
    signal = theta / sqrt(sum(theta * theta)), # signal direction
    leading = sv$v[, 1], # leading right singular vector
    trailing_nonzero = sv$v[, nonzero_rank] # trailing right singular vector
  )
}

# Precompute the quantities shared by all Monte Carlo repetitions. The matrix
# inverse is common to every target direction, so compute it only once.
prepare_directions <- function(X, theta, sigma, directions, lambda_aug) {
  p <- ncol(X)
  XtX <- crossprod(X)
  a_inv <- solve(XtX + lambda_aug)
  s_mat <- a_inv %*% XtX
  bias_operator <- s_mat - diag(p)
  a_inv_Xt <- a_inv %*% t(X)

  lapply(directions, function(v) {
    g <- as.vector(crossprod(bias_operator, v))
    list(
      v = v,
      gamma = as.numeric(crossprod(v, theta)),
      g = g,
      sampling_sd = sigma * sqrt(as.numeric(t(v) %*% a_inv %*% XtX %*% a_inv %*% v)),
      posterior_sd = sigma * sqrt(as.numeric(t(v) %*% a_inv %*% v)),
      exact_bias = abs(as.numeric(crossprod(g, theta))),
      a_inv_Xt = a_inv_Xt
    )
  })
}

interval_stats <- function(
    prep,
    X,
    y,
    sigma,
    z_baseline,
    z_split,
    z_split_zero,
    split_radius,
    center,
    certificate,
    z_eta) {
  theta_hat <- prep$a_inv_Xt %*% y
  gamma_hat <- as.numeric(crossprod(prep$v, theta_hat))

  # Credible and oracle interval
  credible_half <- z_baseline * prep$posterior_sd
  oracle_half <- prep$exact_bias + z_baseline * prep$sampling_sd

  # Zero-certicate split interval
  center_bias <- as.numeric(crossprod(prep$g, center))
  bias_bound <- abs(center_bias) + split_radius * max(abs(prep$g))
  split_zero_half <- bias_bound + z_split_zero * prep$sampling_sd

  # Certified split-Lasso interval
  residual_direction <- prep$g - as.vector(crossprod(X, certificate))
  certified_bias <- abs(center_bias + sum(certificate * (y - X %*% center))) +
    z_eta * sigma * sqrt(sum(certificate * certificate)) +
    split_radius * max(abs(residual_direction))
  split_cert_half <- certified_bias + z_split * prep$sampling_sd

  list(
    Credible = c(covered = as.numeric(abs(gamma_hat - prep$gamma) <= credible_half), width = 2 * credible_half),
    Oracle = c(covered = as.numeric(abs(gamma_hat - prep$gamma) <= oracle_half), width = 2 * oracle_half),
    "Split-Lasso-0" = c(covered = as.numeric(abs(gamma_hat - prep$gamma) <= split_zero_half), width = 2 * split_zero_half),
    "Split-Lasso-cert" = c(covered = as.numeric(abs(gamma_hat - prep$gamma) <= split_cert_half), width = 2 * split_cert_half)
  )
}


alpha_total <- 0.025 # total failure probability
alpha_baseline <- alpha_total
alpha_inf <- alpha_total / 3 # sampling noise
eta <- alpha_total / 3 # certificate noise
delta <- alpha_total / 3 # calibration failure
alpha_inf_zero <- alpha_total - delta
z_baseline <- qnorm(1 - alpha_baseline / 2)
z_alpha_inf <- qnorm(1 - alpha_inf / 2)
z_alpha_inf_zero <- qnorm(1 - alpha_inf_zero / 2)
z_eta <- qnorm(1 - eta / 2)

reps <- 500

# Create dataframe for all four designs
designs <- data.frame(
  design = c(
    "Block iid", "Block corr", "Block high-corr", "Block dense-support"
  ),
  source = rep("Simulated", 4),
  n_cal = rep(120, 4),
  n_inf = rep(120, 4),
  p = rep(300, 4),
  rho = c(0, 0.5, 0.8, 0.5),
  sparsity = c(10, 10, 10, 20),
  stringsAsFactors = FALSE
)

# Create dataframe for all combinations of augmentation settings
settings <- data.frame(
  augmentation = c("Dropout", "Dropout", "Dropout", "Gaussian", "Gaussian", "Gaussian"),
  strength = c(0.10, 0.25, 0.40, 0.05, 0.10, 0.20),
  stringsAsFactors = FALSE
)


methods <- c("Credible", "Oracle", "Split-Lasso-0", "Split-Lasso-cert")
rows <- list()
direction_rows <- list()
design_rows <- list()

# Iiterate over desings
for (di in seq_len(nrow(designs))) {
  design_name <- designs$design[di]
  design_source <- designs$source[di]
  n_cal <- designs$n_cal[di]
  n_inf <- designs$n_inf[di]
  p <- designs$p[di]
  rho <- designs$rho[di]
  sparsity <- designs$sparsity[di]

  cal_design <- compatible_toeplitz_design(n_cal, p, sparsity, rho)
  inf_design <- compatible_toeplitz_design(n_inf, p, sparsity, rho)
  X_cal <- cal_design$X
  X_inf <- inf_design$X
  phi0_sq <- cal_design$phi0_sq
  support_features <- paste0("coordinates 1--", sparsity)


  # Sparse regression vector theta.
  # Only the first s entries are nonzero, with random signs
  support <- seq_len(sparsity)
  theta <- rep(0, p)
  theta[support] <- sample(c(-1, 1), sparsity, replace = TRUE)
  theta <- theta * sqrt(sparsity) / sqrt(sum(theta * theta))

  # Calculate sigma (As described in section 6.1. and 6.2)
  signal <- c(X_cal %*% theta, X_inf %*% theta)
  sigma <- sd(signal) / 2

  # Get prediction directions
  directions <- target_directions(X_inf, theta)

  # Lasso tuning and l1 radius for the split-Lasso uncertainty set
  lambda_lasso <- 2 * sigma * sqrt(2 * log(2 * p / delta) / n_cal)
  radius <- 12 * sparsity * lambda_lasso / phi0_sq

  design_rows[[length(design_rows) + 1L]] <- data.frame(
    design = design_name,
    source = design_source,
    n_cal = n_cal,
    n_inf = n_inf,
    p = p,
    sparsity = sparsity,
    rho = rho,
    sigma = sigma,
    phi0_sq = phi0_sq,
    lambda_lasso = lambda_lasso,
    radius = radius,
    support_features = support_features,
    stringsAsFactors = FALSE
  )

  # Iterate over augmentation settings
  for (si in seq_len(nrow(settings))) {
    aug_kind <- settings$augmentation[si]
    strength <- settings$strength[si]
    lambda_aug <- augmentation_penalty(X_inf, aug_kind, strength)

    # Compute direction-specific deterministic quantities and certificates
    preps <- prepare_directions(X_inf, theta, sigma, directions, lambda_aug)
    certs <- list()
    for (direction_name in names(directions)) {
      certs[[direction_name]] <- choose_certificate(X_inf, preps[[direction_name]]$g, radius, sigma, z_eta)
    }


    covered <- setNames(vector("list", length(methods)), methods)
    widths <- setNames(vector("list", length(methods)), methods)
    for (m in methods) {
      covered[[m]] <- matrix(0, nrow = reps, ncol = length(directions), dimnames = list(NULL, names(directions)))
      widths[[m]] <- matrix(0, nrow = reps, ncol = length(directions), dimnames = list(NULL, names(directions)))
    }

    # 500 reps for each augmentation setting - design combination
    for (rep_id in seq_len(reps)) {
      # Generate calibration and inference responses
      y_cal <- as.vector(X_cal %*% theta + rnorm(n_cal, sd = sigma))
      y_inf <- as.vector(X_inf %*% theta + rnorm(n_inf, sd = sigma))

      # Fit the split Lasso center on the calibration data
      center <- lasso_cd(X_cal, y_cal, lambda_lasso)

      # Iterate over directions
      for (direction_idx in seq_along(directions)) {
        direction_name <- names(directions)[direction_idx]

        stats <- interval_stats(
          preps[[direction_name]],
          X_inf,
          y_inf,
          sigma,
          z_baseline,
          z_alpha_inf,
          z_alpha_inf_zero,
          radius,
          center,
          certs[[direction_name]],
          z_eta
        )

        # store coverage and calculate width with respect to oracle width
        oracle_width <- stats[["Oracle"]]["width"]
        for (method in names(stats)) {
          covered[[method]][rep_id, direction_idx] <- stats[[method]]["covered"]
          widths[[method]][rep_id, direction_idx] <- stats[[method]]["width"] / oracle_width
        }
      }
    }

    # Average the results over the three target directions
    for (method in methods) {
      covered_by_rep <- rowMeans(covered[[method]])
      widths_by_rep <- rowMeans(widths[[method]])
      rows[[length(rows) + 1]] <- data.frame(
        design = design_name,
        source = design_source,
        rho = rho,
        sparsity = sparsity,
        n_cal = n_cal,
        n_inf = n_inf,
        p = p,
        phi0_sq = phi0_sq,
        radius = radius,
        alpha_total = alpha_total,
        alpha_baseline = alpha_baseline,
        alpha_inf = if (method == "Split-Lasso-0") alpha_inf_zero else if (method == "Split-Lasso-cert") alpha_inf else alpha_baseline,
        eta = if (method == "Split-Lasso-cert") eta else 0,
        delta = if (grepl("Split-Lasso", method, fixed = TRUE)) delta else 0,
        augmentation = aug_kind,
        strength = strength,
        method = method,
        coverage = mean(covered_by_rep),
        coverage_se = sd(covered_by_rep) / sqrt(reps),
        rel_width = mean(widths_by_rep),
        rel_width_se = sd(widths_by_rep) / sqrt(reps),
        reps = reps,
        stringsAsFactors = FALSE
      )

      # Create results table (keep results for each prediction direction)
      for (direction_name in names(directions)) {
        direction_rows[[length(direction_rows) + 1]] <- data.frame(
          design = design_name,
          source = design_source,
          rho = rho,
          sparsity = sparsity,
          n_cal = n_cal,
          n_inf = n_inf,
          p = p,
          phi0_sq = phi0_sq,
          radius = radius,
          alpha_total = alpha_total,
          alpha_baseline = alpha_baseline,
          alpha_inf = if (method == "Split-Lasso-0") alpha_inf_zero else if (method == "Split-Lasso-cert") alpha_inf else alpha_baseline,
          eta = if (method == "Split-Lasso-cert") eta else 0,
          delta = if (grepl("Split-Lasso", method, fixed = TRUE)) delta else 0,
          augmentation = aug_kind,
          strength = strength,
          direction = direction_name,
          method = method,
          coverage = mean(covered[[method]][, direction_name]),
          coverage_se = sd(covered[[method]][, direction_name]) / sqrt(reps),
          rel_width = mean(widths[[method]][, direction_name]),
          rel_width_se = sd(widths[[method]][, direction_name]) / sqrt(reps),
          reps = reps,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  message(sprintf("finished %s", design_name))
}

# Create results table (averaged over all three directions)
results <- do.call(rbind, rows)

# Create results table with results for each direction
direction_results <- do.call(rbind, direction_rows)

# summary of designs
design_summary <- do.call(rbind, design_rows)
output_dir <- file.path("results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(output_dir, "high_dim_results.csv"), row.names = FALSE)
write.csv(direction_results, file.path(output_dir, "high_dim_results_by_direction.csv"), row.names = FALSE)
write.csv(design_summary, file.path(output_dir, "high_dim_design_summary.csv"), row.names = FALSE)


# Create statistics for Table 2 in Section 6.2
for (method in methods) {
  sub <- results[results$method == method, ]
  message(sprintf(
    "%s min_cov=%.3f med_cov=%.3f med_width=%.3f",
    method,
    min(sub$coverage),
    median(sub$coverage),
    median(sub$rel_width)
  ))
}
