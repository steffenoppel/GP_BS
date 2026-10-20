################################################################################
# Northern Giant Petrel winter attendance at Archway, Marion Island
# FINAL ANALYSIS SCRIPT
#
# Purpose
# -------
# Reproduce the final detection-corrected winter-attendance analysis used for
# the manuscript, including:
#   1. preparation of bird-by-check encounter histories (2022-2023);
#   2. the marginal two-state dynamic hidden Markov model (HMM);
#   3. MCMC diagnostics and posterior summaries;
#   4. forward-backward smoothing to estimate detection-corrected winter
#      attendance for each bird-winter;
#   5. subsequent-breeding and breeding-success GLMMs;
#   6. propagation of HMM attendance uncertainty through 250 downstream fits;
#   7. final model summaries and manuscript figure code.
#
# IMPORTANT
# ---------
# The HMM is restricted to winters 2022 and 2023 because comparable search-
# effort metadata were not available for 2021. The analysis starts from the
# final analysis-ready file:
#   ngp_archway_detection_model_input_field_weather.csv
#
# The observation categories used by nimbleEcology are:
#   1 = not detected
#   2 = detected
# The latent biological states are:
#   1 = absent from Archway
#   2 = attending Archway
#
# This is a two-stage analysis with propagated posterior attendance uncertainty,
# not a single fully joint likelihood for attendance and reproductive outcomes.
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(1001)

# ==============================================================================
# 0. PACKAGES AND PATHS
# ==============================================================================

required_packages <- c(
  "nimble",
  "nimbleEcology",
  "coda",
  "glmmTMB"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE,
          FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(nimble)
library(nimbleEcology)
library(coda)
library(glmmTMB)

input_file <- "ngp_archway_detection_model_input_field_weather.csv"

if (!file.exists(input_file)) {
  stop(
    "Cannot find '", input_file, "' in the working directory. ",
    "Place the final analysis-ready CSV beside this script or edit input_file."
  )
}

output_dir <- "NGP_final_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir)

# ==============================================================================
# 1. LOAD AND CHECK THE ANALYSIS-READY DATA
# ==============================================================================

dat_det <- read.csv(input_file, stringsAsFactors = FALSE)
dat_det$check_date <- as.Date(dat_det$check_date)

required_columns <- c(
  "id", "winter", "check_index", "check_date", "day_of_year", "detected",
  "prev_success", "experience", "sex", "future_bred", "future_success",
  "search_minutes", "n_observers"
)

stopifnot(all(required_columns %in% names(dat_det)))
stopifnot(all(dat_det$detected %in% c(0, 1)))

cat("\nFull encounter-history dataset:\n")
cat("Rows:", nrow(dat_det), "\n")
cat("Bird-winter records by winter:\n")
print(tapply(dat_det$id, dat_det$winter, function(x) length(unique(x))))

# Primary detection-corrected analysis: 2022 and 2023 only.
dat_det_primary <- subset(dat_det, winter %in% c(2022, 2023))

if (anyNA(dat_det_primary$search_minutes) ||
    anyNA(dat_det_primary$n_observers)) {
  stop("Search effort or observer number is missing in the 2022-2023 primary data.")
}

# ==============================================================================
# 2. BIRD-WINTER COVARIATE TABLE
# ==============================================================================

bird_fields <- c(
  "id", "winter", "prev_success", "experience", "sex",
  "future_bred", "future_success"
)

bird_cov <- unique(dat_det_primary[, bird_fields])
bird_cov <- bird_cov[order(bird_cov$winter, as.character(bird_cov$id)), ]
row.names(bird_cov) <- NULL

# Verify one covariate record per individual per winter.
stopifnot(!anyDuplicated(bird_cov[, c("id", "winter")]))
stopifnot(nrow(bird_cov) == 190)

# Standardise breeding experience across the 190 bird-winter histories.
bird_cov$experience_c <-
  bird_cov$experience - mean(bird_cov$experience, na.rm = TRUE)
bird_cov$experience_z <-
  bird_cov$experience_c / sd(bird_cov$experience_c, na.rm = TRUE)

bird_cov$sex_male <- as.integer(bird_cov$sex == "M")

# Shared individual index across winters.
id_factor <- factor(bird_cov$id)
bird_cov$id_index <- as.integer(id_factor)
J <- nlevels(id_factor)

idx22 <- which(bird_cov$winter == 2022)
idx23 <- which(bird_cov$winter == 2023)

stopifnot(length(idx22) == 97)
stopifnot(length(idx23) == 93)
stopifnot(J == 130)

id22 <- bird_cov$id_index[idx22]
id23 <- bird_cov$id_index[idx23]

cat("\nPrimary HMM sample:\n")
cat("Bird-winter histories:", nrow(bird_cov), "\n")
cat("Unique individuals:", J, "\n")
cat("2022 histories:", length(idx22), "\n")
cat("2023 histories:", length(idx23), "\n")

# ==============================================================================
# 3. BUILD DETECTION MATRICES AND CHECK-LEVEL COVARIATES
# ==============================================================================

make_detection_matrix <- function(year, dat, bird_cov, idx) {
  d <- dat[dat$winter == year, ]
  d <- d[order(d$id, d$check_index), ]

  ids <- as.character(bird_cov$id[idx])
  checks <- sort(unique(d$check_index))

  out <- matrix(
    NA_integer_,
    nrow = length(ids),
    ncol = length(checks),
    dimnames = list(ids, paste0("check_", checks))
  )

  for (i in seq_along(ids)) {
    dd <- d[as.character(d$id) == ids[i], ]
    dd <- dd[match(checks, dd$check_index), ]
    if (nrow(dd) != length(checks) || anyNA(dd$detected)) {
      stop("Incomplete detection history for bird ", ids[i], " in ", year)
    }
    out[i, ] <- as.integer(dd$detected)
  }

  out
}

y22 <- make_detection_matrix(2022, dat_det_primary, bird_cov, idx22)
y23 <- make_detection_matrix(2023, dat_det_primary, bird_cov, idx23)

stopifnot(all(dim(y22) == c(97, 11)))
stopifnot(all(dim(y23) == c(93, 14)))
stopifnot(sum(y22) == 115)
stopifnot(sum(y23) == 169)

# nimbleEcology observation states must be positive integers:
# 1 = not detected, 2 = detected.
x22 <- y22 + 1L
x23 <- y23 + 1L

make_check_metadata <- function(year, dat) {
  d <- dat[dat$winter == year, c(
    "check_index", "check_date", "day_of_year", "search_minutes", "n_observers"
  )]
  d <- unique(d)
  d <- d[order(d$check_index), ]
  row.names(d) <- NULL
  d
}

check22 <- make_check_metadata(2022, dat_det_primary)
check23 <- make_check_metadata(2023, dat_det_primary)

# Standardise date across all 25 survey occasions, not separately by winter.
all_checks <- rbind(
  transform(check22, winter = 2022),
  transform(check23, winter = 2023)
)
all_checks$day_z <- as.numeric(scale(all_checks$day_of_year))

check22$day_z <- all_checks$day_z[all_checks$winter == 2022]
check23$day_z <- all_checks$day_z[all_checks$winter == 2023]

# Detection exposure is search duration in hours.
eff22 <- check22$search_minutes / 60
eff23 <- check23$search_minutes / 60

two22 <- as.integer(check22$n_observers == 2)
two23 <- as.integer(check23$n_observers == 2)

dayz22 <- check22$day_z
dayz23 <- check23$day_z

gap22 <- as.numeric(diff(check22$check_date))
gap23 <- as.numeric(diff(check23$check_date))

cat("\nDetection-history checks:\n")
cat("2022 detections:", sum(y22), "across", ncol(y22), "checks\n")
cat("2023 detections:", sum(y23), "across", ncol(y23), "checks\n")
cat("Gap range:", min(c(gap22, gap23)), "to", max(c(gap22, gap23)), "days\n")

# ==============================================================================
# 4. FINAL MARGINAL DYNAMIC HMM WITH SHARED INDIVIDUAL RANDOM INTERCEPT
# ==============================================================================
#
# Attendance propensity:
# logit(psi_it) = alpha_att
#               + beta_prev * previous breeding success
#               + beta_exp  * standardised breeding experience
#               + beta_sex  * male
#               + beta_year * winter 2023
#               + beta_date * standardised day of year
#               + u_id
#
# Detection model conditional on attending:
# p_t = 1 - exp[-search_hours_t * exp(det0 + det_obs * two_observers_t)]
#
# Gap-adjusted transition model:
# decay = exp(-switch_rate * gap_days)
# P(0->1) = psi_next * (1-decay)
# P(1->1) = psi_next + (1-psi_next) * decay
#
# dDHMMo analytically marginalises the latent state sequence.
# ==============================================================================

re_hmm_code <- nimbleCode({

  # Priors for attendance fixed effects
  alpha_att ~ dnorm(0, sd = 2)
  beta_prev ~ dnorm(0, sd = 1.5)
  beta_exp  ~ dnorm(0, sd = 1.5)
  beta_sex  ~ dnorm(0, sd = 1.5)
  beta_year ~ dnorm(0, sd = 1.5)
  beta_date ~ dnorm(0, sd = 1.5)

  # Persistent individual heterogeneity
  sigma_id ~ dexp(1)
  for (j in 1:J) {
    u_raw[j] ~ dnorm(0, sd = 1)
    u_id[j] <- sigma_id * u_raw[j]
  }

  # State switching rate
  log_switch ~ dnorm(-4, sd = 1.5)
  switch_rate <- exp(log_switch)

  # Detection parameters
  det0    ~ dnorm(-1, sd = 1.5)
  det_obs ~ dnorm(0, sd = 1)

  # --------------------------------------------------------------------------
  # 2022 observation model
  # --------------------------------------------------------------------------
  for (t in 1:T22) {
    det_rate22[t] <- exp(det0 + det_obs * two22[t])
    p22[t] <- 1 - exp(-eff22[t] * det_rate22[t])

    # latent state 1 = absent: detection impossible
    probObs22[1, 1, t] <- 1
    probObs22[1, 2, t] <- 0

    # latent state 2 = attending
    probObs22[2, 1, t] <- 1 - p22[t]
    probObs22[2, 2, t] <- p22[t]
  }

  for (t in 1:(T22 - 1)) {
    decay22[t] <- exp(-switch_rate * gap22[t])
  }

  for (i in 1:N22) {
    for (t in 1:T22) {
      eta22[i, t] <- alpha_att +
                     beta_prev * prev22[i] +
                     beta_exp  * expz22[i] +
                     beta_sex  * male22[i] +
                     beta_date * dayz22[t] +
                     u_id[id22[i]]
      psi22[i, t] <- ilogit(eta22[i, t])
    }

    init22[i, 1] <- 1 - psi22[i, 1]
    init22[i, 2] <- psi22[i, 1]

    for (t in 1:(T22 - 1)) {
      # rows = state at t; columns = state at t+1
      probTrans22[1, 1, t, i] <-
        1 - psi22[i, t + 1] * (1 - decay22[t])
      probTrans22[1, 2, t, i] <-
        psi22[i, t + 1] * (1 - decay22[t])
      probTrans22[2, 1, t, i] <-
        (1 - psi22[i, t + 1]) * (1 - decay22[t])
      probTrans22[2, 2, t, i] <-
        psi22[i, t + 1] + (1 - psi22[i, t + 1]) * decay22[t]
    }

    x22[i, 1:T22] ~ dDHMMo(
      init22[i, 1:2],
      probObs = probObs22[1:2, 1:2, 1:T22],
      probTrans = probTrans22[1:2, 1:2, 1:(T22 - 1), i],
      len = T22,
      checkRowSums = 0
    )
  }

  # --------------------------------------------------------------------------
  # 2023 observation model
  # --------------------------------------------------------------------------
  for (t in 1:T23) {
    det_rate23[t] <- exp(det0 + det_obs * two23[t])
    p23[t] <- 1 - exp(-eff23[t] * det_rate23[t])

    probObs23[1, 1, t] <- 1
    probObs23[1, 2, t] <- 0
    probObs23[2, 1, t] <- 1 - p23[t]
    probObs23[2, 2, t] <- p23[t]
  }

  for (t in 1:(T23 - 1)) {
    decay23[t] <- exp(-switch_rate * gap23[t])
  }

  for (i in 1:N23) {
    for (t in 1:T23) {
      eta23[i, t] <- alpha_att +
                     beta_prev * prev23[i] +
                     beta_exp  * expz23[i] +
                     beta_sex  * male23[i] +
                     beta_year +
                     beta_date * dayz23[t] +
                     u_id[id23[i]]
      psi23[i, t] <- ilogit(eta23[i, t])
    }

    init23[i, 1] <- 1 - psi23[i, 1]
    init23[i, 2] <- psi23[i, 1]

    for (t in 1:(T23 - 1)) {
      probTrans23[1, 1, t, i] <-
        1 - psi23[i, t + 1] * (1 - decay23[t])
      probTrans23[1, 2, t, i] <-
        psi23[i, t + 1] * (1 - decay23[t])
      probTrans23[2, 1, t, i] <-
        (1 - psi23[i, t + 1]) * (1 - decay23[t])
      probTrans23[2, 2, t, i] <-
        psi23[i, t + 1] + (1 - psi23[i, t + 1]) * decay23[t]
    }

    x23[i, 1:T23] ~ dDHMMo(
      init23[i, 1:2],
      probObs = probObs23[1:2, 1:2, 1:T23],
      probTrans = probTrans23[1:2, 1:2, 1:(T23 - 1), i],
      len = T23,
      checkRowSums = 0
    )
  }
})

re_constants <- list(
  J = J,
  N22 = nrow(x22),
  N23 = nrow(x23),
  T22 = ncol(x22),
  T23 = ncol(x23),

  prev22 = as.numeric(bird_cov$prev_success[idx22]),
  prev23 = as.numeric(bird_cov$prev_success[idx23]),

  expz22 = as.numeric(bird_cov$experience_z[idx22]),
  expz23 = as.numeric(bird_cov$experience_z[idx23]),

  male22 = as.numeric(bird_cov$sex_male[idx22]),
  male23 = as.numeric(bird_cov$sex_male[idx23]),

  id22 = as.integer(id22),
  id23 = as.integer(id23),

  eff22 = as.numeric(eff22),
  eff23 = as.numeric(eff23),
  two22 = as.numeric(two22),
  two23 = as.numeric(two23),
  dayz22 = as.numeric(dayz22),
  dayz23 = as.numeric(dayz23),
  gap22 = as.numeric(gap22),
  gap23 = as.numeric(gap23)
)

re_data <- list(
  x22 = x22,
  x23 = x23
)

re_inits <- list(
  alpha_att = 0,
  beta_prev = 0,
  beta_exp = 0,
  beta_sex = 0,
  beta_year = 0,
  beta_date = 0,
  sigma_id = 0.5,
  u_raw = rep(0, J),
  log_switch = -4,
  det0 = -1,
  det_obs = 0
)

# Parameters monitored for the final fit. u_id is required for the subsequent
# forward-backward smoothing step.
re_monitors <- c(
  "alpha_att", "beta_prev", "beta_exp", "beta_sex", "beta_year",
  "beta_date", "sigma_id", "switch_rate", "det0", "det_obs", "u_id"
)

cat("\nBuilding final HMM...\n")
R_re_hmm <- nimbleModel(
  code = re_hmm_code,
  constants = re_constants,
  data = re_data,
  inits = re_inits
)

cat("Initial model log-likelihood:", R_re_hmm$calculate(), "\n")

C_re_hmm <- compileNimble(R_re_hmm)

re_conf <- configureMCMC(
  R_re_hmm,
  monitors = re_monitors,
  enableWAIC = FALSE,
  print = FALSE
)

R_re_mcmc <- buildMCMC(re_conf)
C_re_mcmc <- compileNimble(R_re_mcmc, project = R_re_hmm)

# Final MCMC settings used for inference.
set.seed(1001)
cat("\nRunning 3 HMM chains (30,000 iterations; 5,000 burn-in each)...\n")

re_samples_att <- runMCMC(
  C_re_mcmc,
  niter = 30000,
  nburnin = 5000,
  nchains = 3,
  inits = list(re_inits, re_inits, re_inits),
  samplesAsCodaMCMC = TRUE,
  progressBar = TRUE
)

# ==============================================================================
# 5. HMM DIAGNOSTICS AND POSTERIOR SUMMARIES
# ==============================================================================

key_pars <- c(
  "alpha_att", "beta_prev", "beta_exp", "beta_sex", "beta_year",
  "beta_date", "sigma_id", "switch_rate", "det0", "det_obs"
)

as_mcmc_subset <- function(samples, columns) {
  mcmc.list(lapply(samples, function(ch) {
    as.mcmc(as.matrix(ch)[, columns, drop = FALSE])
  }))
}

key_samples <- as_mcmc_subset(re_samples_att, key_pars)

re_gelman_final <- gelman.diag(key_samples, autoburnin = FALSE)
re_ess_final <- effectiveSize(key_samples)
re_summary_final <- summary(key_samples)

cat("\n=== HMM GELMAN-RUBIN DIAGNOSTICS ===\n")
print(re_gelman_final)

cat("\n=== HMM EFFECTIVE SAMPLE SIZES ===\n")
print(re_ess_final)

cat("\n=== HMM POSTERIOR SUMMARY ===\n")
print(re_summary_final)

att_mat <- do.call(rbind, lapply(re_samples_att, as.matrix))

# Check convergence for all 130 monitored individual effects.
u_cols <- paste0("u_id[", seq_len(J), "]")
u_samples <- as_mcmc_subset(re_samples_att, u_cols)
u_gelman <- gelman.diag(u_samples, autoburnin = FALSE)
u_ess <- effectiveSize(u_samples)

cat("\nMaximum u_id R-hat point estimate:",
    max(u_gelman$psrf[, 1], na.rm = TRUE), "\n")
cat("Maximum u_id R-hat upper CI:",
    max(u_gelman$psrf[, 2], na.rm = TRUE), "\n")
cat("Minimum u_id ESS:", min(u_ess, na.rm = TRUE), "\n")
cat("Median u_id ESS:", median(u_ess, na.rm = TRUE), "\n")

posterior_summary <- function(x) {
  c(
    mean = mean(x),
    median = median(x),
    lower95 = unname(quantile(x, 0.025)),
    upper95 = unname(quantile(x, 0.975))
  )
}

hmm_summary_table <- t(vapply(
  key_pars,
  function(p) posterior_summary(att_mat[, p]),
  FUN.VALUE = c(mean = 0, median = 0, lower95 = 0, upper95 = 0)
))
hmm_summary_table <- data.frame(
  parameter = row.names(hmm_summary_table),
  hmm_summary_table,
  row.names = NULL
)

write.csv(
  hmm_summary_table,
  file.path(output_dir, "HMM_posterior_summary.csv"),
  row.names = FALSE
)

# Main previous-success effect.
beta_prev_draws <- att_mat[, "beta_prev"]
prob_prev_negative <- mean(beta_prev_draws < 0)
prev_OR_draws <- exp(beta_prev_draws)
prev_OR_CrI <- c(
  median = median(prev_OR_draws),
  lower95 = unname(quantile(prev_OR_draws, 0.025)),
  upper95 = unname(quantile(prev_OR_draws, 0.975))
)

cat("\nPrevious breeding success -> latent winter attendance:\n")
cat("P(beta_prev < 0) =", prob_prev_negative, "\n")
print(prev_OR_CrI)

# Conditional detection probabilities for each survey at posterior median
# parameters, useful for describing imperfect detection.
med_det0 <- median(att_mat[, "det0"])
med_det_obs <- median(att_mat[, "det_obs"])

detection_by_check <- rbind(
  data.frame(
    winter = 2022,
    check_index = check22$check_index,
    check_date = check22$check_date,
    search_minutes = check22$search_minutes,
    n_observers = check22$n_observers,
    p_detect_given_attendance =
      1 - exp(-eff22 * exp(med_det0 + med_det_obs * two22))
  ),
  data.frame(
    winter = 2023,
    check_index = check23$check_index,
    check_date = check23$check_date,
    search_minutes = check23$search_minutes,
    n_observers = check23$n_observers,
    p_detect_given_attendance =
      1 - exp(-eff23 * exp(med_det0 + med_det_obs * two23))
  )
)

write.csv(
  detection_by_check,
  file.path(output_dir, "conditional_detection_by_check.csv"),
  row.names = FALSE
)

# Save lean HMM results; compiled NIMBLE objects are intentionally not saved.
save(
  re_hmm_code, re_constants, re_data, re_inits,
  bird_cov, idx22, idx23, id22, id23,
  check22, check23, x22, x23,
  re_samples_att, re_summary_final, re_ess_final, re_gelman_final,
  prob_prev_negative, prev_OR_CrI,
  file = file.path(output_dir, "NGP_random_effect_HMM_FINAL.RData")
)

# ==============================================================================
# 6. FORWARD-BACKWARD SMOOTHING: DETECTION-CORRECTED WINTER ATTENDANCE
# ==============================================================================

smooth_attendance <- function(x, psi, p, gap, switch_rate) {
  Tn <- length(x)
  stopifnot(length(psi) == Tn, length(p) == Tn)
  if (Tn > 1) stopifnot(length(gap) == Tn - 1)

  # Emission probability for each latent state at each survey.
  # Columns: 1 = absent, 2 = attending.
  emit <- matrix(NA_real_, nrow = Tn, ncol = 2)

  for (t in seq_len(Tn)) {
    if (x[t] == 1L) {
      emit[t, 1] <- 1
      emit[t, 2] <- 1 - p[t]
    } else if (x[t] == 2L) {
      emit[t, 1] <- 0
      emit[t, 2] <- p[t]
    } else {
      stop("Observation states must be coded 1 or 2.")
    }
  }

  alpha <- matrix(0, nrow = Tn, ncol = 2)
  beta  <- matrix(0, nrow = Tn, ncol = 2)

  # Initial state distribution.
  alpha[1, ] <- c(1 - psi[1], psi[1]) * emit[1, ]
  s <- sum(alpha[1, ])
  if (!is.finite(s) || s <= 0) stop("Invalid forward probability at t = 1")
  alpha[1, ] <- alpha[1, ] / s

  # Forward pass.
  if (Tn > 1) {
    for (t in 2:Tn) {
      decay <- exp(-switch_rate * gap[t - 1])
      target <- psi[t]

      P <- matrix(c(
        1 - target * (1 - decay),
        target * (1 - decay),
        (1 - target) * (1 - decay),
        target + (1 - target) * decay
      ), nrow = 2, byrow = TRUE)

      alpha[t, ] <- as.numeric(alpha[t - 1, ] %*% P) * emit[t, ]
      s <- sum(alpha[t, ])
      if (!is.finite(s) || s <= 0) stop("Invalid forward probability at t = ", t)
      alpha[t, ] <- alpha[t, ] / s
    }
  }

  # Backward pass.
  beta[Tn, ] <- c(1, 1)

  if (Tn > 1) {
    for (t in (Tn - 1):1) {
      decay <- exp(-switch_rate * gap[t])
      target <- psi[t + 1]

      P <- matrix(c(
        1 - target * (1 - decay),
        target * (1 - decay),
        (1 - target) * (1 - decay),
        target + (1 - target) * decay
      ), nrow = 2, byrow = TRUE)

      beta[t, 1] <- sum(P[1, ] * emit[t + 1, ] * beta[t + 1, ])
      beta[t, 2] <- sum(P[2, ] * emit[t + 1, ] * beta[t + 1, ])

      s <- sum(beta[t, ])
      if (is.finite(s) && s > 0) beta[t, ] <- beta[t, ] / s
    }
  }

  smoothed <- alpha * beta
  smoothed <- smoothed / rowSums(smoothed)

  # Return posterior probability of state 2 = attending.
  smoothed[, 2]
}

# Use 2,000 evenly spaced posterior draws for smoothed attendance estimates.
n_smooth <- min(2000L, nrow(att_mat))
smooth_draw_idx <- unique(round(seq(1, nrow(att_mat), length.out = n_smooth)))
n_smooth <- length(smooth_draw_idx)

A_draws <- matrix(
  NA_real_,
  nrow = n_smooth,
  ncol = nrow(bird_cov),
  dimnames = list(NULL, paste0(as.character(bird_cov$id), "_", bird_cov$winter))
)

cat("\nCalculating detection-corrected attendance over", n_smooth,
    "posterior draws...\n")

for (d in seq_along(smooth_draw_idx)) {
  draw <- att_mat[smooth_draw_idx[d], ]

  alpha_att_d <- draw["alpha_att"]
  beta_prev_d <- draw["beta_prev"]
  beta_exp_d  <- draw["beta_exp"]
  beta_sex_d  <- draw["beta_sex"]
  beta_year_d <- draw["beta_year"]
  beta_date_d <- draw["beta_date"]
  switch_d    <- draw["switch_rate"]
  det0_d      <- draw["det0"]
  det_obs_d   <- draw["det_obs"]

  u_d <- as.numeric(draw[u_cols])

  p22_d <- 1 - exp(-eff22 * exp(det0_d + det_obs_d * two22))
  p23_d <- 1 - exp(-eff23 * exp(det0_d + det_obs_d * two23))

  # 2022 bird-winter histories
  for (i in seq_along(idx22)) {
    k <- idx22[i]
    eta <- alpha_att_d +
           beta_prev_d * bird_cov$prev_success[k] +
           beta_exp_d  * bird_cov$experience_z[k] +
           beta_sex_d  * bird_cov$sex_male[k] +
           beta_date_d * dayz22 +
           u_d[id22[i]]

    psi <- plogis(eta)

    sm <- smooth_attendance(
      x = x22[i, ],
      psi = psi,
      p = p22_d,
      gap = gap22,
      switch_rate = switch_d
    )

    A_draws[d, k] <- mean(sm)
  }

  # 2023 bird-winter histories
  for (i in seq_along(idx23)) {
    k <- idx23[i]
    eta <- alpha_att_d +
           beta_prev_d * bird_cov$prev_success[k] +
           beta_exp_d  * bird_cov$experience_z[k] +
           beta_sex_d  * bird_cov$sex_male[k] +
           beta_year_d +
           beta_date_d * dayz23 +
           u_d[id23[i]]

    psi <- plogis(eta)

    sm <- smooth_attendance(
      x = x23[i, ],
      psi = psi,
      p = p23_d,
      gap = gap23,
      switch_rate = switch_d
    )

    A_draws[d, k] <- mean(sm)
  }

  if (d %% 100 == 0 || d == n_smooth) {
    cat("  completed", d, "of", n_smooth, "draws\n")
  }
}

stopifnot(!anyNA(A_draws))
stopifnot(all(A_draws >= 0 & A_draws <= 1))

attendance_dc_mean <- colMeans(A_draws)
attendance_dc_median <- apply(A_draws, 2, median)
attendance_dc_lwr <- apply(A_draws, 2, quantile, probs = 0.025)
attendance_dc_upr <- apply(A_draws, 2, quantile, probs = 0.975)

attendance_summary <- cbind(
  bird_cov,
  attendance_dc_mean = attendance_dc_mean,
  attendance_dc_median = attendance_dc_median,
  attendance_dc_lower95 = attendance_dc_lwr,
  attendance_dc_upper95 = attendance_dc_upr
)

cat("\nPosterior-mean detection-corrected attendance:\n")
print(summary(attendance_dc_mean))

write.csv(
  attendance_summary,
  file.path(output_dir, "bird_winter_detection_corrected_attendance.csv"),
  row.names = FALSE
)

save(
  A_draws, smooth_draw_idx, attendance_summary, smooth_attendance,
  file = file.path(output_dir, "NGP_detection_corrected_attendance.RData")
)

# ==============================================================================
# 7. SUBSEQUENT BREEDING MODEL
# ==============================================================================
# All 190 bird-winter histories are included.
# The attendance coefficient is scaled per +10 percentage points.
# ==============================================================================

breed_dat <- bird_cov
breed_dat$attendance_dc <- attendance_dc_mean
breed_dat$attendance10_dc <- breed_dat$attendance_dc * 10
breed_dat$id <- factor(breed_dat$id)
breed_dat$sex <- relevel(factor(breed_dat$sex), ref = "F")
breed_dat$winter <- relevel(factor(breed_dat$winter), ref = "2022")

m_future_breed_dc <- glmmTMB(
  future_bred ~ attendance10_dc + prev_success + experience_z + sex + winter +
    (1 | id),
  family = binomial(link = "logit"),
  data = breed_dat
)

cat("\n=== Diagnostic model using posterior-mean attendance: future breeding ===\n")
print(summary(m_future_breed_dc))

# Propagate attendance uncertainty using 250 representative A_draws.
n_imp <- min(250L, nrow(A_draws))
imp_idx <- unique(round(seq(1, nrow(A_draws), length.out = n_imp)))
n_imp <- length(imp_idx)

breed_imp_results <- data.frame(
  draw = imp_idx,
  beta = NA_real_,
  se = NA_real_,
  sd_id = NA_real_,
  pdHess = NA
)

fit_control <- glmmTMBControl(
  optCtrl = list(iter.max = 10000, eval.max = 10000)
)

cat("\nPropagating attendance uncertainty through future-breeding model...\n")

for (k in seq_along(imp_idx)) {
  d <- imp_idx[k]
  dd <- breed_dat
  dd$attendance10_dc <- A_draws[d, ] * 10

  fit <- try(
    glmmTMB(
      future_bred ~ attendance10_dc + prev_success + experience_z + sex + winter +
        (1 | id),
      family = binomial(link = "logit"),
      data = dd,
      control = fit_control
    ),
    silent = TRUE
  )

  if (!inherits(fit, "try-error")) {
    cf <- fixef(fit)$cond
    vc <- vcov(fit)$cond

    breed_imp_results$beta[k] <- unname(cf["attendance10_dc"])
    breed_imp_results$se[k] <-
      unname(sqrt(diag(vc))["attendance10_dc"])
    breed_imp_results$sd_id[k] <-
      as.numeric(attr(VarCorr(fit)$cond$id, "stddev"))[1]
    breed_imp_results$pdHess[k] <- isTRUE(fit$sdr$pdHess)
  }
}

pool_rubin <- function(beta, se) {
  ok <- is.finite(beta) & is.finite(se)
  Q <- beta[ok]
  U <- se[ok]^2
  m <- length(Q)

  if (m < 2) stop("Fewer than two valid fits available for pooling.")

  Qbar <- mean(Q)
  Ubar <- mean(U)
  B <- var(Q)
  Tvar <- Ubar + (1 + 1 / m) * B
  SE_total <- sqrt(Tvar)

  if (B > 0) {
    df <- (m - 1) * (1 + Ubar / ((1 + 1 / m) * B))^2
  } else {
    df <- Inf
  }

  crit <- if (is.finite(df)) qt(0.975, df = df) else qnorm(0.975)
  CI_beta <- Qbar + c(-1, 1) * crit * SE_total
  test_stat <- Qbar / SE_total
  p_value <- if (is.finite(df)) {
    2 * pt(-abs(test_stat), df = df)
  } else {
    2 * pnorm(-abs(test_stat))
  }

  between_component <- (1 + 1 / m) * B

  c(
    n_fits = m,
    beta = Qbar,
    SE = SE_total,
    OR = exp(Qbar),
    lower95 = exp(CI_beta[1]),
    upper95 = exp(CI_beta[2]),
    p = p_value,
    df = df,
    within_var = Ubar,
    between_var = B,
    attendance_uncertainty_fraction = between_component / Tvar
  )
}

combined_breed <- pool_rubin(
  beta = breed_imp_results$beta,
  se = breed_imp_results$se
)

cat("\n=== FINAL uncertainty-propagated future-breeding effect ===\n")
print(combined_breed)
cat("Proportion of valid posterior-attendance fits with beta > 0:",
    mean(breed_imp_results$beta[is.finite(breed_imp_results$beta)] > 0), "\n")

write.csv(
  breed_imp_results,
  file.path(output_dir, "future_breeding_250_attendance_draw_fits.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(t(combined_breed)),
  file.path(output_dir, "future_breeding_pooled_result.csv"),
  row.names = FALSE
)

save(
  breed_dat, m_future_breed_dc, breed_imp_results, combined_breed, imp_idx,
  file = file.path(output_dir, "NGP_future_breeding_detection_corrected_FINAL.RData")
)

# ==============================================================================
# 8. SUBSEQUENT BREEDING SUCCESS AMONG BIRDS THAT BRED
# ==============================================================================
# This model is explicitly conditional on subsequent breeding. It must not be
# interpreted as the effect of attendance on total reproductive success among
# all cohort birds.
# ==============================================================================

success_idx <- which(
  breed_dat$future_bred == 1 & !is.na(breed_dat$future_success)
)

success_dat <- breed_dat[success_idx, ]
row.names(success_dat) <- NULL

stopifnot(nrow(success_dat) == 138)

m_future_success_dc <- glmmTMB(
  future_success ~ attendance10_dc + prev_success + experience_z + sex + winter +
    (1 | id),
  family = binomial(link = "logit"),
  data = success_dat
)

cat("\n=== Diagnostic model using posterior-mean attendance: breeding success ===\n")
print(summary(m_future_success_dc))

success_imp_results <- data.frame(
  draw = imp_idx,
  beta = NA_real_,
  se = NA_real_,
  sd_id = NA_real_,
  pdHess = NA
)

cat("\nPropagating attendance uncertainty through breeding-success model...\n")

for (k in seq_along(imp_idx)) {
  d <- imp_idx[k]
  dd <- success_dat
  dd$attendance10_dc <- A_draws[d, success_idx] * 10

  fit <- try(
    glmmTMB(
      future_success ~ attendance10_dc + prev_success + experience_z + sex + winter +
        (1 | id),
      family = binomial(link = "logit"),
      data = dd,
      control = fit_control
    ),
    silent = TRUE
  )

  if (!inherits(fit, "try-error")) {
    cf <- fixef(fit)$cond
    vc <- vcov(fit)$cond

    success_imp_results$beta[k] <- unname(cf["attendance10_dc"])
    success_imp_results$se[k] <-
      unname(sqrt(diag(vc))["attendance10_dc"])
    success_imp_results$sd_id[k] <-
      as.numeric(attr(VarCorr(fit)$cond$id, "stddev"))[1]
    success_imp_results$pdHess[k] <- isTRUE(fit$sdr$pdHess)
  }
}

combined_success <- pool_rubin(
  beta = success_imp_results$beta,
  se = success_imp_results$se
)

cat("\n=== FINAL uncertainty-propagated breeding-success effect ===\n")
print(combined_success)
cat("Proportion of valid posterior-attendance fits with beta < 0:",
    mean(success_imp_results$beta[is.finite(success_imp_results$beta)] < 0), "\n")

write.csv(
  success_imp_results,
  file.path(output_dir, "future_success_250_attendance_draw_fits.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(t(combined_success)),
  file.path(output_dir, "future_success_pooled_result.csv"),
  row.names = FALSE
)

save(
  success_dat, success_idx, m_future_success_dc,
  success_imp_results, combined_success, imp_idx,
  file = file.path(output_dir, "NGP_future_success_detection_corrected_FINAL.RData")
)

# ==============================================================================
# 9. COMPACT FINAL RESULTS TABLE
# ==============================================================================

final_results <- data.frame(
  analysis = c(
    "Previous breeding success -> winter attendance",
    "Persistent among-individual heterogeneity in attendance",
    "Winter attendance -> subsequent breeding",
    "Winter attendance -> success among subsequent breeders"
  ),
  metric = c(
    "Posterior OR (successful vs failed)",
    "Posterior SD on logit scale",
    "Pooled OR per +10 percentage points attendance",
    "Pooled OR per +10 percentage points attendance"
  ),
  estimate = c(
    prev_OR_CrI["median"],
    median(att_mat[, "sigma_id"]),
    combined_breed["OR"],
    combined_success["OR"]
  ),
  lower95 = c(
    prev_OR_CrI["lower95"],
    quantile(att_mat[, "sigma_id"], 0.025),
    combined_breed["lower95"],
    combined_success["lower95"]
  ),
  upper95 = c(
    prev_OR_CrI["upper95"],
    quantile(att_mat[, "sigma_id"], 0.975),
    combined_breed["upper95"],
    combined_success["upper95"]
  ),
  evidence = c(
    paste0("P(beta<0)=", round(prob_prev_negative, 3)),
    "Bayesian 95% credible interval",
    paste0("p=", signif(combined_breed["p"], 3)),
    paste0("p=", signif(combined_success["p"], 3))
  ),
  row.names = NULL
)

write.csv(
  final_results,
  file.path(output_dir, "FINAL_model_results_table.csv"),
  row.names = FALSE
)

cat("\n=== FINAL RESULTS TABLE ===\n")
print(final_results)

# ==============================================================================
# 10. MANUSCRIPT FIGURE: ATTENDANCE AND SUBSEQUENT REPRODUCTION
# ==============================================================================
#
# These curves use the two diagnostic GLMMs fitted to posterior-MEAN attendance.
# Therefore, the grey bands are 95% confidence intervals from those single
# downstream GLMMs; they are NOT uncertainty bands pooled over the 250 HMM
# attendance draws. Formal inference in the manuscript is based on the pooled
# results above.
#
# Population-level predictions are produced with individual random effects
# excluded (re.form = NA). Experience is held at its mean (experience_z = 0),
# previous success = 0, sex = female, and winter = 2022.
# ==============================================================================

prediction_grid <- seq(0, 1, length.out = 201)

new_breed <- data.frame(
  attendance10_dc = prediction_grid * 10,
  prev_success = 0,
  experience_z = 0,
  sex = factor("F", levels = levels(breed_dat$sex)),
  winter = factor("2022", levels = levels(breed_dat$winter)),
  id = factor(levels(breed_dat$id)[1], levels = levels(breed_dat$id))
)

pred_breed <- predict(
  m_future_breed_dc,
  newdata = new_breed,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

new_breed$fit_prob <- plogis(pred_breed$fit)
new_breed$lwr <- plogis(pred_breed$fit - 1.96 * pred_breed$se.fit)
new_breed$upr <- plogis(pred_breed$fit + 1.96 * pred_breed$se.fit)

new_success <- data.frame(
  attendance10_dc = prediction_grid * 10,
  prev_success = 0,
  experience_z = 0,
  sex = factor("F", levels = levels(success_dat$sex)),
  winter = factor("2022", levels = levels(success_dat$winter)),
  id = factor(levels(success_dat$id)[1], levels = levels(success_dat$id))
)

pred_success <- predict(
  m_future_success_dc,
  newdata = new_success,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

new_success$fit_prob <- plogis(pred_success$fit)
new_success$lwr <- plogis(pred_success$fit - 1.96 * pred_success$se.fit)
new_success$upr <- plogis(pred_success$fit + 1.96 * pred_success$se.fit)

png(
  file.path(output_dir, "Figure_subsequent_reproduction_vs_winter_attendance.png"),
  width = 2200,
  height = 1050,
  res = 220
)

op <- par(
  mfrow = c(1, 2),
  mar = c(4.5, 5, 1.5, 1),
  oma = c(3, 0, 0, 0),
  las = 1,
  xaxs = "i",
  yaxs = "i"
)

set.seed(1001)

# Panel a: subsequent breeding
yjit_breed <- pmin(
  pmax(breed_dat$future_bred + runif(nrow(breed_dat), -0.025, 0.025), 0.01),
  0.99
)

plot(
  breed_dat$attendance_dc,
  yjit_breed,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "",
  ylab = "Probability of recorded subsequent breeding",
  pch = 16,
  cex = 0.65,
  col = adjustcolor("black", alpha.f = 0.35),
  bty = "l"
)

polygon(
  c(prediction_grid, rev(prediction_grid)),
  c(new_breed$lwr, rev(new_breed$upr)),
  border = NA,
  col = "grey85"
)
lines(prediction_grid, new_breed$fit_prob, lwd = 2)
points(
  breed_dat$attendance_dc,
  yjit_breed,
  pch = 16,
  cex = 0.65,
  col = adjustcolor("black", alpha.f = 0.35)
)
mtext("a", side = 3, adj = 0, line = 0.2, font = 2)

# Panel b: breeding success among subsequent breeders
yjit_success <- pmin(
  pmax(success_dat$future_success + runif(nrow(success_dat), -0.025, 0.025), 0.01),
  0.99
)

plot(
  success_dat$attendance_dc,
  yjit_success,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "",
  ylab = "Probability of breeding success among subsequent breeders",
  pch = 16,
  cex = 0.65,
  col = adjustcolor("black", alpha.f = 0.35),
  bty = "l"
)

polygon(
  c(prediction_grid, rev(prediction_grid)),
  c(new_success$lwr, rev(new_success$upr)),
  border = NA,
  col = "grey85"
)
lines(prediction_grid, new_success$fit_prob, lwd = 2)
points(
  success_dat$attendance_dc,
  yjit_success,
  pch = 16,
  cex = 0.65,
  col = adjustcolor("black", alpha.f = 0.35)
)
mtext("b", side = 3, adj = 0, line = 0.2, font = 2)

mtext("Estimated winter attendance", side = 1, outer = TRUE, line = 1.6)

par(op)
dev.off()

# Vector version for manuscript production.
pdf(
  file.path(output_dir, "Figure_subsequent_reproduction_vs_winter_attendance.pdf"),
  width = 10,
  height = 4.8
)

op <- par(
  mfrow = c(1, 2),
  mar = c(4.5, 5, 1.5, 1),
  oma = c(3, 0, 0, 0),
  las = 1,
  xaxs = "i",
  yaxs = "i"
)

plot(
  breed_dat$attendance_dc, yjit_breed,
  xlim = c(0, 1), ylim = c(0, 1), xlab = "",
  ylab = "Probability of recorded subsequent breeding",
  pch = 16, cex = 0.65,
  col = adjustcolor("black", alpha.f = 0.35), bty = "l"
)
polygon(
  c(prediction_grid, rev(prediction_grid)),
  c(new_breed$lwr, rev(new_breed$upr)),
  border = NA, col = "grey85"
)
lines(prediction_grid, new_breed$fit_prob, lwd = 2)
points(
  breed_dat$attendance_dc, yjit_breed,
  pch = 16, cex = 0.65,
  col = adjustcolor("black", alpha.f = 0.35)
)
mtext("a", side = 3, adj = 0, line = 0.2, font = 2)

plot(
  success_dat$attendance_dc, yjit_success,
  xlim = c(0, 1), ylim = c(0, 1), xlab = "",
  ylab = "Probability of breeding success among subsequent breeders",
  pch = 16, cex = 0.65,
  col = adjustcolor("black", alpha.f = 0.35), bty = "l"
)
polygon(
  c(prediction_grid, rev(prediction_grid)),
  c(new_success$lwr, rev(new_success$upr)),
  border = NA, col = "grey85"
)
lines(prediction_grid, new_success$fit_prob, lwd = 2)
points(
  success_dat$attendance_dc, yjit_success,
  pch = 16, cex = 0.65,
  col = adjustcolor("black", alpha.f = 0.35)
)
mtext("b", side = 3, adj = 0, line = 0.2, font = 2)

mtext("Estimated winter attendance", side = 1, outer = TRUE, line = 1.6)

par(op)
dev.off()

# ==============================================================================
# 11. SAVE FINAL WORKSPACE COMPONENTS AND SESSION INFORMATION
# ==============================================================================

# Save only the analysis objects needed to reproduce tables/figures; do not use
# save.image(), because compiled NIMBLE objects can make the workspace very large.
save(
  bird_cov, check22, check23, x22, x23,
  re_samples_att, att_mat,
  A_draws, attendance_summary,
  breed_dat, success_dat,
  m_future_breed_dc, m_future_success_dc,
  breed_imp_results, success_imp_results,
  combined_breed, combined_success,
  final_results,
  file = file.path(output_dir, "NGP_FINAL_ANALYSIS_OBJECTS.RData")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "R_sessionInfo.txt")
)

cat("\n============================================================\n")
cat("FINAL ANALYSIS COMPLETE\n")
cat("Outputs written to:", normalizePath(output_dir), "\n")
cat("============================================================\n")










# ==============================================================================
# 12. ALTERNATIVE JOLLY-SEBER RECAPTURE MODEL
# ==============================================================================


GP_JS <- nimbleCode({
  # -------------------------------------------------
  # Parameters:
  # phi: annual survival probability
  # p: recapture probability
  # N[i]: size of observed population in year i
  # gamma: probability to enter marked population
  # -------------------------------------------------
  
  # Priors and constraints
  for (i in 1:M) {
    for (t in 1:(n.years - 1)) {
      phi[i, t] <- mean.phi       ## this can be linked to the multi-state model phi
    } #t
    for (t in 1:n.years) {
      #expanded in effort corrected model
      for (j in 1:n.occ) {
        #expanded in effort corrected model for occasion specific capture probability
        logit(p[i, j, t]) <- mean.lp + b.effort * hours[j, t]
      } #j
    } #t
  } #i
  
  mean.phi ~ dbeta(80,15)              # Prior for mean adult survival - should be higher than 70%
  mean.lp <- log(mean.p / (1 - mean.p))
  mean.p ~ dunif(0, 1)                # Prior for mean capture
  b.effort ~ dunif(-10, 10)
  
  for (t in 1:n.years) {
    gamma[t] ~ dbeta(1.5, 25)           ## THIS IS THE PROBABILITY TO ENTER THE MARKED POPULATION - either fix or use informative prior
  } #t
  
  # Likelihood
  for (i in 1:M) {
    # First occasion
    # State process
    z[i, 1] ~ dbern(gamma[1])
    
    # Observation process
    for (j in 1:n.occ) {
      mu1[i, j] <- z[i, 1] * p[i, j, 1]
      y[i, j, 1] ~ dbern(mu1[i, j])
    } #j
    
    
    # Subsequent occasions
    for (t in 2:n.years) {
      # State process
      q[i, t - 1] <- 1 - z[i, t - 1]		# Availability for recruitment
      mu2[i, t] <- phi[i, t - 1] * z[i, t - 1] + gamma[t] * prod(q[i, 1:(t - 1)])
      z[i, t] ~ dbern(mu2[i, t])
      
      # Observation process
      for (j in 1:n.occ) {
        mu3[i, j, t] <- z[i, t] * p[i, j, t]
        y[i, j, t] ~ dbern(mu3[i, j, t])
      } #j
    } #t
  } #i
  
  # Calculate derived population parameters
  for (t in 1:n.years) {
    N[t] <- sum(z[1:M, t])        # Actual population size
    prop.observed[t] <- N[t]/M
  } #t
})  #### END OF THE MODEL STATEMENT
