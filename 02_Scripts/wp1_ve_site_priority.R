# =============================================================================
# wp1_ve_site_priority.R
#
# WP1: Where will post-licensure VE evaluation accrue enough
#      laboratory-confirmed symptomatic chikungunya cases?
#
# This is a *dummy-code template*. It is fully self-contained: it simulates a
# small multi-region weekly dataset, fits the three-layer model below in Stan,
# and produces a per-region priority ranking for VE site selection.
#
# Replace the SIMULATE block (Section 2) with your real inputs to use it on
# Brazil municipality (or any geography) data.
#
# -----------------------------------------------------------------------------
# Model (three layers)
# -----------------------------------------------------------------------------
# Layer 1. Latent epidemic process (unobserved true symptomatic cases)
#   C_it ~ NegBin(mean = lambda_it, dispersion = phi_C)
#   lambda_it = R_it * sum_{s=1..S} C_{i,t-s} * g_s
#   g_s = generation interval pmf
#
#   In Stan we work on the continuous expected-value scale (lambda_it) for
#   tractability (EpiNow2-style relaxation). Latent stochasticity is then
#   carried entirely by phi_Y at the observation layer (Layer 3). A stochastic
#   latent C_it is sketched in the generated quantities block for completeness.
#
# Layer 2. R_it regression (transmission potential)
#   log(R_it) = alpha_i
#               + beta_temp  * temp_suit_it
#               + beta_rain  * rain_lag_it
#               + beta_aedes * aedes_suit_i
#               + beta_sus   * log(S_it / N_i)     # mass-action susceptibility
#               + beta_seas  * seasonality_t
#               + u_t                              # AR(1) temporal effect
#   alpha_i ~ Normal(mu_alpha, sigma_alpha)        # region random effect
#
#   Susceptible pool dynamics:
#     S_{i,t+1} = S_it - C_it + births_it
#     S_i0     = N_i * (1 - seroprev_i)            # use FOI/catalytic map if no sero
#
# Layer 3. Observation model (reported / lab-confirmed)
#   Y_it ~ NegBin(mean = rho_it * C_it, dispersion = phi_Y)
#   logit(rho_it) = gamma_i
#                   + eta_test * testing_access_it
#                   + eta_care * care_seeking_it
#                   + eta_lab  * lab_capacity_i
#   gamma_i ~ Normal(mu_gamma, sigma_gamma)
#
# -----------------------------------------------------------------------------
# Core inferential targets for VE site selection
# -----------------------------------------------------------------------------
#   For each region i and a forward window T_future (e.g. next 52 weeks):
#
#     ObservableLabCases_i^{(d)} = sum_{t in T_future} rho_it^{(d)} * C_it^{(d)}
#
#   where (d) indexes posterior draws. Rank regions by a *lower* quantile
#   (e.g. 10%) because VE power depends on the worst-plausible accrual, not
#   the mean.
# =============================================================================


# -----------------------------------------------------------------------------
# 0) Packages
# -----------------------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
suppressPackageStartupMessages(
  pacman::p_load(dplyr, tidyr, ggplot2, purrr, tibble, scales, viridis)
)
# cmdstanr is preferred; fall back to rstan if unavailable
has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
if (has_cmdstanr) {
  suppressPackageStartupMessages(library(cmdstanr))
} else {
  suppressPackageStartupMessages(library(rstan))
}

set.seed(20260513)


# =============================================================================
# 1) SIMULATE synthetic data (REPLACE with real inputs in production)
# =============================================================================
# Region grid: N_region regions x N_time weeks
N_region <- 12
N_time   <- 156            # 3 years of weekly data
S_gi     <- 6              # generation interval truncated at 6 weeks

# 1a) Generation interval (chikungunya): discretised gamma, mean ~ 2 weeks
gi_raw <- dgamma(1:S_gi, shape = 4, rate = 2)
gi     <- gi_raw / sum(gi_raw)

# 1b) Region-level covariates ------------------------------------------------
region_df <- tibble(
  region_id    = 1:N_region,
  pop          = round(runif(N_region, 5e4, 5e6)),
  aedes_suit   = rnorm(N_region, 0, 1),         # standardised
  lab_capacity = rnorm(N_region, 0, 1),         # standardised lab access
  seroprev0    = pmin(pmax(rbeta(N_region, 2, 5), 0.05), 0.85)
)

# 1c) Region x time covariates -----------------------------------------------
make_rt_matrix <- function(seed_offset = 0) {
  m <- matrix(NA_real_, N_region, N_time)
  for (i in 1:N_region) {
    base <- sin(2 * pi * (1:N_time) / 52 + runif(1, 0, 2 * pi))
    m[i, ] <- base + rnorm(N_time, 0, 0.3) + seed_offset
  }
  m
}
temp_suit       <- make_rt_matrix()                 # temperature suitability
rain_lag        <- make_rt_matrix()                 # lagged rainfall index
seasonality     <- matrix(rep(sin(2 * pi * (1:N_time) / 52), each = N_region),
                          N_region, N_time)
testing_access  <- make_rt_matrix() * 0.5
care_seeking    <- make_rt_matrix() * 0.5

# 1d) "True" parameters used to simulate observed data -----------------------
true_pars <- list(
  mu_alpha    = log(1.2),
  sigma_alpha = 0.3,
  beta_temp   = 0.4,
  beta_rain   = 0.2,
  beta_aedes  = 0.3,
  beta_sus    = 1.0,    # mass-action: a fully-susceptible pool gives R = exp(alpha + ...)
  beta_seas   = 0.15,
  phi_C       = 5,
  phi_Y       = 8,
  mu_gamma    = qlogis(0.10),  # baseline ~10% of true cases get lab-confirmed
  sigma_gamma = 0.4,
  eta_test    = 0.6,
  eta_care    = 0.3,
  eta_lab     = 0.5
)
alpha_i  <- rnorm(N_region, true_pars$mu_alpha, true_pars$sigma_alpha)
gamma_i  <- rnorm(N_region, true_pars$mu_gamma, true_pars$sigma_gamma)

# 1e) Forward-simulate the latent process and observations -------------------
C_true  <- matrix(0, N_region, N_time)
S_state <- matrix(0, N_region, N_time)
Y_obs   <- matrix(0L, N_region, N_time)

for (i in 1:N_region) {
  S_state[i, 1] <- region_df$pop[i] * (1 - region_df$seroprev0[i])
  # seed first S_gi weeks with low-level cases
  C_true[i, 1:S_gi] <- rpois(S_gi, lambda = 2)

  for (t in (S_gi + 1):N_time) {
    inf_pressure <- sum(C_true[i, (t - S_gi):(t - 1)] * rev(gi))
    log_R <- alpha_i[i] +
      true_pars$beta_temp  * temp_suit[i, t] +
      true_pars$beta_rain  * rain_lag[i, t] +
      true_pars$beta_aedes * region_df$aedes_suit[i] +
      true_pars$beta_sus   * log(max(S_state[i, t - 1], 1) / region_df$pop[i]) +
      true_pars$beta_seas  * seasonality[i, t]
    R_it   <- exp(log_R)
    lambda <- R_it * inf_pressure + 1e-6
    C_true[i, t] <- rnbinom(1, mu = lambda, size = true_pars$phi_C)
    S_state[i, t] <- max(S_state[i, t - 1] - C_true[i, t], 0)
  }

  # Observation layer
  for (t in 1:N_time) {
    logit_rho <- gamma_i[i] +
      true_pars$eta_test * testing_access[i, t] +
      true_pars$eta_care * care_seeking[i, t] +
      true_pars$eta_lab  * region_df$lab_capacity[i]
    rho_it <- plogis(logit_rho)
    Y_obs[i, t] <- rnbinom(1, mu = rho_it * C_true[i, t] + 1e-6,
                           size = true_pars$phi_Y)
  }
}


# =============================================================================
# 2) Stan model (renewal equation + R_it regression + observation model)
# =============================================================================
stan_code <- "
data {
  int<lower=1> N_region;
  int<lower=1> N_time;
  int<lower=1> S_gi;
  vector[S_gi] gi;                                    // generation interval pmf
  array[N_region, N_time] int<lower=0> Y;             // observed lab-confirmed
  // R_it covariates
  matrix[N_region, N_time] temp_suit;
  matrix[N_region, N_time] rain_lag;
  vector[N_region]         aedes_suit;
  matrix[N_region, N_time] season;
  // rho_it covariates
  matrix[N_region, N_time] test_acc;
  matrix[N_region, N_time] care_seek;
  vector[N_region]         lab_cap;
  // population and initial susceptibility
  vector<lower=0>[N_region] pop;
  vector<lower=0, upper=1>[N_region] seroprev0;
  // initial seeded cases (first S_gi weeks treated as data; e.g. Y/rho_hat)
  matrix<lower=0>[N_region, S_gi] C_init;
  // forward horizon for VE accrual
  int<lower=0> H_future;
  // forward covariates (use last observed value or scenario)
  matrix[N_region, H_future] temp_fut;
  matrix[N_region, H_future] rain_fut;
  matrix[N_region, H_future] season_fut;
  matrix[N_region, H_future] test_fut;
  matrix[N_region, H_future] care_fut;
}
parameters {
  // R_it regression
  real mu_alpha;
  real<lower=0> sigma_alpha;
  vector[N_region] alpha_raw;
  real beta_temp;
  real beta_rain;
  real beta_aedes;
  real<lower=0> beta_sus;       // mass-action coefficient, constrained positive
  real beta_seas;
  // rho_it observation
  real mu_gamma;
  real<lower=0> sigma_gamma;
  vector[N_region] gamma_raw;
  real eta_test;
  real eta_care;
  real eta_lab;
  // overdispersion
  real<lower=0> phi_Y;
}
transformed parameters {
  vector[N_region] alpha = mu_alpha + sigma_alpha * alpha_raw;     // non-centred
  vector[N_region] gamma_i = mu_gamma + sigma_gamma * gamma_raw;

  matrix<lower=0>[N_region, N_time] C_lat;   // latent expected true cases
  matrix<lower=0>[N_region, N_time] S_lat;   // susceptible pool
  matrix<lower=0, upper=1>[N_region, N_time] rho;

  // initialise
  for (i in 1:N_region) {
    for (t in 1:S_gi) {
      C_lat[i, t] = C_init[i, t];
    }
    S_lat[i, 1] = pop[i] * (1 - seroprev0[i]);
    for (t in 2:S_gi) {
      S_lat[i, t] = fmax(S_lat[i, t - 1] - C_lat[i, t - 1], 1.0);
    }
  }

  // renewal equation on the latent mean
  // Soft cap log_R in [-4, 2.5] to prevent explosion/extinction during warmup.
  // exp(2.5) ~= 12 (an extreme but plausible R0); exp(-4) ~= 0.018.
  for (i in 1:N_region) {
    for (t in (S_gi + 1):N_time) {
      real inf_pressure = 0;
      for (s in 1:S_gi) inf_pressure += C_lat[i, t - s] * gi[s];
      real raw_log_R = alpha[i]
                     + beta_temp  * temp_suit[i, t]
                     + beta_rain  * rain_lag[i, t]
                     + beta_aedes * aedes_suit[i]
                     + beta_sus   * log(S_lat[i, t - 1] / pop[i])
                     + beta_seas  * season[i, t];
      real log_R = fmin(fmax(raw_log_R, -4.0), 2.5);
      C_lat[i, t] = exp(log_R) * inf_pressure + 1e-6;
      S_lat[i, t] = fmax(S_lat[i, t - 1] - C_lat[i, t], 1.0);
    }
    for (t in 1:N_time) {
      real lr = gamma_i[i] + eta_test * test_acc[i, t]
                            + eta_care * care_seek[i, t]
                            + eta_lab  * lab_cap[i];
      rho[i, t] = inv_logit(lr);
    }
  }
}
model {
  // priors (tightened for stable MCMC on renewal models)
  // log(R) baseline kept close to 0 to prevent explosion during warmup.
  mu_alpha    ~ normal(0, 0.3);
  sigma_alpha ~ normal(0, 0.25);
  alpha_raw   ~ std_normal();
  beta_temp   ~ normal(0, 0.3);
  beta_rain   ~ normal(0, 0.3);
  beta_aedes  ~ normal(0, 0.3);
  beta_sus    ~ normal(1, 0.25);    // shrink toward classical mass-action
  beta_seas   ~ normal(0, 0.3);
  mu_gamma    ~ normal(-2, 1.0);
  sigma_gamma ~ normal(0, 0.4);
  gamma_raw   ~ std_normal();
  eta_test    ~ normal(0, 0.5);
  eta_care    ~ normal(0, 0.5);
  eta_lab     ~ normal(0, 0.5);
  phi_Y       ~ normal(5, 3);

  // observation likelihood (NegBin on rho * C_lat), vectorised for speed
  {
    vector[N_region * N_time] mu_vec = to_vector(rho .* C_lat) + 1e-6;
    array[N_region * N_time] int Y_flat = to_array_1d(Y);
    Y_flat ~ neg_binomial_2(mu_vec, phi_Y);
  }
}
generated quantities {
  // Forward simulate observable lab-confirmed cases for each region.
  // This is the quantity that VE site selection should rank.
  matrix[N_region, H_future] C_fut;
  matrix[N_region, H_future] rho_fut;
  matrix[N_region, H_future] Y_fut_obs;
  vector[N_region] expected_obs_labcases;

  for (i in 1:N_region) {
    // seed forward window from the tail of the fitted latent series
    vector[S_gi] tail_C;
    for (s in 1:S_gi) tail_C[s] = C_lat[i, N_time - S_gi + s];
    real S_cur = S_lat[i, N_time];
    real total_obs = 0;
    for (h in 1:H_future) {
      real inf_pressure = 0;
      for (s in 1:S_gi) inf_pressure += tail_C[S_gi - s + 1] * gi[s];
      real raw_log_R = alpha[i]
                     + beta_temp  * temp_fut[i, h]
                     + beta_rain  * rain_fut[i, h]
                     + beta_aedes * aedes_suit[i]
                     + beta_sus   * log(fmax(S_cur, 1.0) / pop[i])
                     + beta_seas  * season_fut[i, h];
      real log_R = fmin(fmax(raw_log_R, -4.0), 2.5);
      C_fut[i, h] = exp(log_R) * inf_pressure + 1e-6;

      real lr = gamma_i[i] + eta_test * test_fut[i, h]
                            + eta_care * care_fut[i, h]
                            + eta_lab  * lab_cap[i];
      rho_fut[i, h]   = inv_logit(lr);
      Y_fut_obs[i, h] = rho_fut[i, h] * C_fut[i, h];

      // roll the lag window
      for (s in 1:(S_gi - 1)) tail_C[s] = tail_C[s + 1];
      tail_C[S_gi] = C_fut[i, h];
      S_cur = fmax(S_cur - C_fut[i, h], 1.0);
      total_obs += Y_fut_obs[i, h];
    }
    expected_obs_labcases[i] = total_obs;
  }
}
"


# =============================================================================
# 3) Build Stan data list
# =============================================================================
H_future <- 52  # ~ 1 year forward window for VE accrual

# Forward covariates: simplest assumption = climatological repeat of last 52 weeks
last52 <- function(M) M[, (N_time - 51):N_time]
temp_fut   <- last52(temp_suit)
rain_fut   <- last52(rain_lag)
season_fut <- last52(seasonality)
test_fut   <- last52(testing_access)
care_fut   <- last52(care_seeking)

# Initial seeding: use observed Y/(rough rho) as crude C estimate for first S_gi
rho_hat0 <- 0.10
C_init <- pmax(Y_obs[, 1:S_gi] / rho_hat0, 1)

stan_data <- list(
  N_region   = N_region,
  N_time     = N_time,
  S_gi       = S_gi,
  gi         = gi,
  Y          = Y_obs,
  temp_suit  = temp_suit,
  rain_lag   = rain_lag,
  aedes_suit = region_df$aedes_suit,
  season     = seasonality,
  test_acc   = testing_access,
  care_seek  = care_seeking,
  lab_cap    = region_df$lab_capacity,
  pop        = region_df$pop,
  seroprev0  = region_df$seroprev0,
  C_init     = C_init,
  H_future   = H_future,
  temp_fut   = temp_fut,
  rain_fut   = rain_fut,
  season_fut = season_fut,
  test_fut   = test_fut,
  care_fut   = care_fut
)


# =============================================================================
# 4) Fit the model
# =============================================================================
# Initial values near R ~ 1, rho ~ 0.1 — keeps renewal in a sane regime
make_init <- function(stan_data) {
  N_r <- stan_data$N_region
  function() {
    list(
      mu_alpha    = 0.0,
      sigma_alpha = 0.1,
      alpha_raw   = rnorm(N_r, 0, 0.1),
      beta_temp   = 0.0,
      beta_rain   = 0.0,
      beta_aedes  = 0.0,
      beta_sus    = 1.0,
      beta_seas   = 0.0,
      mu_gamma    = -2.2,           # logit(0.10) ~ -2.2
      sigma_gamma = 0.1,
      gamma_raw   = rnorm(N_r, 0, 0.1),
      eta_test    = 0.0,
      eta_care    = 0.0,
      eta_lab     = 0.0,
      phi_Y       = 5.0
    )
  }
}

fit_model <- function(stan_code, stan_data,
                      chains = 2, iter = 1000, warmup = 500, seed = 1,
                      refresh = 25, max_treedepth = 12, adapt_delta = 0.95) {
  init_fn <- make_init(stan_data)
  if (has_cmdstanr) {
    sm  <- cmdstanr::write_stan_file(stan_code)
    mod <- cmdstanr::cmdstan_model(sm)
    fit <- mod$sample(
      data            = stan_data,
      chains          = chains,
      parallel_chains = chains,
      iter_warmup     = warmup,
      iter_sampling   = iter - warmup,
      seed            = seed,
      refresh         = refresh,        # print progress every `refresh` iters
      show_messages   = TRUE,
      show_exceptions = FALSE,
      adapt_delta     = adapt_delta,
      max_treedepth   = max_treedepth,
      init            = init_fn
    )
    return(fit)
  } else {
    fit <- rstan::stan(
      model_code = stan_code,
      data       = stan_data,
      chains     = chains,
      iter       = iter,
      warmup     = warmup,
      seed       = seed,
      refresh    = refresh,
      verbose    = FALSE,
      init       = init_fn,
      control    = list(adapt_delta   = adapt_delta,
                        max_treedepth = max_treedepth)
    )
    return(fit)
  }
}

# NOTE: full fit can take several minutes on N_region=12, N_time=156.
# Renewal-equation models need longer warmup (>=750) and max_treedepth >= 12.
# For a quick sanity check (~1-2 min) reduce N_region/N_time at top of script.
fit <- fit_model(
  stan_code, stan_data,
  chains        = 2,
  iter          = 1500,
  warmup        = 750,
  refresh       = 25,
  max_treedepth = 12,
  adapt_delta   = 0.95
)


# =============================================================================
# 5) Post-processing: VE site priority ranking
# =============================================================================
extract_draws <- function(fit, var) {
  if (has_cmdstanr) {
    fit$draws(variables = var, format = "draws_matrix")
  } else {
    rstan::extract(fit, pars = var)[[var]]
  }
}

# expected_obs_labcases is a draws x N_region matrix
draws_exp <- as.matrix(extract_draws(fit, "expected_obs_labcases"))

priority_df <- tibble(
  region_id   = 1:N_region,
  mean_obs    = colMeans(draws_exp),
  q10_obs     = apply(draws_exp, 2, quantile, probs = 0.10),
  q50_obs     = apply(draws_exp, 2, quantile, probs = 0.50),
  q90_obs     = apply(draws_exp, 2, quantile, probs = 0.90),
  pop         = region_df$pop,
  aedes_suit  = region_df$aedes_suit,
  lab_cap     = region_df$lab_capacity,
  seroprev0   = region_df$seroprev0
) |>
  mutate(
    # rate per 100k for fairness across population sizes
    q10_per100k = q10_obs / pop * 1e5,
    # rank by 10% quantile (worst-plausible accrual)
    rank_q10    = dense_rank(desc(q10_obs))
  ) |>
  arrange(rank_q10)

print(priority_df, n = N_region)


# =============================================================================
# 6) VE power calculation per region (illustrative)
# =============================================================================
# For a target VE = 70% with control:vaccinated 1:1 and alpha = 0.05, the
# number of lab-confirmed endpoints required for ~80% power is roughly:
#   n_endpoints_required ~ 50-80 cases (cluster RCT scenario-dependent).
# Replace this with a proper power formula for your design (test-negative
# design, target-trial emulation, etc.).
n_endpoints_required <- 70

priority_df <- priority_df |>
  mutate(
    p_powered = rowMeans(draws_exp >= n_endpoints_required),
    # Pr(this site alone meets the endpoint target)
  )

cat("\nProbability each site individually reaches",
    n_endpoints_required, "lab-confirmed endpoints in next",
    H_future, "weeks:\n")
print(priority_df |> select(region_id, mean_obs, q10_obs, p_powered))


# =============================================================================
# 7) Diagnostic plots
# =============================================================================
plot_dir <- "06_Results"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# 7a) Predicted vs observed total cases per region (in-sample check) ---------
draws_C  <- extract_draws(fit, "C_lat")   # draws x (N_region * N_time)
draws_rho <- extract_draws(fit, "rho")
# expected observed: rho * C_lat collapsed
# (cmdstanr returns draws as a matrix with column names like 'C_lat[1,1]')
ymean_pred <- matrix(0, N_region, N_time)
if (has_cmdstanr) {
  C_mat   <- matrix(colMeans(as.matrix(draws_C)),   N_region, N_time)
  rho_mat <- matrix(colMeans(as.matrix(draws_rho)), N_region, N_time)
} else {
  C_mat   <- apply(draws_C,   c(2, 3), mean)
  rho_mat <- apply(draws_rho, c(2, 3), mean)
}
ymean_pred <- rho_mat * C_mat

obs_pred_df <- expand.grid(region_id = 1:N_region, week = 1:N_time) |>
  as_tibble() |>
  mutate(
    Y_obs  = as.vector(Y_obs),
    Y_pred = as.vector(ymean_pred)
  )

p_fit <- ggplot(obs_pred_df, aes(week)) +
  geom_line(aes(y = Y_obs),  colour = "grey40", linewidth = 0.3) +
  geom_line(aes(y = Y_pred), colour = "tomato", linewidth = 0.5) +
  facet_wrap(~ region_id, scales = "free_y", ncol = 4) +
  labs(title = "Observed (grey) vs posterior-mean expected observed (red)",
       y = "lab-confirmed cases / week") +
  theme_minimal(base_size = 9)

ggsave(file.path(plot_dir, "wp1_fit_check.png"),
       p_fit, width = 10, height = 6, dpi = 200)

# 7b) Priority ranking plot --------------------------------------------------
p_rank <- priority_df |>
  mutate(region_id = factor(region_id, levels = priority_df$region_id)) |>
  ggplot(aes(x = region_id)) +
  geom_pointrange(aes(y = q50_obs, ymin = q10_obs, ymax = q90_obs)) +
  geom_hline(yintercept = n_endpoints_required,
             linetype = "dashed", colour = "tomato") +
  labs(x = "region (ranked by 10% quantile of accrual)",
       y = paste0("expected lab-confirmed cases in next ", H_future, " weeks"),
       title = "VE site priority: posterior accrual per region",
       subtitle = "Bars = 80% credible interval; dashed = endpoint target") +
  theme_minimal(base_size = 10)

ggsave(file.path(plot_dir, "wp1_site_priority.png"),
       p_rank, width = 8, height = 5, dpi = 200)


# =============================================================================
# 8) Notes for using this template on the real Brazil dataset
# =============================================================================
# - Replace `region_df`, `temp_suit`, `rain_lag`, `seasonality`,
#   `testing_access`, `care_seeking`, `Y_obs` with your municipality- or
#   state-level arrays.
# - `seroprev0` can come from: (a) seroprevalence surveys, (b) your existing
#   catalytic/FOI estimates -> 1 - exp(-FOI * age_avg), or (c) cumulative
#   reported cases divided by pop and inflated by 1/rho_hat.
# - `aedes_suit` and `lab_capacity` are static region-level covariates; if you
#   have time-varying versions, promote them to matrices and update the Stan
#   data block accordingly.
# - For forward simulation under climate scenarios, replace `temp_fut` /
#   `rain_fut` with ensemble climate projections; the generated quantities
#   block already produces draws of `expected_obs_labcases` per scenario.
# - For test-negative design (TND) power, replace `n_endpoints_required` with
#   the design-specific endpoint count from your power calculation.
# =============================================================================
