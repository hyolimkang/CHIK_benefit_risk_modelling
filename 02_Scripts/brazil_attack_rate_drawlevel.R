## Draw-keyed attack-rate construction (baseline stratum aligned)
## ------------------------------------------------------------------
## Assumes setup.R has been sourced and the following objects exist:
## sim_results_vc_ixchiq_model, all_draws_ix_true, lhs_sample_young,
## pop_by_state, foi_daily_by_state

# Baseline stratum used to avoid duplicated pre-burden rows across VE/Coverage/Scenario
baseline_cov <- "cov50"
baseline_ve <- "VE0"
baseline_scenario <- 1L

# Primary AR pathway: linear I/S0 with a safe cap for probability interfaces.
use_linear_ar_primary <- TRUE
ar_linear_cap <- 0.99
# Default majority rule used before data-driven threshold fitting.
prob_majority_pct <- 50
prob_round_digits <- 0

# If TRUE, use legacy-compatible fixed symptomatic fraction for downstream summaries.
use_fixed_symp_const <- TRUE
symp_const_fixed <- 0.524

get_total_S0_by_stratum <- function(sim_results, t0 = 1) {
  out <- data.frame()

  for (region in names(sim_results)) {
    for (ve in names(sim_results[[region]])) {
      for (vc in names(sim_results[[region]][[ve]])) {
        scenario_list <- sim_results[[region]][[ve]][[vc]]
        for (sc in seq_along(scenario_list)) {
          scen <- scenario_list[[sc]]

          # Preferred: draw-level S0 saved in sim_result$total_S0_by_draw
          if (!is.null(scen$sim_result$total_S0_by_draw)) {
            s0_vec <- as.numeric(scen$sim_result$total_S0_by_draw)
            out <- rbind(
              out,
              data.frame(
                region = region,
                VE = ve,
                VC = vc,
                Scenario = as.integer(sc),
                draw_id = seq_along(s0_vec),
                t0 = t0,
                total_S0 = s0_vec
              )
            )
          } else {
            # Backward compatibility: older objects only keep last-draw S matrix
            S <- scen$sim_out$S
            total_S0 <- sum(S[, t0], na.rm = TRUE)
            out <- rbind(
              out,
              data.frame(
                region = region,
                VE = ve,
                VC = vc,
                Scenario = as.integer(sc),
                draw_id = NA_integer_,
                t0 = t0,
                total_S0 = total_S0
              )
            )
          }
        }
      }
    }
  }

  out
}

S0_strata <- get_total_S0_by_stratum(sim_results_vc_ixchiq_model, t0 = 1) %>%
  dplyr::filter(
    VE == baseline_ve,
    VC == baseline_cov,
    Scenario == baseline_scenario
  ) %>%
  dplyr::select(region, draw_id, total_S0)

# If duplicates remain by region x draw_id, keep first and warn.
dup_region_draw <- S0_strata %>%
  dplyr::count(region, draw_id, name = "n_rows") %>%
  dplyr::filter(n_rows > 1)
if (nrow(dup_region_draw) > 0) {
  warning("Multiple baseline S0 rows per region/draw_id found; keeping first row per key.")
}

S0_lookup <- S0_strata %>%
  dplyr::group_by(region, draw_id) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

build_drawkey_preburden <- function(all_draws_df, baseline_cov, baseline_ve, baseline_scenario) {
  pre_check <- all_draws_df %>%
    dplyr::group_by(.data$Region, .data$draw_id) %>%
    dplyr::summarise(n_pre = dplyr::n_distinct(.data$total_pre), .groups = "drop")

  if (any(pre_check$n_pre != 1L)) {
    stop("Non-unique total_pre detected for some Region/draw_id; cannot build deterministic pre-burden.")
  }

  row_key <- all_draws_df %>%
    dplyr::mutate(
      scenario_chr = as.character(.data$Scenario),
      is_baseline = .data$Coverage == baseline_cov &
        .data$VE == baseline_ve &
        .data$scenario_chr == as.character(baseline_scenario)
    ) %>%
    dplyr::arrange(.data$Region, .data$draw_id, dplyr::desc(.data$is_baseline), .data$Coverage, .data$VE, .data$scenario_chr) %>%
    dplyr::group_by(.data$Region, .data$draw_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  pre_key <- row_key %>%
    dplyr::transmute(Region, draw_id, total_pre_key = .data$total_pre)

  if ("rho" %in% names(row_key)) {
    rho_key <- row_key %>%
      dplyr::transmute(Region, draw_id, rho_key = as.numeric(.data$rho))
  } else if ("total_pre_true" %in% names(row_key)) {
    rho_key <- row_key %>%
      dplyr::transmute(
        Region,
        draw_id,
        rho_key = .data$total_pre / pmax(as.numeric(.data$total_pre_true), 1e-12)
      )
  } else {
    stop("Neither rho nor total_pre_true available to reconstruct draw-key pre-burden.")
  }

  out <- pre_key %>%
    dplyr::left_join(rho_key, by = c("Region", "draw_id")) %>%
    dplyr::mutate(total_pre_true_key = .data$total_pre_key / .data$rho_key)

  if (any(!is.finite(out$total_pre_true_key) | out$total_pre_true_key <= 0, na.rm = TRUE)) {
    stop("Invalid total_pre_true_key generated for some Region/draw_id.")
  }

  out
}

## Deterministic draw-keyed symptomatic burden (rho-adjusted)
drawkey_pre <- build_drawkey_preburden(
  all_draws_df = all_draws_ix_true,
  baseline_cov = baseline_cov,
  baseline_ve = baseline_ve,
  baseline_scenario = baseline_scenario
)

tot_symp_by_region_draw <- drawkey_pre %>%
  dplyr::transmute(Region, draw_id, tot_symp = .data$total_pre_true_key)

## Convert symptomatic burden to inferred infections using draw-wise symptomatic fraction
symp_by_draw <- lhs_sample_young %>%
  dplyr::transmute(
    draw_id = dplyr::row_number(),
    symp_prop = as.numeric(symp_overall)
  ) %>%
  dplyr::mutate(symp_prop = pmin(pmax(symp_prop, 1e-6), 1 - 1e-6))

true_inf_by_region_draw <- tot_symp_by_region_draw %>%
  dplyr::left_join(symp_by_draw, by = "draw_id") %>%
  dplyr::mutate(tot_inf = tot_symp / symp_prop)

## Draw-keyed S0 table: region baseline S0 explicitly keyed to each draw.
## (S0 is region-level from the selected baseline stratum; keyed by draw to avoid
## accidental many-to-many joins and to keep arithmetic explicitly draw-indexed.)
S0_by_region_draw <- true_inf_by_region_draw %>%
  dplyr::distinct(Region, draw_id) %>%
  dplyr::left_join(S0_lookup, by = c("Region" = "region", "draw_id" = "draw_id"))

# Fallback for older sim_results objects that do not contain draw-level S0.
if (any(is.na(S0_by_region_draw$total_S0))) {
  S0_region_fallback <- S0_lookup %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(total_S0 = dplyr::first(total_S0), .groups = "drop")
  S0_by_region_draw <- S0_by_region_draw %>%
    dplyr::select(-total_S0) %>%
    dplyr::left_join(S0_region_fallback, by = c("Region" = "region"))
}

if (any(is.na(S0_by_region_draw$total_S0))) {
  stop("Missing baseline total_S0 for some Region/draw_id after stratum filter.")
}

AR_draw_by_region <- true_inf_by_region_draw %>%
  dplyr::left_join(S0_by_region_draw, by = c("Region", "draw_id")) %>%
  dplyr::left_join(pop_by_state, by = c("Region" = "region")) %>%
  dplyr::mutate(
    # Linear I/S0
    AR_lin_S0 = .data$tot_inf / .data$total_S0,
    AR_lin_S0_pct = 100 * .data$AR_lin_S0,
    flag_inf_gt_S0 = .data$tot_inf > .data$total_S0,
    # Bounded probability transform
    H_S0 = .data$AR_lin_S0,
    AR_S0 = 1 - exp(-.data$H_S0),
    AR_S0_pct = 100 * .data$AR_S0,
    gap_H_minus_AR_S0 = .data$H_S0 - .data$AR_S0,
    # Population denominator variants
    AR_lin_pop = .data$tot_inf / .data$tot_pop,
    AR_lin_pop_pct = 100 * .data$AR_lin_pop,
    AR_pop = 1 - exp(-.data$AR_lin_pop),
    AR_pop_pct = 100 * .data$AR_pop
  )

## Draw lists
eps <- 1e-12
draws_list_lin_S0 <- split(AR_draw_by_region$AR_lin_S0, AR_draw_by_region$Region)
draws_list_lin_S0 <- lapply(draws_list_lin_S0, function(x) {
  x <- x[is.finite(x)]
  x <- x[x > 0]
  pmin(x, ar_linear_cap)
})

draws_list_prob_S0 <- split(AR_draw_by_region$AR_S0, AR_draw_by_region$Region)
draws_list_prob_S0 <- lapply(draws_list_prob_S0, function(x) {
  x <- x[is.finite(x)]
  x <- x[x > 0]
  pmin(pmax(x, eps), 1 - eps)
})

## Backward-compatible default for downstream FOI scaling
## Primary pathway now defaults to linear I/S0 (capped only for probability interfaces).
draws_list <- if (use_linear_ar_primary) draws_list_lin_S0 else draws_list_prob_S0

## Region summaries
AR_summary_by_region <- AR_draw_by_region %>%
  dplyr::filter(is.finite(.data$AR_S0), .data$AR_S0 > 0, .data$AR_S0 < 1) %>%
  dplyr::group_by(.data$Region) %>%
  dplyr::summarise(
    AR_lin_mid_S0 = stats::median(.data$AR_lin_S0, na.rm = TRUE),
    AR_lin_lo_S0 = stats::quantile(.data$AR_lin_S0, 0.025, na.rm = TRUE, names = FALSE),
    AR_lin_hi_S0 = stats::quantile(.data$AR_lin_S0, 0.975, na.rm = TRUE, names = FALSE),
    AR_mid_S0 = stats::median(.data$AR_S0, na.rm = TRUE),
    AR_lo_S0 = stats::quantile(.data$AR_S0, 0.025, na.rm = TRUE, names = FALSE),
    AR_hi_S0 = stats::quantile(.data$AR_S0, 0.975, na.rm = TRUE, names = FALSE),
    frac_inf_gt_S0 = mean(.data$flag_inf_gt_S0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(region = Region) %>%
  dplyr::mutate(
    category_med_AR_lin_S0 = dplyr::case_when(
      .data$AR_lin_mid_S0 >= 0.10 ~ "High",
      .data$AR_lin_mid_S0 > 0.05 ~ "Moderate",
      TRUE ~ "Low"
    ),
    category_med_AR_S0 = dplyr::case_when(
      .data$AR_mid_S0 >= 0.10 ~ "High",
      .data$AR_mid_S0 > 0.05 ~ "Moderate",
      TRUE ~ "Low"
    )
  )

## Region-level probability of AR categories (current threshold scheme)
## - high:     AR > 10%
## - moderate: 5% < AR <= 10%
## - low:      AR <= 5%
AR_prob_by_region <- AR_draw_by_region %>%
  dplyr::group_by(.data$Region) %>%
  dplyr::summarise(
    n_draws = dplyr::n(),

    # Bounded AR definition (AR_S0 = 1-exp(-I/S0))
    p_high_AR_S0 = mean(.data$AR_S0 > 0.10, na.rm = TRUE),
    p_moderate_AR_S0 = mean(.data$AR_S0 > 0.05 & .data$AR_S0 <= 0.10, na.rm = TRUE),
    p_low_AR_S0 = mean(.data$AR_S0 <= 0.05, na.rm = TRUE),
    p_ge_5_AR_S0 = mean(.data$AR_S0 > 0.05, na.rm = TRUE),
    p_ge_10_AR_S0 = mean(.data$AR_S0 > 0.10, na.rm = TRUE),

    # Linear AR definition (AR_lin_S0 = I/S0)
    p_high_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.10, na.rm = TRUE),
    p_moderate_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.05 & .data$AR_lin_S0 <= 0.10, na.rm = TRUE),
    p_low_AR_lin_S0 = mean(.data$AR_lin_S0 <= 0.05, na.rm = TRUE),
    p_ge_5_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.05, na.rm = TRUE),
    p_ge_10_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.10, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  dplyr::rename(region = Region) %>%
  dplyr::left_join(
    AR_summary_by_region %>%
      dplyr::select(region, AR_mid_S0, AR_lin_mid_S0),
    by = "region"
  ) %>%
  dplyr::mutate(
    p_high_AR_S0_pct = 100 * .data$p_high_AR_S0,
    p_moderate_AR_S0_pct = 100 * .data$p_moderate_AR_S0,
    p_low_AR_S0_pct = 100 * .data$p_low_AR_S0,
    p_ge_5_AR_S0_pct = 100 * .data$p_ge_5_AR_S0,
    p_ge_10_AR_S0_pct = 100 * .data$p_ge_10_AR_S0,
    p_high_AR_lin_S0_pct = 100 * .data$p_high_AR_lin_S0,
    p_moderate_AR_lin_S0_pct = 100 * .data$p_moderate_AR_lin_S0,
    p_low_AR_lin_S0_pct = 100 * .data$p_low_AR_lin_S0,
    p_ge_5_AR_lin_S0_pct = 100 * .data$p_ge_5_AR_lin_S0,
    p_ge_10_AR_lin_S0_pct = 100 * .data$p_ge_10_AR_lin_S0,
    p_ge_5_AR_S0_round_pct = round(.data$p_ge_5_AR_S0_pct, digits = prob_round_digits),
    p_ge_10_AR_S0_round_pct = round(.data$p_ge_10_AR_S0_pct, digits = prob_round_digits),
    p_ge_5_AR_lin_S0_round_pct = round(.data$p_ge_5_AR_lin_S0_pct, digits = prob_round_digits),
    p_ge_10_AR_lin_S0_round_pct = round(.data$p_ge_10_AR_lin_S0_pct, digits = prob_round_digits),

    # Temporary majority-rule category; overwritten below by fitted probability thresholds.
    category_AR_S0 = dplyr::case_when(
      .data$p_ge_10_AR_S0_round_pct >= prob_majority_pct ~ "High",
      .data$p_ge_5_AR_S0_round_pct >= prob_majority_pct ~ "Moderate",
      TRUE ~ "Low"
    ),
    category_AR_lin_S0 = dplyr::case_when(
      .data$p_ge_10_AR_lin_S0_round_pct >= prob_majority_pct ~ "High",
      .data$p_ge_5_AR_lin_S0_round_pct >= prob_majority_pct ~ "Moderate",
      TRUE ~ "Low"
    ),

    # Median-based category using fixed cutoffs:
    # low: <=5%, moderate: (5%,10%], high: >10%
    category_med_AR_S0 = dplyr::case_when(
      .data$AR_mid_S0 > 0.10 ~ "High",
      .data$AR_mid_S0 > 0.05 ~ "Moderate",
      TRUE ~ "Low"
    ),
    category_med_AR_lin_S0 = dplyr::case_when(
      .data$AR_lin_mid_S0 > 0.10 ~ "High",
      .data$AR_lin_mid_S0 > 0.05 ~ "Moderate",
      TRUE ~ "Low"
    )
  )

## ------------------------------------------------------------------
## Method B: fixed symptomatic fraction (legacy-compatible)
## ------------------------------------------------------------------
symp_const <- symp_const_fixed

true_inf_by_region_draw_B <- tot_symp_by_region_draw %>%
  dplyr::mutate(tot_inf = .data$tot_symp / symp_const)

S0_by_region_draw_B <- true_inf_by_region_draw_B %>%
  dplyr::distinct(Region, draw_id) %>%
  dplyr::left_join(S0_by_region_draw %>% dplyr::distinct(Region, draw_id, total_S0),
                   by = c("Region", "draw_id"))

AR_draw_by_region_B <- true_inf_by_region_draw_B %>%
  dplyr::left_join(S0_by_region_draw_B, by = c("Region", "draw_id")) %>%
  dplyr::left_join(pop_by_state, by = c("Region" = "region")) %>%
  dplyr::mutate(
    AR_lin_S0 = .data$tot_inf / .data$total_S0,
    AR_S0 = 1 - exp(-.data$AR_lin_S0),
    AR_lin_pop = .data$tot_inf / .data$tot_pop,
    AR_pop = 1 - exp(-.data$AR_lin_pop)
  )

AR_summary_by_region_B <- AR_draw_by_region_B %>%
  dplyr::group_by(.data$Region) %>%
  dplyr::summarise(
    AR_lin_mid_S0 = stats::median(.data$AR_lin_S0, na.rm = TRUE),
    AR_lin_lo_S0 = stats::quantile(.data$AR_lin_S0, 0.025, na.rm = TRUE, names = FALSE),
    AR_lin_hi_S0 = stats::quantile(.data$AR_lin_S0, 0.975, na.rm = TRUE, names = FALSE),
    AR_mid_S0 = stats::median(.data$AR_S0, na.rm = TRUE),
    AR_lo_S0 = stats::quantile(.data$AR_S0, 0.025, na.rm = TRUE, names = FALSE),
    AR_hi_S0 = stats::quantile(.data$AR_S0, 0.975, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  dplyr::rename(region = Region) %>%
  dplyr::mutate(
    category_med_AR_lin_S0 = dplyr::case_when(
      .data$AR_lin_mid_S0 >= 0.10 ~ "High",
      .data$AR_lin_mid_S0 > 0.05 ~ "Moderate",
      TRUE ~ "Low"
    ),
    category_med_AR_S0 = dplyr::case_when(
      .data$AR_mid_S0 >= 0.10 ~ "High",
      .data$AR_mid_S0 > 0.05 ~ "Moderate",
      TRUE ~ "Low"
    )
  )

AR_prob_by_region_B <- AR_draw_by_region_B %>%
  dplyr::group_by(.data$Region) %>%
  dplyr::summarise(
    n_draws = dplyr::n(),
    p_high_AR_S0 = mean(.data$AR_S0 > 0.10, na.rm = TRUE),
    p_moderate_AR_S0 = mean(.data$AR_S0 > 0.05 & .data$AR_S0 <= 0.10, na.rm = TRUE),
    p_low_AR_S0 = mean(.data$AR_S0 <= 0.05, na.rm = TRUE),
    p_ge_5_AR_S0 = mean(.data$AR_S0 > 0.05, na.rm = TRUE),
    p_ge_10_AR_S0 = mean(.data$AR_S0 > 0.10, na.rm = TRUE),
    p_high_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.10, na.rm = TRUE),
    p_moderate_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.05 & .data$AR_lin_S0 <= 0.10, na.rm = TRUE),
    p_low_AR_lin_S0 = mean(.data$AR_lin_S0 <= 0.05, na.rm = TRUE),
    p_ge_5_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.05, na.rm = TRUE),
    p_ge_10_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(region = Region) %>%
  dplyr::mutate(
    p_high_AR_S0_pct = 100 * .data$p_high_AR_S0,
    p_moderate_AR_S0_pct = 100 * .data$p_moderate_AR_S0,
    p_low_AR_S0_pct = 100 * .data$p_low_AR_S0,
    p_ge_5_AR_S0_pct = 100 * .data$p_ge_5_AR_S0,
    p_ge_10_AR_S0_pct = 100 * .data$p_ge_10_AR_S0,
    p_high_AR_lin_S0_pct = 100 * .data$p_high_AR_lin_S0,
    p_moderate_AR_lin_S0_pct = 100 * .data$p_moderate_AR_lin_S0,
    p_low_AR_lin_S0_pct = 100 * .data$p_low_AR_lin_S0,
    p_ge_5_AR_lin_S0_pct = 100 * .data$p_ge_5_AR_lin_S0,
    p_ge_10_AR_lin_S0_pct = 100 * .data$p_ge_10_AR_lin_S0,
    p_ge_5_AR_S0_round_pct = round(.data$p_ge_5_AR_S0_pct, digits = prob_round_digits),
    p_ge_10_AR_S0_round_pct = round(.data$p_ge_10_AR_S0_pct, digits = prob_round_digits),
    p_ge_5_AR_lin_S0_round_pct = round(.data$p_ge_5_AR_lin_S0_pct, digits = prob_round_digits),
    p_ge_10_AR_lin_S0_round_pct = round(.data$p_ge_10_AR_lin_S0_pct, digits = prob_round_digits),
    # Temporary majority-rule category; overwritten below by fitted probability thresholds.
    category_AR_S0 = dplyr::case_when(
      .data$p_ge_10_AR_S0_round_pct >= prob_majority_pct ~ "High",
      .data$p_ge_5_AR_S0_round_pct >= prob_majority_pct ~ "Moderate",
      TRUE ~ "Low"
    ),
    category_AR_lin_S0 = dplyr::case_when(
      .data$p_ge_10_AR_lin_S0_round_pct >= prob_majority_pct ~ "High",
      .data$p_ge_5_AR_lin_S0_round_pct >= prob_majority_pct ~ "Moderate",
      TRUE ~ "Low"
    )
  )

## ------------------------------------------------------------------
## Legacy setting alignment + data-driven cutoff fitting
## ------------------------------------------------------------------
legacy_setting_key <- c(
  "Ceará" = "High",
  "Piauí" = "High",
  "Paraíba" = "High",
  "Paraiba" = "High",
  "Alagoas" = "High",
  "Tocantins" = "Moderate",
  "Pernambuco" = "Moderate",
  "Bahia" = "Low",
  "Rio Grande do Norte" = "Low",
  "Minas Gerais" = "Low",
  "Sergipe" = "Low",
  "Goiás" = "Low"
)

classify_cutoffs <- function(x, cut_low, cut_high) {
  dplyr::case_when(
    x > cut_high ~ "High",
    x > cut_low ~ "Moderate",
    TRUE ~ "Low"
  )
}

classify_prob_thresholds <- function(p_gt10_pct, p_gt5_pct, th_high_pct, th_moderate_pct) {
  dplyr::case_when(
    p_gt10_pct >= th_high_pct ~ "High",
    p_gt5_pct >= th_moderate_pct ~ "Moderate",
    TRUE ~ "Low"
  )
}

fit_prob_thresholds_to_legacy <- function(
    prob_df,
    legacy_df,
    region_col = "region",
    p_gt10_col = "p_ge_10_AR_lin_S0_pct",
    p_gt5_col = "p_ge_5_AR_lin_S0_pct",
    legacy_col = "setting_legacy",
    required_matches = NULL
) {
  df <- prob_df %>%
    dplyr::left_join(
      legacy_df %>% dplyr::select(dplyr::all_of(c(region_col, legacy_col))),
      by = region_col
    )

  p10 <- as.numeric(df[[p_gt10_col]])
  p5 <- as.numeric(df[[p_gt5_col]])
  y <- as.character(df[[legacy_col]])
  reg <- as.character(df[[region_col]])

  th_seq <- seq(0, 100, by = 1)
  grid <- expand.grid(th_high_pct = th_seq, th_moderate_pct = th_seq)

  required_idx <- integer(0)
  required_target <- character(0)
  if (!is.null(required_matches) && length(required_matches) > 0) {
    req_regions <- names(required_matches)
    required_idx <- match(req_regions, reg)
    keep <- !is.na(required_idx)
    required_idx <- required_idx[keep]
    required_target <- as.character(required_matches[keep])
  }

  eval_grid <- function(r) {
    pred <- classify_prob_thresholds(
      p_gt10_pct = p10,
      p_gt5_pct = p5,
      th_high_pct = r[["th_high_pct"]],
      th_moderate_pct = r[["th_moderate_pct"]]
    )
    acc <- mean(pred == y, na.rm = TRUE)
    req_ok <- TRUE
    if (length(required_idx) > 0) {
      req_ok <- all(pred[required_idx] == required_target)
    }
    list(acc = acc, req_ok = req_ok)
  }

  eval_res <- apply(grid, 1, eval_grid)
  acc <- vapply(eval_res, function(z) z$acc, numeric(1))
  req_ok <- vapply(eval_res, function(z) z$req_ok, logical(1))

  if (any(req_ok)) {
    candidate <- which(req_ok)
    best <- candidate[which.max(acc[candidate])]
    return(list(
      th_high_pct = grid$th_high_pct[best],
      th_moderate_pct = grid$th_moderate_pct[best],
      accuracy = acc[best],
      constrained = TRUE
    ))
  }

  best <- which.max(acc)
  if (length(required_idx) > 0) {
    warning("No probability threshold pair satisfies required region classes; using unconstrained best-accuracy thresholds.")
  }
  list(
    th_high_pct = grid$th_high_pct[best],
    th_moderate_pct = grid$th_moderate_pct[best],
    accuracy = acc[best],
    constrained = FALSE
  )
}

fit_cutoffs_to_legacy <- function(
    df,
    value_col = "AR_lin_mid_S0",
    legacy_col = "setting_legacy",
    region_col = "region",
    required_matches = NULL
) {
  x <- df[[value_col]]
  y <- as.character(df[[legacy_col]])
  reg <- as.character(df[[region_col]])

  mids <- sort(unique((sort(unique(x))[-1] + sort(unique(x))[-length(unique(x))]) / 2))
  if (length(mids) == 0) {
    return(list(cut_low = NA_real_, cut_high = NA_real_, accuracy = NA_real_, constrained = NA))
  }

  grid <- expand.grid(cut_low = mids, cut_high = mids)
  grid <- grid[grid$cut_low < grid$cut_high, , drop = FALSE]

  if (nrow(grid) == 0) {
    return(list(cut_low = NA_real_, cut_high = NA_real_, accuracy = NA_real_, constrained = NA))
  }

  required_idx <- integer(0)
  required_target <- character(0)
  if (!is.null(required_matches) && length(required_matches) > 0) {
    req_regions <- names(required_matches)
    required_idx <- match(req_regions, reg)
    keep <- !is.na(required_idx)
    required_idx <- required_idx[keep]
    required_target <- as.character(required_matches[keep])
  }

  eval_grid <- function(r) {
    pred <- classify_cutoffs(x, r[["cut_low"]], r[["cut_high"]])
    acc <- mean(pred == y, na.rm = TRUE)
    req_ok <- TRUE
    if (length(required_idx) > 0) {
      req_ok <- all(pred[required_idx] == required_target)
    }
    list(acc = acc, req_ok = req_ok)
  }

  eval_res <- apply(grid, 1, eval_grid)
  acc <- vapply(eval_res, function(z) z$acc, numeric(1))
  req_ok <- vapply(eval_res, function(z) z$req_ok, logical(1))

  if (any(req_ok)) {
    candidate <- which(req_ok)
    best <- candidate[which.max(acc[candidate])]
    return(list(
      cut_low = grid$cut_low[best],
      cut_high = grid$cut_high[best],
      accuracy = acc[best],
      constrained = TRUE
    ))
  }

  best <- which.max(acc)
  if (length(required_idx) > 0) {
    warning("No cutoff pair satisfies required region classes; using unconstrained best-accuracy cutoffs.")
  }
  list(
    cut_low = grid$cut_low[best],
    cut_high = grid$cut_high[best],
    accuracy = acc[best],
    constrained = FALSE
  )
}

AR_summary_by_region <- AR_summary_by_region %>%
  dplyr::mutate(setting_legacy = unname(legacy_setting_key[region]))

AR_summary_by_region_B <- AR_summary_by_region_B %>%
  dplyr::mutate(setting_legacy = unname(legacy_setting_key[region]))

# Fit cutoffs to legacy labels using median linear AR (Method A / Method B)
required_legacy_classes <- c(
  "Pernambuco" = "Moderate",
  "Paraíba" = "High",
  "Paraiba" = "High",
  "Tocantins" = "Moderate"
)

legacy_fit_A <- fit_cutoffs_to_legacy(
  AR_summary_by_region,
  value_col = "AR_lin_mid_S0",
  legacy_col = "setting_legacy",
  region_col = "region",
  required_matches = required_legacy_classes
)
legacy_fit_B <- fit_cutoffs_to_legacy(
  AR_summary_by_region_B,
  value_col = "AR_lin_mid_S0",
  legacy_col = "setting_legacy",
  region_col = "region",
  required_matches = required_legacy_classes
)

# Fit probability thresholds (Pr(AR>10%), Pr(AR>5%)) to legacy labels.
prob_fit_A <- fit_prob_thresholds_to_legacy(
  prob_df = AR_prob_by_region,
  legacy_df = AR_summary_by_region,
  region_col = "region",
  p_gt10_col = "p_ge_10_AR_lin_S0_pct",
  p_gt5_col = "p_ge_5_AR_lin_S0_pct",
  legacy_col = "setting_legacy",
  required_matches = required_legacy_classes
)
prob_fit_B <- fit_prob_thresholds_to_legacy(
  prob_df = AR_prob_by_region_B,
  legacy_df = AR_summary_by_region_B,
  region_col = "region",
  p_gt10_col = "p_ge_10_AR_lin_S0_pct",
  p_gt5_col = "p_ge_5_AR_lin_S0_pct",
  legacy_col = "setting_legacy",
  required_matches = required_legacy_classes
)

AR_summary_by_region <- AR_summary_by_region %>%
  dplyr::mutate(
    category_fit_AR_lin_S0 = classify_cutoffs(.data$AR_lin_mid_S0, legacy_fit_A$cut_low, legacy_fit_A$cut_high),
    fit_accuracy_AR_lin_S0 = mean(.data$category_fit_AR_lin_S0 == .data$setting_legacy, na.rm = TRUE)
  )

AR_summary_by_region_B <- AR_summary_by_region_B %>%
  dplyr::mutate(
    category_fit_AR_lin_S0 = classify_cutoffs(.data$AR_lin_mid_S0, legacy_fit_B$cut_low, legacy_fit_B$cut_high),
    fit_accuracy_AR_lin_S0 = mean(.data$category_fit_AR_lin_S0 == .data$setting_legacy, na.rm = TRUE)
  )

AR_prob_by_region <- AR_prob_by_region %>%
  dplyr::mutate(
    category_AR_S0 = classify_prob_thresholds(
      p_gt10_pct = .data$p_ge_10_AR_S0_pct,
      p_gt5_pct = .data$p_ge_5_AR_S0_pct,
      th_high_pct = prob_fit_A$th_high_pct,
      th_moderate_pct = prob_fit_A$th_moderate_pct
    ),
    category_AR_lin_S0 = classify_prob_thresholds(
      p_gt10_pct = .data$p_ge_10_AR_lin_S0_pct,
      p_gt5_pct = .data$p_ge_5_AR_lin_S0_pct,
      th_high_pct = prob_fit_A$th_high_pct,
      th_moderate_pct = prob_fit_A$th_moderate_pct
    ),
    prob_threshold_high_pct = prob_fit_A$th_high_pct,
    prob_threshold_moderate_pct = prob_fit_A$th_moderate_pct,
    prob_fit_accuracy = prob_fit_A$accuracy
  )

AR_prob_by_region_B <- AR_prob_by_region_B %>%
  dplyr::mutate(
    category_AR_S0 = classify_prob_thresholds(
      p_gt10_pct = .data$p_ge_10_AR_S0_pct,
      p_gt5_pct = .data$p_ge_5_AR_S0_pct,
      th_high_pct = prob_fit_B$th_high_pct,
      th_moderate_pct = prob_fit_B$th_moderate_pct
    ),
    category_AR_lin_S0 = classify_prob_thresholds(
      p_gt10_pct = .data$p_ge_10_AR_lin_S0_pct,
      p_gt5_pct = .data$p_ge_5_AR_lin_S0_pct,
      th_high_pct = prob_fit_B$th_high_pct,
      th_moderate_pct = prob_fit_B$th_moderate_pct
    ),
    prob_threshold_high_pct = prob_fit_B$th_high_pct,
    prob_threshold_moderate_pct = prob_fit_B$th_moderate_pct,
    prob_fit_accuracy = prob_fit_B$accuracy
  )

# Optional compact table to compare legacy vs fitted categories
setting_compare_A <- AR_summary_by_region %>%
  dplyr::left_join(
    AR_prob_by_region %>%
      dplyr::select(region, category_AR_lin_S0, p_ge_5_AR_lin_S0_pct, p_ge_10_AR_lin_S0_pct, p_low_AR_lin_S0_pct),
    by = "region"
  ) %>%
  dplyr::mutate(
    setting_threshold_rule = dplyr::case_when(
      .data$setting_legacy == "High" ~ "AR > 10%",
      .data$setting_legacy == "Moderate" ~ "AR > 5%",
      TRUE ~ "AR <= 5%"
    ),
    p_meet_setting_AR_lin_S0_pct = dplyr::case_when(
      .data$setting_legacy == "High" ~ .data$p_ge_10_AR_lin_S0_pct,
      .data$setting_legacy == "Moderate" ~ .data$p_ge_5_AR_lin_S0_pct,
      TRUE ~ .data$p_low_AR_lin_S0_pct
    )
  ) %>%
  dplyr::rename(
    category_prob_AR_lin_S0 = category_AR_lin_S0,
    pr_ar_gt_5_pct = p_ge_5_AR_lin_S0_pct,
    pr_ar_gt_10_pct = p_ge_10_AR_lin_S0_pct,
    pr_meet_setting_pct = p_meet_setting_AR_lin_S0_pct
  ) %>%
  dplyr::mutate(
    category_final_AR_lin_S0 = .data$category_prob_AR_lin_S0
  ) %>%
  dplyr::select(
    region, category_final_AR_lin_S0, category_prob_AR_lin_S0, pr_ar_gt_10_pct, pr_ar_gt_5_pct,
    setting_legacy, setting_threshold_rule, pr_meet_setting_pct, AR_lin_mid_S0,
    category_med_AR_lin_S0, category_fit_AR_lin_S0, fit_accuracy_AR_lin_S0
  )

setting_compare_B <- AR_summary_by_region_B %>%
  dplyr::left_join(
    AR_prob_by_region_B %>%
      dplyr::select(region, category_AR_lin_S0, p_ge_5_AR_lin_S0_pct, p_ge_10_AR_lin_S0_pct, p_low_AR_lin_S0_pct),
    by = "region"
  ) %>%
  dplyr::mutate(
    setting_threshold_rule = dplyr::case_when(
      .data$setting_legacy == "High" ~ "AR > 10%",
      .data$setting_legacy == "Moderate" ~ "AR > 5%",
      TRUE ~ "AR <= 5%"
    ),
    p_meet_setting_AR_lin_S0_pct = dplyr::case_when(
      .data$setting_legacy == "High" ~ .data$p_ge_10_AR_lin_S0_pct,
      .data$setting_legacy == "Moderate" ~ .data$p_ge_5_AR_lin_S0_pct,
      TRUE ~ .data$p_low_AR_lin_S0_pct
    )
  ) %>%
  dplyr::rename(
    category_prob_AR_lin_S0 = category_AR_lin_S0,
    pr_ar_gt_5_pct = p_ge_5_AR_lin_S0_pct,
    pr_ar_gt_10_pct = p_ge_10_AR_lin_S0_pct,
    pr_meet_setting_pct = p_meet_setting_AR_lin_S0_pct
  ) %>%
  dplyr::mutate(
    category_final_AR_lin_S0 = .data$category_prob_AR_lin_S0
  ) %>%
  dplyr::select(
    region, category_final_AR_lin_S0, category_prob_AR_lin_S0, pr_ar_gt_10_pct, pr_ar_gt_5_pct,
    setting_legacy, setting_threshold_rule, pr_meet_setting_pct, AR_lin_mid_S0,
    category_med_AR_lin_S0, category_fit_AR_lin_S0, fit_accuracy_AR_lin_S0
  )

AR_summary_selected <- if (use_fixed_symp_const) AR_summary_by_region_B else AR_summary_by_region
AR_prob_selected <- if (use_fixed_symp_const) AR_prob_by_region_B else AR_prob_by_region
setting_compare <- if (use_fixed_symp_const) setting_compare_B else setting_compare_A
legacy_fit_selected <- if (use_fixed_symp_const) legacy_fit_B else legacy_fit_A
prob_fit_selected <- if (use_fixed_symp_const) prob_fit_B else prob_fit_A
symp_setting_used <- if (use_fixed_symp_const) {
  paste0("fixed_", symp_const_fixed)
} else {
  "drawwise"
}

# Useful to print/log which fitted cutoffs are being used for classification alignment.
legacy_fit_cutoffs_pct <- c(
  cut_low_pct = 100 * legacy_fit_selected$cut_low,
  cut_high_pct = 100 * legacy_fit_selected$cut_high,
  fit_accuracy = legacy_fit_selected$accuracy
)

prob_fit_thresholds_pct <- c(
  high_threshold_pct = prob_fit_selected$th_high_pct,
  moderate_threshold_pct = prob_fit_selected$th_moderate_pct,
  fit_accuracy = prob_fit_selected$accuracy
)

AR_by_state <- list(
  mid = setNames(pmin(AR_summary_selected$AR_lin_mid_S0, ar_linear_cap), AR_summary_selected$region),
  lo = setNames(pmin(AR_summary_selected$AR_lin_lo_S0, ar_linear_cap), AR_summary_selected$region),
  hi = setNames(pmin(AR_summary_selected$AR_lin_hi_S0, ar_linear_cap), AR_summary_selected$region)
)

AR_by_state_prob <- lapply(AR_by_state, function(v) 1 - exp(-v))

## Optional plot (bounded AR_S0)
p <- ggplot(AR_draw_by_region, aes(x = .data$AR_S0)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30) +
  geom_density(linewidth = 0.8) +
  facet_wrap(~Region, scales = "free_y") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Attack rate (AR_S0 = 1 - exp(-tot_inf / total_S0))",
    y = "Density",
    caption = paste0(
      "Baseline stratum: ", baseline_cov, ", ", baseline_ve, ", Scenario ", baseline_scenario,
      ". Infections inferred as total_pre_true / symp_prop."
    )
  ) +
  theme_bw(base_size = 10)

ggsave("06_Results/suppl_attack_rates_drawlevel.pdf", plot = p, width = 8, height = 6, device = cairo_pdf)

scale_foi_to_target_ar <- function(foi_daily, ar_target) {
  H0 <- sum(foi_daily, na.rm = TRUE)
  Ht <- -log(1 - ar_target)
  m <- Ht / H0
  foi_daily * m
}

foi_daily_adj_by_state <- lapply(names(foi_daily_by_state), function(st) {
  foi0 <- foi_daily_by_state[[st]]
  list(
    mid = scale_foi_to_target_ar(foi0, AR_by_state_prob$mid[[st]]),
    lo = scale_foi_to_target_ar(foi0, AR_by_state_prob$lo[[st]]),
    hi = scale_foi_to_target_ar(foi0, AR_by_state_prob$hi[[st]])
  )
})
names(foi_daily_adj_by_state) <- names(foi_daily_by_state)
foi_daily_by_state_mid <- lapply(foi_daily_adj_by_state, `[[`, "mid")
states_to_run <- names(foi_daily_by_state_mid)
