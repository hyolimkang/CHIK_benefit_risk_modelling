## ======================================================================
## Draw-keyed attack-rate construction (clean integrated version)
## - Primary classification rule updated:
##   classify each state by the attack-rate band with the highest
##   posterior probability among:
##     Low       : AR <= 5%
##     Moderate  : 5% < AR <= 10%
##     High      : AR > 10%
## - No grey zone / uncertainty category
## - Keeps both Method A (draw-wise symptomatic fraction)
##   and Method B (fixed symptomatic fraction)
## ======================================================================

## Assumes setup.R has already been sourced and these objects exist:
## sim_results_vc_ixchiq_model, all_draws_ix_true, lhs_sample_young,
## pop_by_state, foi_daily_by_state

## ----------------------------------------------------------------------
## Global settings
## ----------------------------------------------------------------------
baseline_cov <- "cov50"
baseline_ve <- "VE0"
baseline_scenario <- 1L

use_linear_ar_primary <- TRUE
ar_linear_cap <- 0.99

use_fixed_symp_const <- TRUE
symp_const_fixed <- 0.524

## legacy setting labels (optional alignment / comparison only)
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

required_legacy_classes <- c(
  "Pernambuco" = "Moderate",
  "Paraíba" = "High",
  "Paraiba" = "High",
  "Tocantins" = "Moderate"
)

eps <- 1e-12

## ----------------------------------------------------------------------
## Helper functions
## ----------------------------------------------------------------------
get_total_S0_by_stratum <- function(sim_results, t0 = 1) {
  out <- data.frame()
  
  for (region in names(sim_results)) {
    for (ve in names(sim_results[[region]])) {
      for (vc in names(sim_results[[region]][[ve]])) {
        scenario_list <- sim_results[[region]][[ve]][[vc]]
        
        for (sc in seq_along(scenario_list)) {
          scen <- scenario_list[[sc]]
          
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
    dplyr::arrange(
      .data$Region, .data$draw_id,
      dplyr::desc(.data$is_baseline),
      .data$Coverage, .data$VE, .data$scenario_chr
    ) %>%
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

classify_cutoffs <- function(x, cut_low, cut_high) {
  dplyr::case_when(
    x > cut_high ~ "High",
    x > cut_low ~ "Moderate",
    TRUE ~ "Low"
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
  
  ux <- sort(unique(x))
  if (length(ux) < 2) {
    return(list(cut_low = NA_real_, cut_high = NA_real_, accuracy = NA_real_, constrained = NA))
  }
  
  mids <- sort((ux[-1] + ux[-length(ux)]) / 2)
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

classify_by_max_prob <- function(p_low, p_moderate, p_high) {
  dplyr::case_when(
    p_high >= p_moderate & p_high >= p_low ~ "High",
    p_moderate >= p_high & p_moderate >= p_low ~ "Moderate",
    TRUE ~ "Low"
  )
}

summarise_ar_by_region <- function(ar_draw_df) {
  ar_draw_df %>%
    dplyr::filter(is.finite(.data$AR_S0), .data$AR_S0 > 0, .data$AR_S0 < 1) %>%
    dplyr::group_by(.data$Region) %>%
    dplyr::summarise(
      AR_lin_mid_S0 = stats::median(.data$AR_lin_S0, na.rm = TRUE),
      AR_lin_lo_S0  = stats::quantile(.data$AR_lin_S0, 0.025, na.rm = TRUE, names = FALSE),
      AR_lin_hi_S0  = stats::quantile(.data$AR_lin_S0, 0.975, na.rm = TRUE, names = FALSE),
      
      AR_mid_S0 = stats::median(.data$AR_S0, na.rm = TRUE),
      AR_lo_S0  = stats::quantile(.data$AR_S0, 0.025, na.rm = TRUE, names = FALSE),
      AR_hi_S0  = stats::quantile(.data$AR_S0, 0.975, na.rm = TRUE, names = FALSE),
      
      frac_inf_gt_S0 = mean(.data$flag_inf_gt_S0, na.rm = TRUE),
      
      AR_lin_mid_pop = stats::median(.data$AR_lin_pop, na.rm = TRUE),
      AR_lin_lo_pop  = stats::quantile(.data$AR_lin_pop, 0.025, na.rm = TRUE, names = FALSE),
      AR_lin_hi_pop  = stats::quantile(.data$AR_lin_pop, 0.975, na.rm = TRUE, names = FALSE),
      
      AR_mid_pop = stats::median(.data$AR_pop, na.rm = TRUE),
      AR_lo_pop  = stats::quantile(.data$AR_pop, 0.025, na.rm = TRUE, names = FALSE),
      AR_hi_pop  = stats::quantile(.data$AR_pop, 0.975, na.rm = TRUE, names = FALSE),
      
      .groups = "drop"
    ) %>%
    dplyr::rename(region = Region) %>%
    dplyr::mutate(
      category_med_AR_lin_S0 = dplyr::case_when(
        .data$AR_lin_mid_S0 > 0.10 ~ "High",
        .data$AR_lin_mid_S0 > 0.05 ~ "Moderate",
        TRUE ~ "Low"
      ),
      category_med_AR_S0 = dplyr::case_when(
        .data$AR_mid_S0 > 0.10 ~ "High",
        .data$AR_mid_S0 > 0.05 ~ "Moderate",
        TRUE ~ "Low"
      ),
      category_med_AR_lin_pop = dplyr::case_when(
        .data$AR_lin_mid_pop > 0.10 ~ "High",
        .data$AR_lin_mid_pop > 0.05 ~ "Moderate",
        TRUE ~ "Low"
      ),
      category_med_AR_pop = dplyr::case_when(
        .data$AR_mid_pop > 0.10 ~ "High",
        .data$AR_mid_pop > 0.05 ~ "Moderate",
        TRUE ~ "Low"
      )
    )
}

summarise_prob_by_region <- function(ar_draw_df, ar_summary_df) {
  ar_draw_df %>%
    dplyr::group_by(.data$Region) %>%
    dplyr::summarise(
      n_draws = dplyr::n(),
      
      ## bounded AR based on S0
      p_high_AR_S0 = mean(.data$AR_S0 > 0.10, na.rm = TRUE),
      p_moderate_AR_S0 = mean(.data$AR_S0 > 0.05 & .data$AR_S0 <= 0.10, na.rm = TRUE),
      p_low_AR_S0 = mean(.data$AR_S0 <= 0.05, na.rm = TRUE),
      
      ## linear AR based on S0
      p_high_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.10, na.rm = TRUE),
      p_moderate_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.05 & .data$AR_lin_S0 <= 0.10, na.rm = TRUE),
      p_low_AR_lin_S0 = mean(.data$AR_lin_S0 <= 0.05, na.rm = TRUE),
      
      ## bounded AR based on population
      p_high_AR_pop = mean(.data$AR_pop > 0.10, na.rm = TRUE),
      p_moderate_AR_pop = mean(.data$AR_pop > 0.05 & .data$AR_pop <= 0.10, na.rm = TRUE),
      p_low_AR_pop = mean(.data$AR_pop <= 0.05, na.rm = TRUE),
      
      ## linear AR based on population
      p_high_AR_lin_pop = mean(.data$AR_lin_pop > 0.10, na.rm = TRUE),
      p_moderate_AR_lin_pop = mean(.data$AR_lin_pop > 0.05 & .data$AR_lin_pop <= 0.10, na.rm = TRUE),
      p_low_AR_lin_pop = mean(.data$AR_lin_pop <= 0.05, na.rm = TRUE),
      
      ## optional exceedance summaries retained for diagnostics
      p_ge_5_AR_S0 = mean(.data$AR_S0 > 0.05, na.rm = TRUE),
      p_ge_10_AR_S0 = mean(.data$AR_S0 > 0.10, na.rm = TRUE),
      p_ge_5_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.05, na.rm = TRUE),
      p_ge_10_AR_lin_S0 = mean(.data$AR_lin_S0 > 0.10, na.rm = TRUE),
      p_ge_5_AR_pop = mean(.data$AR_pop > 0.05, na.rm = TRUE),
      p_ge_10_AR_pop = mean(.data$AR_pop > 0.10, na.rm = TRUE),
      p_ge_5_AR_lin_pop = mean(.data$AR_lin_pop > 0.05, na.rm = TRUE),
      p_ge_10_AR_lin_pop = mean(.data$AR_lin_pop > 0.10, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    dplyr::rename(region = Region) %>%
    dplyr::left_join(
      ar_summary_df %>%
        dplyr::select(
          region,
          AR_mid_S0, AR_lin_mid_S0,
          AR_mid_pop, AR_lin_mid_pop
        ),
      by = "region"
    ) %>%
    dplyr::mutate(
      p_high_AR_S0_pct = 100 * .data$p_high_AR_S0,
      p_moderate_AR_S0_pct = 100 * .data$p_moderate_AR_S0,
      p_low_AR_S0_pct = 100 * .data$p_low_AR_S0,
      
      p_high_AR_lin_S0_pct = 100 * .data$p_high_AR_lin_S0,
      p_moderate_AR_lin_S0_pct = 100 * .data$p_moderate_AR_lin_S0,
      p_low_AR_lin_S0_pct = 100 * .data$p_low_AR_lin_S0,
      
      p_high_AR_pop_pct = 100 * .data$p_high_AR_pop,
      p_moderate_AR_pop_pct = 100 * .data$p_moderate_AR_pop,
      p_low_AR_pop_pct = 100 * .data$p_low_AR_pop,
      
      p_high_AR_lin_pop_pct = 100 * .data$p_high_AR_lin_pop,
      p_moderate_AR_lin_pop_pct = 100 * .data$p_moderate_AR_lin_pop,
      p_low_AR_lin_pop_pct = 100 * .data$p_low_AR_lin_pop,
      
      p_ge_5_AR_S0_pct = 100 * .data$p_ge_5_AR_S0,
      p_ge_10_AR_S0_pct = 100 * .data$p_ge_10_AR_S0,
      p_ge_5_AR_lin_S0_pct = 100 * .data$p_ge_5_AR_lin_S0,
      p_ge_10_AR_lin_S0_pct = 100 * .data$p_ge_10_AR_lin_S0,
      p_ge_5_AR_pop_pct = 100 * .data$p_ge_5_AR_pop,
      p_ge_10_AR_pop_pct = 100 * .data$p_ge_10_AR_pop,
      p_ge_5_AR_lin_pop_pct = 100 * .data$p_ge_5_AR_lin_pop,
      p_ge_10_AR_lin_pop_pct = 100 * .data$p_ge_10_AR_lin_pop,
      
      ## NEW primary probability classification:
      category_AR_S0 = classify_by_max_prob(
        .data$p_low_AR_S0,
        .data$p_moderate_AR_S0,
        .data$p_high_AR_S0
      ),
      category_AR_lin_S0 = classify_by_max_prob(
        .data$p_low_AR_lin_S0,
        .data$p_moderate_AR_lin_S0,
        .data$p_high_AR_lin_S0
      ),
      category_AR_pop = classify_by_max_prob(
        .data$p_low_AR_pop,
        .data$p_moderate_AR_pop,
        .data$p_high_AR_pop
      ),
      category_AR_lin_pop = classify_by_max_prob(
        .data$p_low_AR_lin_pop,
        .data$p_moderate_AR_lin_pop,
        .data$p_high_AR_lin_pop
      ),
      
      ## median-based reference categories
      category_med_AR_S0 = dplyr::case_when(
        .data$AR_mid_S0 > 0.10 ~ "High",
        .data$AR_mid_S0 > 0.05 ~ "Moderate",
        TRUE ~ "Low"
      ),
      category_med_AR_lin_S0 = dplyr::case_when(
        .data$AR_lin_mid_S0 > 0.10 ~ "High",
        .data$AR_lin_mid_S0 > 0.05 ~ "Moderate",
        TRUE ~ "Low"
      ),
      category_med_AR_pop = dplyr::case_when(
        .data$AR_mid_pop > 0.10 ~ "High",
        .data$AR_mid_pop > 0.05 ~ "Moderate",
        TRUE ~ "Low"
      ),
      category_med_AR_lin_pop = dplyr::case_when(
        .data$AR_lin_mid_pop > 0.10 ~ "High",
        .data$AR_lin_mid_pop > 0.05 ~ "Moderate",
        TRUE ~ "Low"
      )
    )
}

build_ar_draws <- function(true_inf_by_region_draw, S0_by_region_draw, pop_by_state) {
  true_inf_by_region_draw %>%
    dplyr::left_join(S0_by_region_draw, by = c("Region", "draw_id")) %>%
    dplyr::left_join(pop_by_state, by = c("Region" = "region")) %>%
    dplyr::mutate(
      AR_lin_S0 = .data$tot_inf / .data$total_S0,
      AR_lin_S0_pct = 100 * .data$AR_lin_S0,
      flag_inf_gt_S0 = .data$tot_inf > .data$total_S0,
      
      H_S0 = .data$AR_lin_S0,
      AR_S0 = 1 - exp(-.data$H_S0),
      AR_S0_pct = 100 * .data$AR_S0,
      gap_H_minus_AR_S0 = .data$H_S0 - .data$AR_S0,
      
      AR_lin_pop = .data$tot_inf / .data$tot_pop,
      AR_lin_pop_pct = 100 * .data$AR_lin_pop,
      AR_pop = 1 - exp(-.data$AR_lin_pop),
      AR_pop_pct = 100 * .data$AR_pop
    )
}

scale_foi_to_target_ar <- function(foi_daily, ar_target) {
  H0 <- sum(foi_daily, na.rm = TRUE)
  Ht <- -log(1 - ar_target)
  m <- Ht / H0
  foi_daily * m
}

## ----------------------------------------------------------------------
## 1) Build baseline S0 lookup
## ----------------------------------------------------------------------
S0_strata <- get_total_S0_by_stratum(sim_results_vc_ixchiq_model, t0 = 1) %>%
  dplyr::filter(
    VE == baseline_ve,
    VC == baseline_cov,
    Scenario == baseline_scenario
  ) %>%
  dplyr::select(region, draw_id, total_S0)

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

## ----------------------------------------------------------------------
## 2) Build draw-keyed symptomatic burden
## ----------------------------------------------------------------------
drawkey_pre <- build_drawkey_preburden(
  all_draws_df = all_draws_ix_true,
  baseline_cov = baseline_cov,
  baseline_ve = baseline_ve,
  baseline_scenario = baseline_scenario
)

tot_symp_by_region_draw <- drawkey_pre %>%
  dplyr::transmute(Region, draw_id, tot_symp = .data$total_pre_true_key)

## ----------------------------------------------------------------------
## 3) Method A: draw-wise symptomatic fraction
## ----------------------------------------------------------------------
symp_by_draw <- lhs_sample_young %>%
  dplyr::transmute(
    draw_id = dplyr::row_number(),
    symp_prop = as.numeric(symp_overall)
  ) %>%
  dplyr::mutate(symp_prop = pmin(pmax(symp_prop, 1e-6), 1 - 1e-6))

true_inf_by_region_draw_A <- tot_symp_by_region_draw %>%
  dplyr::left_join(symp_by_draw, by = "draw_id") %>%
  dplyr::mutate(tot_inf = .data$tot_symp / .data$symp_prop)

S0_by_region_draw_A <- true_inf_by_region_draw_A %>%
  dplyr::distinct(Region, draw_id) %>%
  dplyr::left_join(S0_lookup, by = c("Region" = "region", "draw_id" = "draw_id"))

if (any(is.na(S0_by_region_draw_A$total_S0))) {
  S0_region_fallback <- S0_lookup %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(total_S0 = dplyr::first(total_S0), .groups = "drop")
  
  S0_by_region_draw_A <- S0_by_region_draw_A %>%
    dplyr::select(-total_S0) %>%
    dplyr::left_join(S0_region_fallback, by = c("Region" = "region"))
}

if (any(is.na(S0_by_region_draw_A$total_S0))) {
  stop("Missing baseline total_S0 for some Region/draw_id after stratum filter (Method A).")
}

AR_draw_by_region_A <- build_ar_draws(
  true_inf_by_region_draw = true_inf_by_region_draw_A,
  S0_by_region_draw = S0_by_region_draw_A,
  pop_by_state = pop_by_state
)

AR_summary_by_region_A <- summarise_ar_by_region(AR_draw_by_region_A) %>%
  dplyr::mutate(setting_legacy = unname(legacy_setting_key[region]))

AR_prob_by_region_A <- summarise_prob_by_region(AR_draw_by_region_A, AR_summary_by_region_A)

legacy_fit_A <- fit_cutoffs_to_legacy(
  AR_summary_by_region_A,
  value_col = "AR_lin_mid_S0",
  legacy_col = "setting_legacy",
  region_col = "region",
  required_matches = required_legacy_classes
)

AR_summary_by_region_A <- AR_summary_by_region_A %>%
  dplyr::mutate(
    category_fit_AR_lin_S0 = classify_cutoffs(.data$AR_lin_mid_S0, legacy_fit_A$cut_low, legacy_fit_A$cut_high),
    fit_accuracy_AR_lin_S0 = mean(.data$category_fit_AR_lin_S0 == .data$setting_legacy, na.rm = TRUE)
  )

setting_compare_A <- AR_summary_by_region_A %>%
  dplyr::left_join(
    AR_prob_by_region_A %>%
      dplyr::select(
        region,
        category_AR_lin_S0,
        p_low_AR_lin_S0_pct,
        p_moderate_AR_lin_S0_pct,
        p_high_AR_lin_S0_pct,
        p_ge_5_AR_lin_S0_pct,
        p_ge_10_AR_lin_S0_pct
      ),
    by = "region"
  ) %>%
  dplyr::rename(
    category_maxprob_AR_lin_S0 = category_AR_lin_S0,
    pr_low_pct = p_low_AR_lin_S0_pct,
    pr_moderate_pct = p_moderate_AR_lin_S0_pct,
    pr_high_pct = p_high_AR_lin_S0_pct,
    pr_ar_gt_5_pct = p_ge_5_AR_lin_S0_pct,
    pr_ar_gt_10_pct = p_ge_10_AR_lin_S0_pct
  ) %>%
  dplyr::mutate(
    setting_threshold_rule = dplyr::case_when(
      .data$setting_legacy == "High" ~ "AR > 10%",
      .data$setting_legacy == "Moderate" ~ "5% < AR <= 10%",
      TRUE ~ "AR <= 5%"
    ),
    pr_meet_setting_pct = dplyr::case_when(
      .data$setting_legacy == "High" ~ .data$pr_high_pct,
      .data$setting_legacy == "Moderate" ~ .data$pr_moderate_pct,
      TRUE ~ .data$pr_low_pct
    ),
    category_final_AR_lin_S0 = .data$category_maxprob_AR_lin_S0
  ) %>%
  dplyr::select(
    region,
    category_final_AR_lin_S0,
    category_maxprob_AR_lin_S0,
    pr_low_pct, pr_moderate_pct, pr_high_pct,
    pr_ar_gt_5_pct, pr_ar_gt_10_pct,
    setting_legacy, setting_threshold_rule, pr_meet_setting_pct,
    AR_lin_mid_S0, category_med_AR_lin_S0,
    category_fit_AR_lin_S0, fit_accuracy_AR_lin_S0
  )

## ----------------------------------------------------------------------
## 4) Method B: fixed symptomatic fraction
## ----------------------------------------------------------------------
symp_const <- symp_const_fixed

true_inf_by_region_draw_B <- tot_symp_by_region_draw %>%
  dplyr::mutate(tot_inf = .data$tot_symp / symp_const)

S0_by_region_draw_B <- true_inf_by_region_draw_B %>%
  dplyr::distinct(Region, draw_id) %>%
  dplyr::left_join(
    S0_by_region_draw_A %>% dplyr::distinct(Region, draw_id, total_S0),
    by = c("Region", "draw_id")
  )

if (any(is.na(S0_by_region_draw_B$total_S0))) {
  stop("Missing baseline total_S0 for some Region/draw_id (Method B).")
}

AR_draw_by_region_B <- build_ar_draws(
  true_inf_by_region_draw = true_inf_by_region_draw_B,
  S0_by_region_draw = S0_by_region_draw_B,
  pop_by_state = pop_by_state
)

AR_summary_by_region_B <- summarise_ar_by_region(AR_draw_by_region_B) %>%
  dplyr::mutate(setting_legacy = unname(legacy_setting_key[region]))

AR_prob_by_region_B <- summarise_prob_by_region(AR_draw_by_region_B, AR_summary_by_region_B)

legacy_fit_B <- fit_cutoffs_to_legacy(
  AR_summary_by_region_B,
  value_col = "AR_lin_mid_S0",
  legacy_col = "setting_legacy",
  region_col = "region",
  required_matches = required_legacy_classes
)

AR_summary_by_region_B <- AR_summary_by_region_B %>%
  dplyr::mutate(
    category_fit_AR_lin_S0 = classify_cutoffs(.data$AR_lin_mid_S0, legacy_fit_B$cut_low, legacy_fit_B$cut_high),
    fit_accuracy_AR_lin_S0 = mean(.data$category_fit_AR_lin_S0 == .data$setting_legacy, na.rm = TRUE)
  )

setting_compare_B <- AR_summary_by_region_B %>%
  dplyr::left_join(
    AR_prob_by_region_B %>%
      dplyr::select(
        region,
        category_AR_lin_S0,
        p_low_AR_lin_S0_pct,
        p_moderate_AR_lin_S0_pct,
        p_high_AR_lin_S0_pct,
        p_ge_5_AR_lin_S0_pct,
        p_ge_10_AR_lin_S0_pct
      ),
    by = "region"
  ) %>%
  dplyr::rename(
    category_maxprob_AR_lin_S0 = category_AR_lin_S0,
    pr_low_pct = p_low_AR_lin_S0_pct,
    pr_moderate_pct = p_moderate_AR_lin_S0_pct,
    pr_high_pct = p_high_AR_lin_S0_pct,
    pr_ar_gt_5_pct = p_ge_5_AR_lin_S0_pct,
    pr_ar_gt_10_pct = p_ge_10_AR_lin_S0_pct
  ) %>%
  dplyr::mutate(
    setting_threshold_rule = dplyr::case_when(
      .data$setting_legacy == "High" ~ "AR > 10%",
      .data$setting_legacy == "Moderate" ~ "5% < AR <= 10%",
      TRUE ~ "AR <= 5%"
    ),
    pr_meet_setting_pct = dplyr::case_when(
      .data$setting_legacy == "High" ~ .data$pr_high_pct,
      .data$setting_legacy == "Moderate" ~ .data$pr_moderate_pct,
      TRUE ~ .data$pr_low_pct
    ),
    category_final_AR_lin_S0 = .data$category_maxprob_AR_lin_S0
  ) %>%
  dplyr::select(
    region,
    category_final_AR_lin_S0,
    category_maxprob_AR_lin_S0,
    pr_low_pct, pr_moderate_pct, pr_high_pct,
    pr_ar_gt_5_pct, pr_ar_gt_10_pct,
    setting_legacy, setting_threshold_rule, pr_meet_setting_pct,
    AR_lin_mid_S0, category_med_AR_lin_S0,
    category_fit_AR_lin_S0, fit_accuracy_AR_lin_S0
  )

## ----------------------------------------------------------------------
## 5) Select primary method outputs
## ----------------------------------------------------------------------
AR_draw_selected <- if (use_fixed_symp_const) AR_draw_by_region_B else AR_draw_by_region_A
AR_summary_selected <- if (use_fixed_symp_const) AR_summary_by_region_B else AR_summary_by_region_A
AR_prob_selected <- if (use_fixed_symp_const) AR_prob_by_region_B else AR_prob_by_region_A
setting_compare <- if (use_fixed_symp_const) setting_compare_B else setting_compare_A
legacy_fit_selected <- if (use_fixed_symp_const) legacy_fit_B else legacy_fit_A
symp_setting_used <- if (use_fixed_symp_const) paste0("fixed_", symp_const_fixed) else "drawwise"

legacy_fit_cutoffs_pct <- c(
  cut_low_pct = 100 * legacy_fit_selected$cut_low,
  cut_high_pct = 100 * legacy_fit_selected$cut_high,
  fit_accuracy = legacy_fit_selected$accuracy
)

## ----------------------------------------------------------------------
## 6) Draw lists for downstream FOI scaling
## ----------------------------------------------------------------------
draws_list_lin_S0 <- split(AR_draw_selected$AR_lin_S0, AR_draw_selected$Region)
draws_list_lin_S0 <- lapply(draws_list_lin_S0, function(x) {
  x <- x[is.finite(x)]
  x <- x[x > 0]
  pmin(x, ar_linear_cap)
})

draws_list_prob_S0 <- split(AR_draw_selected$AR_S0, AR_draw_selected$Region)
draws_list_prob_S0 <- lapply(draws_list_prob_S0, function(x) {
  x <- x[is.finite(x)]
  x <- x[x > 0]
  pmin(pmax(x, eps), 1 - eps)
})

draws_list <- if (use_linear_ar_primary) draws_list_lin_S0 else draws_list_prob_S0

## ----------------------------------------------------------------------
## 7) AR_by_state summaries for downstream use
## ----------------------------------------------------------------------
AR_by_state <- list(
  mid = setNames(pmin(AR_summary_selected$AR_lin_mid_S0, ar_linear_cap), AR_summary_selected$region),
  lo  = setNames(pmin(AR_summary_selected$AR_lin_lo_S0,  ar_linear_cap), AR_summary_selected$region),
  hi  = setNames(pmin(AR_summary_selected$AR_lin_hi_S0,  ar_linear_cap), AR_summary_selected$region)
)

AR_by_state_prob <- lapply(AR_by_state, function(v) 1 - exp(-v))

## ----------------------------------------------------------------------
## 8) Optional plot
## ----------------------------------------------------------------------
p <- ggplot(AR_draw_selected, aes(x = .data$AR_S0)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30) +
  geom_density(linewidth = 0.8) +
  facet_wrap(~Region, scales = "free_y") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Attack rate (AR_S0 = 1 - exp(-tot_inf / total_S0))",
    y = "Density",
    caption = paste0(
      "Baseline stratum: ", baseline_cov, ", ", baseline_ve, ", Scenario ", baseline_scenario,
      ". Infections inferred as total_pre_true / symptomatic fraction.",
      " Classification uses maximum posterior probability over Low/Moderate/High bands."
    )
  ) +
  theme_bw(base_size = 10)

ggsave(
  "06_Results/suppl_attack_rates_drawlevel.pdf",
  plot = p,
  width = 8,
  height = 6,
  device = cairo_pdf
)

## ----------------------------------------------------------------------
## 9) FOI scaling
## ----------------------------------------------------------------------
foi_daily_adj_by_state <- lapply(names(foi_daily_by_state), function(st) {
  foi0 <- foi_daily_by_state[[st]]
  list(
    mid = scale_foi_to_target_ar(foi0, AR_by_state_prob$mid[[st]]),
    lo  = scale_foi_to_target_ar(foi0, AR_by_state_prob$lo[[st]]),
    hi  = scale_foi_to_target_ar(foi0, AR_by_state_prob$hi[[st]])
  )
})
names(foi_daily_adj_by_state) <- names(foi_daily_by_state)

foi_daily_by_state_mid <- lapply(foi_daily_adj_by_state, `[[`, "mid")
states_to_run <- names(foi_daily_by_state_mid)

## ----------------------------------------------------------------------
## 10) Useful outputs to inspect
## ----------------------------------------------------------------------
print(symp_setting_used)
print(legacy_fit_cutoffs_pct)

print(
  AR_prob_selected %>%
    dplyr::select(
      region,
      category_AR_lin_S0,
      p_low_AR_lin_S0_pct,
      p_moderate_AR_lin_S0_pct,
      p_high_AR_lin_S0_pct,
      AR_lin_mid_S0
    ) %>%
    dplyr::arrange(dplyr::desc(p_high_AR_lin_S0_pct), dplyr::desc(p_moderate_AR_lin_S0_pct))
)

print(setting_compare)