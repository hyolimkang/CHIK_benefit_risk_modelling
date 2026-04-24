# =============================================================================
# brazil_all_draws_ori_v2_cov_alt.R
#
# Replay Section 8 (BRR pipeline and reporting tables) from
# brazil_all_draws_ori_v2.R for *alternative* coverage levels
# (cov10 and cov90), using the exact same pipeline structure as cov50.
#
# The upstream objects that are already coverage-agnostic
# (q_all_regions, risk_components_all, risk_join_all, risk_draw_df,
#  daly_pars_true, all_draws_*_true, combined_nnv_df_region_coverage_model,
#  setting_key, map_scenario_agecat_int, age_map, rr_vals, etc.) are
# reused directly. Only the coverage-dependent pieces
# (tot_vacc_map, benefit_draws, benefit_base, benefit_draw_df,
#  draw_level_xy_serostatus, and everything downstream) are rebuilt
# per coverage.
#
# Prereq: run brazil_all_draws_ori_v2.R first so the following
# objects exist in the global workspace:
#   - all_draws_daly_true, all_draws_sae_true, all_draws_fatal_true
#   - combined_nnv_df_region_coverage_model
#   - postsim_vc_ixchiq_model
#   - lhs_sample, daly_pars_true, risk_draw_df
#   - q_all_regions, risk_components_all, risk_join_all
#   - setting_key, age_map, rr_vals, map_scenario_agecat_int
#   - helpers: compute_daly_one, make_averted_draws_true,
#              fmt_ci, make_pr_gt1_wide, make_brr_ceac_outbreak,
#              plot_brr_ceac_outbreak_ve
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(flextable)
library(officer)

# -----------------------------------------------------------------------------
# Ensure the serostatus risk table (q_all_regions / risk_components_all /
# risk_join_all) covers every coverage level present in postsim_vc_ixchiq_model.
# If not, rebuild them from scratch using the same logic as Section 8B.
# -----------------------------------------------------------------------------
ensure_risk_tables_all_coverages <- function() {
  stopifnot(
    exists("postsim_vc_ixchiq_model", envir = globalenv()),
    exists("calc_q_seromix_for_scenarios_agecat", envir = globalenv()),
    exists("age_map", envir = globalenv()),
    exists("risk_draw_df", envir = globalenv()),
    exists("rr_vals", envir = globalenv())
  )

  postsim <- get("postsim_vc_ixchiq_model", envir = globalenv())
  cov_levels <- unique(unlist(lapply(postsim, function(region_list) {
    unlist(lapply(region_list, names))
  })))
  cov_levels <- sort(unique(cov_levels))

  current_cov <- if (exists("risk_join_all", envir = globalenv())) {
    sort(unique(as.character(get("risk_join_all", envir = globalenv())$Coverage)))
  } else {
    character(0)
  }

  if (length(setdiff(cov_levels, current_cov)) == 0 && length(current_cov) > 0) {
    message("[ensure_risk_tables] risk_join_all already covers: ",
            paste(current_cov, collapse = ", "))
    return(invisible(NULL))
  }

  message("[ensure_risk_tables] Rebuilding q_all_regions, risk_components_all, ",
          "risk_join_all to include all coverages: ",
          paste(cov_levels, collapse = ", "),
          " (current: ",
          if (length(current_cov)) paste(current_cov, collapse = ", ") else "none",
          ")")

  q_all_regions_new <- purrr::imap_dfr(postsim, function(region_list, region_name) {
    purrr::imap_dfr(region_list, function(ve_list, ve_name) {
      purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
        q_tmp <- calc_q_seromix_for_scenarios_agecat(
          sim_region_ve_cov = cov_list$scenario_result,
          age_map = age_map
        )
        if (is.null(q_tmp) || nrow(q_tmp) == 0) {
          warning("calc_q_seromix returned empty rows for ",
                  region_name, " / ", ve_name, " / ", cov_name)
          return(NULL)
        }
        q_tmp %>%
          dplyr::mutate(Region = region_name, VE = ve_name, Coverage = cov_name) %>%
          dplyr::select(
            Region, VE, Coverage, Scenario, draw_id, AgeCat,
            total_vacc_age, seroneg_vacc_age, seropos_vacc_age,
            q_seroneg_vacc, q_seropos_vacc
          )
      })
    })
  })

  risk_components_all_new <- q_all_regions_new %>%
    tidyr::crossing(RR_seropos = rr_vals) %>%
    dplyr::left_join(risk_draw_df, by = c("draw_id", "AgeCat")) %>%
    dplyr::mutate(
      p_sae_vacc_seroneg   = p_sae_vacc_base,
      p_death_vacc_seroneg = p_death_vacc_base,
      p_sae_vacc_seropos   = p_sae_vacc_base   * RR_seropos,
      p_death_vacc_seropos = p_death_vacc_base * RR_seropos,
      p_sae_vacc_seroneg_contrib   = dplyr::if_else(total_vacc_age > 0, q_seroneg_vacc * p_sae_vacc_seroneg,   0),
      p_death_vacc_seroneg_contrib = dplyr::if_else(total_vacc_age > 0, q_seroneg_vacc * p_death_vacc_seroneg, 0),
      p_sae_vacc_seropos_contrib   = dplyr::if_else(total_vacc_age > 0, q_seropos_vacc * p_sae_vacc_seropos,   0),
      p_death_vacc_seropos_contrib = dplyr::if_else(total_vacc_age > 0, q_seropos_vacc * p_death_vacc_seropos, 0),
      p_sae_vacc_adj   = p_sae_vacc_seroneg_contrib   + p_sae_vacc_seropos_contrib,
      p_death_vacc_adj = p_death_vacc_seroneg_contrib + p_death_vacc_seropos_contrib,
      sae_10k_base   = 1e4 * p_sae_vacc_base   + 1e4 * p_death_vacc_base,
      death_10k_base = 1e4 * p_death_vacc_base,
      sae_10k_seroneg   = 1e4 * p_sae_vacc_seroneg_contrib   + 1e4 * p_death_vacc_seroneg_contrib,
      death_10k_seroneg = 1e4 * p_death_vacc_seroneg_contrib,
      sae_10k_seropos   = 1e4 * p_sae_vacc_seropos_contrib   + 1e4 * p_death_vacc_seropos_contrib,
      death_10k_seropos = 1e4 * p_death_vacc_seropos_contrib,
      sae_10k_adj   = sae_10k_seroneg   + sae_10k_seropos,
      death_10k_adj = death_10k_seroneg + death_10k_seropos
    )

  risk_join_all_new <- risk_components_all_new %>%
    dplyr::select(
      Region, VE, Coverage, Scenario, draw_id, AgeCat, RR_seropos,
      total_vacc_age, q_seroneg_vacc, q_seropos_vacc,
      p_sae_vacc_base, p_death_vacc_base,
      p_sae_vacc_seroneg, p_death_vacc_seroneg,
      p_sae_vacc_seropos, p_death_vacc_seropos,
      p_sae_vacc_seroneg_contrib, p_death_vacc_seroneg_contrib,
      p_sae_vacc_seropos_contrib, p_death_vacc_seropos_contrib,
      p_sae_vacc_adj, p_death_vacc_adj,
      sae_10k_base, death_10k_base,
      sae_10k_seroneg, death_10k_seroneg,
      sae_10k_seropos, death_10k_seropos,
      sae_10k_adj, death_10k_adj
    )

  assign("q_all_regions",       q_all_regions_new,       envir = globalenv())
  assign("risk_components_all", risk_components_all_new, envir = globalenv())
  assign("risk_join_all",       risk_join_all_new,       envir = globalenv())

  message("[ensure_risk_tables] Rebuilt. risk_join_all Coverage now: ",
          paste(sort(unique(risk_join_all_new$Coverage)), collapse = ", "))
  invisible(NULL)
}

ensure_risk_tables_all_coverages()

# -----------------------------------------------------------------------------
# Pipeline wrapper: runs Section 8B-8D for a single coverage level
# -----------------------------------------------------------------------------
run_brr_pipeline_for_coverage <- function(coverage_code,
                                          save_docx   = TRUE,
                                          results_dir = "06_Results") {
  stopifnot(coverage_code %in% c("cov10", "cov50", "cov90"))

  # ---- sanity check: required objects in global env ------------------------
  need <- c(
    "combined_nnv_df_region_coverage_model",
    "all_draws_daly_true", "all_draws_sae_true", "all_draws_fatal_true",
    "lhs_sample", "daly_pars_true", "risk_join_all",
    "setting_key", "rr_vals", "map_scenario_agecat_int",
    "make_averted_draws_true", "compute_daly_one",
    "fmt_ci", "make_pr_gt1_wide",
    "make_brr_ceac_outbreak", "plot_brr_ceac_outbreak_ve"
  )
  missing_obj <- need[!vapply(need, exists, logical(1), envir = globalenv())]
  if (length(missing_obj) > 0) {
    stop(
      "Missing required objects. Source brazil_all_draws_ori_v2.R first.\n",
      "  Missing: ", paste(missing_obj, collapse = ", ")
    )
  }

  # ---- diagnostic: does this coverage exist in the upstream objects? -------
  vc_vals <- unique(as.character(combined_nnv_df_region_coverage_model$VC))
  cov_vals_daly  <- unique(as.character(all_draws_daly_true$Coverage))
  cov_vals_sae   <- unique(as.character(all_draws_sae_true$Coverage))
  cov_vals_fatal <- unique(as.character(all_draws_fatal_true$Coverage))
  cov_vals_risk  <- unique(as.character(risk_join_all$Coverage))

  message(
    "[run_brr_pipeline_for_coverage] coverage_code = ", coverage_code, "\n",
    "  combined_nnv_df_region_coverage_model VC values: ",
      paste(vc_vals, collapse = ", "), "\n",
    "  all_draws_daly_true  Coverage values: ",
      paste(cov_vals_daly,  collapse = ", "), "\n",
    "  all_draws_sae_true   Coverage values: ",
      paste(cov_vals_sae,   collapse = ", "), "\n",
    "  all_draws_fatal_true Coverage values: ",
      paste(cov_vals_fatal, collapse = ", "), "\n",
    "  risk_join_all        Coverage values: ",
      paste(cov_vals_risk,  collapse = ", ")
  )

  if (!(coverage_code %in% vc_vals)) {
    stop("combined_nnv_df_region_coverage_model has no rows with VC == '",
         coverage_code, "'. Available VCs: ", paste(vc_vals, collapse = ", "))
  }
  if (!(coverage_code %in% cov_vals_daly)) {
    stop("all_draws_daly_true has no rows with Coverage == '",
         coverage_code, "'. Available Coverage: ",
         paste(cov_vals_daly, collapse = ", "))
  }
  if (!(coverage_code %in% cov_vals_risk)) {
    stop("risk_join_all has no rows with Coverage == '",
         coverage_code, "'. Available Coverage: ",
         paste(cov_vals_risk, collapse = ", "))
  }

  # -------------------------------------------------------------------------
  # Step 8.2) tot_vacc_map for this coverage
  # -------------------------------------------------------------------------
  tot_vacc_map_cov <- combined_nnv_df_region_coverage_model %>%
    dplyr::filter(VC == coverage_code) %>%
    dplyr::mutate(
      AgeCat = dplyr::case_when(
        AgeGroup %in% 2:4   ~ "1-11",
        AgeGroup == 5       ~ "12-17",
        AgeGroup %in% 6:15  ~ "18-64",
        AgeGroup %in% 16:20 ~ "65+",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(AgeCat)) %>%
    dplyr::group_by(scenario, region, AgeCat, VE) %>%
    dplyr::summarise(tot_vacc_grp = sum(tot_vacc, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(
      target = dplyr::case_when(
        scenario == "Scenario_1" & AgeCat == "1-11"  ~ 1L,
        scenario == "Scenario_2" & AgeCat == "12-17" ~ 1L,
        scenario == "Scenario_3" & AgeCat == "18-64" ~ 1L,
        scenario == "Scenario_4" & AgeCat == "65+"   ~ 1L,
        TRUE ~ 0L
      )
    ) %>%
    dplyr::filter(target == 1L) %>%
    dplyr::transmute(
      scenario = as.integer(gsub("Scenario_", "", scenario)),
      AgeCat, VE, tot_vacc_grp, region
    ) %>%
    dplyr::rename(Region = region, Scenario = scenario)

  # -------------------------------------------------------------------------
  # Step 8.3) Benefit (averted burden) draws at this coverage
  # -------------------------------------------------------------------------
  benefit_draws_cov <- dplyr::bind_rows(
    make_averted_draws_true(all_draws_daly_true,  "DALY",  coverage_keep = coverage_code),
    make_averted_draws_true(all_draws_sae_true,   "SAE",   coverage_keep = coverage_code),
    make_averted_draws_true(all_draws_fatal_true, "Death", coverage_keep = coverage_code)
  )

  if (nrow(benefit_draws_cov) == 0) {
    stop("benefit_draws is empty for coverage '", coverage_code,
         "'. Check that all_draws_*_true contain rows with this Coverage.")
  }
  if (nrow(tot_vacc_map_cov) == 0) {
    stop("tot_vacc_map is empty for coverage '", coverage_code,
         "'. Check combined_nnv_df_region_coverage_model$VC values.")
  }

  # -------------------------------------------------------------------------
  # Step 8.6) Benefit base (attach DALY pars and age category)
  # -------------------------------------------------------------------------
  benefit_base_cov <- benefit_draws_cov %>%
    dplyr::mutate(
      Scenario  = as.integer(Scenario),
      AgeCat    = map_scenario_agecat_int(Scenario),
      risk_band = dplyr::if_else(AgeCat == "65+", "65+", "u65")
    ) %>%
    dplyr::left_join(daly_pars_true, by = "draw_id")

  # -------------------------------------------------------------------------
  # Step 8.10.6) Draw-level burden table + join with serostatus risk
  # -------------------------------------------------------------------------
  benefit_draw_df_cov <- benefit_base_cov %>%
    dplyr::left_join(
      tot_vacc_map_cov,
      by = c("Region", "Scenario", "VE", "AgeCat")
    ) %>%
    dplyr::mutate(
      Scenario      = as.integer(Scenario),
      baseline_10k  = (baseline / tot_vacc_grp) * 1e4,
      post_10k      = (post     / tot_vacc_grp) * 1e4,
      averted_10k   = (averted  / tot_vacc_grp) * 1e4,
      pct_reduction = dplyr::if_else(
        baseline > 0,
        100 * averted / baseline,
        NA_real_
      )
    ) %>%
    tidyr::crossing(RR_seropos = rr_vals)

  draw_level_xy_serostatus_cov <- benefit_draw_df_cov %>%
    dplyr::left_join(
      risk_join_all,
      by = c("Region", "VE", "Coverage", "Scenario", "draw_id",
             "AgeCat", "RR_seropos")
    )

  # -------------------------------------------------------------------------
  # Step 8.10.7) Compute vaccine-caused DALY components
  # -------------------------------------------------------------------------
  draw_level_xy_serostatus_cov <- draw_level_xy_serostatus_cov %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      daly_10k_base = if (outcome == "DALY") {
        compute_daly_one_age_specific(
          age_group      = AgeCat,
          sae_10k        = sae_10k_base,
          deaths_sae_10k = death_10k_base,
          draw_id        = draw_id,
          draw_pars = list(
            le_lost_1_11  = le_lost_1_11,
            le_lost_12_17 = le_lost_12_17,
            le_lost_18_64 = le_lost_18_64,
            le_lost_65    = le_lost_65,
            dw_hosp       = dw_hosp,
            dw_nonhosp    = dw_nonhosp,
            dw_subac      = dw_subac,
            dw_chronic    = dw_chronic,
            dur_acute     = dur_acute,
            dur_nonhosp   = dur_nonhosp,
            dur_subac     = dur_subac,
            dur_6m        = dur_6m,
            dur_12m       = dur_12m,
            dur_30m       = dur_30m,
            acute         = acute,
            subac         = subac,
            chr6m         = chr6m,
            chr12m        = chr12m,
            chr30m        = chr30m
          )
        )$daly_sae
      } else { NA_real_ },

      daly_10k_seroneg = if (outcome == "DALY") {
        compute_daly_one_age_specific(
          age_group      = AgeCat,
          sae_10k        = sae_10k_seroneg,
          deaths_sae_10k = death_10k_seroneg,
          draw_id        = draw_id,
          draw_pars = list(
            le_lost_1_11  = le_lost_1_11,
            le_lost_12_17 = le_lost_12_17,
            le_lost_18_64 = le_lost_18_64,
            le_lost_65    = le_lost_65,
            dw_hosp       = dw_hosp,
            dw_nonhosp    = dw_nonhosp,
            dw_subac      = dw_subac,
            dw_chronic    = dw_chronic,
            dur_acute     = dur_acute,
            dur_nonhosp   = dur_nonhosp,
            dur_subac     = dur_subac,
            dur_6m        = dur_6m,
            dur_12m       = dur_12m,
            dur_30m       = dur_30m,
            acute         = acute,
            subac         = subac,
            chr6m         = chr6m,
            chr12m        = chr12m,
            chr30m        = chr30m
          )
        )$daly_sae
      } else { NA_real_ },

      daly_10k_seropos = if (outcome == "DALY") {
        compute_daly_one_age_specific(
          age_group      = AgeCat,
          sae_10k        = sae_10k_seropos,
          deaths_sae_10k = death_10k_seropos,
          draw_id        = draw_id,
          draw_pars = list(
            le_lost_1_11  = le_lost_1_11,
            le_lost_12_17 = le_lost_12_17,
            le_lost_18_64 = le_lost_18_64,
            le_lost_65    = le_lost_65,
            dw_hosp       = dw_hosp,
            dw_nonhosp    = dw_nonhosp,
            dw_subac      = dw_subac,
            dw_chronic    = dw_chronic,
            dur_acute     = dur_acute,
            dur_nonhosp   = dur_nonhosp,
            dur_subac     = dur_subac,
            dur_6m        = dur_6m,
            dur_12m       = dur_12m,
            dur_30m       = dur_30m,
            acute         = acute,
            subac         = subac,
            chr6m         = chr6m,
            chr12m        = chr12m,
            chr30m        = chr30m
          )
        )$daly_sae
      } else { NA_real_ },

      daly_10k_adj = if (outcome == "DALY") {
        daly_10k_seroneg + daly_10k_seropos
      } else { NA_real_ }
    ) %>%
    dplyr::ungroup()

  # -------------------------------------------------------------------------
  # Step 8.10.8) Reporting labels + BRR columns
  # -------------------------------------------------------------------------
  draw_level_xy_serostatus_cov <- draw_level_xy_serostatus_cov %>%
    dplyr::mutate(
      VE_label = factor(
        VE,
        levels = c("VE0", "VE98.9"),
        labels = c("Disease blocking only", "Disease and infection blocking")
      ),
      setting = unname(setting_key[Region])
    ) %>%
    dplyr::mutate(
      x_10k_base = dplyr::case_when(
        outcome == "SAE"   ~ sae_10k_base,
        outcome == "Death" ~ death_10k_base,
        outcome == "DALY"  ~ daly_10k_base,
        TRUE ~ NA_real_
      ),
      x_10k_adj = dplyr::case_when(
        outcome == "SAE"   ~ sae_10k_adj,
        outcome == "Death" ~ death_10k_adj,
        outcome == "DALY"  ~ daly_10k_adj,
        TRUE ~ NA_real_
      ),
      brr_base = dplyr::if_else(
        is.na(x_10k_base) | x_10k_base == 0,
        NA_real_,
        averted_10k / x_10k_base
      ),
      brr_adj = dplyr::if_else(
        is.na(x_10k_adj) | x_10k_adj == 0,
        NA_real_,
        averted_10k / x_10k_adj
      )
    )

  # -------------------------------------------------------------------------
  # Step 8.10.9) Summary across all age groups
  # -------------------------------------------------------------------------
  brr_draw_summary_cov <- draw_level_xy_serostatus_cov %>%
    dplyr::mutate(
      brr_base = ifelse(is.infinite(brr_base), NA, brr_base),
      brr_adj  = ifelse(is.infinite(brr_adj),  NA, brr_adj),
      setting  = factor(setting, levels = c("Low", "Moderate", "High"))
    ) %>%
    dplyr::group_by(outcome, Scenario, AgeCat, VE_label, RR_seropos, setting) %>%
    dplyr::summarise(
      brr_base_med = quantile(brr_base, 0.50,  na.rm = TRUE),
      brr_base_lo  = quantile(brr_base, 0.025, na.rm = TRUE),
      brr_base_hi  = quantile(brr_base, 0.975, na.rm = TRUE),

      brr_adj_med  = quantile(brr_adj,  0.50,  na.rm = TRUE),
      brr_adj_lo   = quantile(brr_adj,  0.025, na.rm = TRUE),
      brr_adj_hi   = quantile(brr_adj,  0.975, na.rm = TRUE),

      av_med = quantile(averted_10k, 0.50,  na.rm = TRUE),
      av_lo  = quantile(averted_10k, 0.025, na.rm = TRUE),
      av_hi  = quantile(averted_10k, 0.975, na.rm = TRUE),

      ca_sae_base_med   = quantile(sae_10k_base,   0.50,  na.rm = TRUE),
      ca_sae_base_lo    = quantile(sae_10k_base,   0.025, na.rm = TRUE),
      ca_sae_base_hi    = quantile(sae_10k_base,   0.975, na.rm = TRUE),

      ca_death_base_med = quantile(death_10k_base, 0.50,  na.rm = TRUE),
      ca_death_base_lo  = quantile(death_10k_base, 0.025, na.rm = TRUE),
      ca_death_base_hi  = quantile(death_10k_base, 0.975, na.rm = TRUE),

      ca_daly_base_med  = quantile(daly_10k_base,  0.50,  na.rm = TRUE),
      ca_daly_base_lo   = quantile(daly_10k_base,  0.025, na.rm = TRUE),
      ca_daly_base_hi   = quantile(daly_10k_base,  0.975, na.rm = TRUE),

      ca_sae_adj_med    = quantile(sae_10k_adj,    0.50,  na.rm = TRUE),
      ca_sae_adj_lo     = quantile(sae_10k_adj,    0.025, na.rm = TRUE),
      ca_sae_adj_hi     = quantile(sae_10k_adj,    0.975, na.rm = TRUE),

      ca_death_adj_med  = quantile(death_10k_adj,  0.50,  na.rm = TRUE),
      ca_death_adj_lo   = quantile(death_10k_adj,  0.025, na.rm = TRUE),
      ca_death_adj_hi   = quantile(death_10k_adj,  0.975, na.rm = TRUE),

      ca_daly_adj_med   = quantile(daly_10k_adj,   0.50,  na.rm = TRUE),
      ca_daly_adj_lo    = quantile(daly_10k_adj,   0.025, na.rm = TRUE),
      ca_daly_adj_hi    = quantile(daly_10k_adj,   0.975, na.rm = TRUE),

      .groups = "drop"
    ) %>%
    dplyr::mutate(
      scenario  = Scenario,
      age_group = AgeCat
    ) %>%
    dplyr::arrange(outcome, Scenario, setting, AgeCat, VE_label, RR_seropos)

  # Step 8.10.10) Filter to >=18 and RR_seropos = 0
  brr_draw_summary_filtered_cov <- brr_draw_summary_cov %>%
    dplyr::filter(
      AgeCat %in% c("18-64", "65+"),
      RR_seropos == 0.0
    ) %>%
    dplyr::select(-scenario, -age_group)

  # -------------------------------------------------------------------------
  # Step 8.10.11-12) Wide-format tables
  # -------------------------------------------------------------------------
  ca_summary_cov <- brr_draw_summary_filtered_cov %>%
    dplyr::filter(outcome == "DALY") %>%
    dplyr::select(Scenario, AgeCat, VE_label, RR_seropos, setting,
                  dplyr::starts_with("ca_"))

  brr_table_long_serostatus_cov <- brr_draw_summary_filtered_cov %>%
    dplyr::select(outcome, Scenario, AgeCat, VE_label, RR_seropos, setting,
                  brr_base_med, brr_base_lo, brr_base_hi,
                  brr_adj_med,  brr_adj_lo,  brr_adj_hi,
                  av_med, av_lo, av_hi) %>%
    dplyr::left_join(
      ca_summary_cov,
      by = c("Scenario", "AgeCat", "VE_label", "RR_seropos", "setting")
    ) %>%
    dplyr::mutate(
      brr_base_formatted      = sprintf("%.2f (%.2f\u2013%.2f)", brr_base_med, brr_base_lo, brr_base_hi),
      brr_adj_formatted       = sprintf("%.2f (%.2f\u2013%.2f)", brr_adj_med,  brr_adj_lo,  brr_adj_hi),
      av_formatted            = sprintf("%.2f (%.2f\u2013%.2f)", av_med,       av_lo,       av_hi),
      ca_sae_base_formatted   = sprintf("%.2f (%.2f\u2013%.2f)", ca_sae_base_med,   ca_sae_base_lo,   ca_sae_base_hi),
      ca_death_base_formatted = sprintf("%.2f (%.2f\u2013%.2f)", ca_death_base_med, ca_death_base_lo, ca_death_base_hi),
      ca_daly_base_formatted  = sprintf("%.2f (%.2f\u2013%.2f)", ca_daly_base_med,  ca_daly_base_lo,  ca_daly_base_hi),
      ca_sae_adj_formatted    = sprintf("%.2f (%.2f\u2013%.2f)", ca_sae_adj_med,    ca_sae_adj_lo,    ca_sae_adj_hi),
      ca_death_adj_formatted  = sprintf("%.2f (%.2f\u2013%.2f)", ca_death_adj_med,  ca_death_adj_lo,  ca_death_adj_hi),
      ca_daly_adj_formatted   = sprintf("%.2f (%.2f\u2013%.2f)", ca_daly_adj_med,   ca_daly_adj_lo,   ca_daly_adj_hi)
    )

  # BRR long-form for CEAC
  brr_long_serostatus_cov <- draw_level_xy_serostatus_cov %>%
    dplyr::mutate(
      outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death", "daly" = "DALY"),
      setting = factor(setting, levels = c("Low", "Moderate", "High")),
      AgeCat  = factor(AgeCat,  levels = c("18-64", "65+"))
    ) %>%
    dplyr::filter(
      is.finite(brr_base), brr_base > 0,
      is.finite(brr_adj),  brr_adj  > 0
    ) %>%
    dplyr::transmute(
      Region, setting, VE_label, AgeCat, outcome, RR_seropos,
      brr_base, brr_adj
    )

  brr_long_serostatus_long_cov <- brr_long_serostatus_cov %>%
    dplyr::filter(!is.na(AgeCat)) %>%
    tidyr::pivot_longer(
      cols      = c(brr_base, brr_adj),
      names_to  = "brr_type",
      values_to = "brr"
    ) %>%
    dplyr::filter(is.finite(brr), brr > 0)

  brr_rr0 <- brr_long_serostatus_long_cov %>%
    dplyr::filter(RR_seropos == 0)

  if (nrow(brr_rr0) == 0) {
    warning(
      "[coverage = ", coverage_code,
      "] brr_long_serostatus_long has no finite positive BRR rows with ",
      "RR_seropos == 0. Falling back to default threshold range ",
      "10^-1 .. 10^1 for CEAC. Check that tot_vacc_map joined correctly ",
      "and that averted/risk draws are non-zero."
    )
    brr_range_cov <- data.frame(min_brr = 0.1, max_brr = 10)
    thresholds_auto_cov <- 10^seq(-1, 1, by = 0.02)
  } else {
    brr_range_cov <- brr_rr0 %>%
      dplyr::summarise(
        min_brr = min(brr, na.rm = TRUE),
        max_brr = max(brr, na.rm = TRUE)
      )
    if (!is.finite(brr_range_cov$min_brr) ||
        !is.finite(brr_range_cov$max_brr)) {
      warning("[coverage = ", coverage_code,
              "] non-finite BRR range; using default 10^-1..10^1.")
      brr_range_cov <- data.frame(min_brr = 0.1, max_brr = 10)
      thresholds_auto_cov <- 10^seq(-1, 1, by = 0.02)
    } else {
      thresholds_auto_cov <- 10^seq(
        floor(log10(brr_range_cov$min_brr)),
        ceiling(log10(brr_range_cov$max_brr)),
        by = 0.02
      )
    }
  }

  ceac_ob_cov <- make_brr_ceac_outbreak(
    brr_long_serostatus_long_cov,
    group_vars = c("setting", "VE_label", "AgeCat", "outcome",
                   "RR_seropos", "brr_type")
  )

  pr_gt1_wide_setting_base_cov <- make_pr_gt1_wide(ceac_ob_cov, "brr_base")
  pr_gt1_wide_setting_adj_cov  <- make_pr_gt1_wide(ceac_ob_cov, "brr_adj")

  # Wide-format reporting table (RR_seropos == 0)
  ca_summary_wide <- brr_draw_summary_filtered_cov %>%
    dplyr::filter(RR_seropos == 0, outcome == "DALY") %>%
    dplyr::select(Scenario, AgeCat, VE_label, setting,
                  ca_sae_base_med,   ca_sae_base_lo,   ca_sae_base_hi,
                  ca_death_base_med, ca_death_base_lo, ca_death_base_hi,
                  ca_daly_base_med,  ca_daly_base_lo,  ca_daly_base_hi,
                  ca_sae_adj_med,    ca_sae_adj_lo,    ca_sae_adj_hi,
                  ca_death_adj_med,  ca_death_adj_lo,  ca_death_adj_hi,
                  ca_daly_adj_med,   ca_daly_adj_lo,   ca_daly_adj_hi)

  brr_draw_for_wide_cov <- brr_draw_summary_filtered_cov %>%
    dplyr::filter(RR_seropos == 0) %>%
    dplyr::select(outcome, Scenario, AgeCat, VE_label, setting,
                  av_med,       av_lo,       av_hi,
                  brr_base_med, brr_base_lo, brr_base_hi,
                  brr_adj_med,  brr_adj_lo,  brr_adj_hi) %>%
    dplyr::left_join(ca_summary_wide,
                     by = c("Scenario", "AgeCat", "VE_label", "setting"))

  brr_table_wide_serostatus_cov <- brr_draw_for_wide_cov %>%
    dplyr::mutate(
      setting = factor(setting, levels = c("Low", "Moderate", "High")),

      av_formatted            = fmt_ci(av_med,       av_lo,       av_hi),
      brr_base_formatted      = fmt_ci(brr_base_med, brr_base_lo, brr_base_hi),
      brr_adj_formatted       = fmt_ci(brr_adj_med,  brr_adj_lo,  brr_adj_hi),
      ca_sae_base_formatted   = fmt_ci(ca_sae_base_med,   ca_sae_base_lo,   ca_sae_base_hi),
      ca_death_base_formatted = fmt_ci(ca_death_base_med, ca_death_base_lo, ca_death_base_hi),
      ca_daly_base_formatted  = fmt_ci(ca_daly_base_med,  ca_daly_base_lo,  ca_daly_base_hi),
      ca_sae_adj_formatted    = fmt_ci(ca_sae_adj_med,    ca_sae_adj_lo,    ca_sae_adj_hi),
      ca_death_adj_formatted  = fmt_ci(ca_death_adj_med,  ca_death_adj_lo,  ca_death_adj_hi),
      ca_daly_adj_formatted   = fmt_ci(ca_daly_adj_med,   ca_daly_adj_lo,   ca_daly_adj_hi)
    ) %>%
    dplyr::select(outcome, Scenario, AgeCat, setting, VE_label,
                  av_formatted, brr_base_formatted, brr_adj_formatted,
                  ca_sae_base_formatted, ca_death_base_formatted, ca_daly_base_formatted,
                  ca_sae_adj_formatted,  ca_death_adj_formatted,  ca_daly_adj_formatted) %>%
    tidyr::pivot_wider(
      names_from  = VE_label,
      values_from = c(av_formatted, brr_base_formatted, brr_adj_formatted,
                      ca_sae_base_formatted, ca_death_base_formatted, ca_daly_base_formatted,
                      ca_sae_adj_formatted,  ca_death_adj_formatted,  ca_daly_adj_formatted),
      names_glue  = "{VE_label}_{.value}"
    ) %>%
    dplyr::arrange(outcome, Scenario, setting, AgeCat) %>%
    dplyr::rename(
      Outcome     = outcome,
      `Age group` = AgeCat,
      Setting     = setting
    )

  brr_table_wide_serostatus_cov <- brr_table_wide_serostatus_cov %>%
    dplyr::rename(
      `DB (Averted)`            = `Disease blocking only_av_formatted`,
      `DB (BRR base)`           = `Disease blocking only_brr_base_formatted`,
      `DB (BRR adj)`            = `Disease blocking only_brr_adj_formatted`,
      `DB (SAE caused base)`    = `Disease blocking only_ca_sae_base_formatted`,
      `DB (SAE caused adj)`     = `Disease blocking only_ca_sae_adj_formatted`,
      `DB (Death caused base)`  = `Disease blocking only_ca_death_base_formatted`,
      `DB (Death caused adj)`   = `Disease blocking only_ca_death_adj_formatted`,
      `DB (DALY caused base)`   = `Disease blocking only_ca_daly_base_formatted`,
      `DB (DALY caused adj)`    = `Disease blocking only_ca_daly_adj_formatted`,
      `DIB (Averted)`           = `Disease and infection blocking_av_formatted`,
      `DIB (BRR base)`          = `Disease and infection blocking_brr_base_formatted`,
      `DIB (BRR adj)`           = `Disease and infection blocking_brr_adj_formatted`,
      `DIB (SAE caused base)`   = `Disease and infection blocking_ca_sae_base_formatted`,
      `DIB (SAE caused adj)`    = `Disease and infection blocking_ca_sae_adj_formatted`,
      `DIB (Death caused base)` = `Disease and infection blocking_ca_death_base_formatted`,
      `DIB (Death caused adj)`  = `Disease and infection blocking_ca_death_adj_formatted`,
      `DIB (DALY caused base)`  = `Disease and infection blocking_ca_daly_base_formatted`,
      `DIB (DALY caused adj)`   = `Disease and infection blocking_ca_daly_adj_formatted`
    )

  pr_gt1_wide_serostatus_base_cov <- pr_gt1_wide_setting_base_cov %>%
    dplyr::rename_with(
      ~ gsub("Pr(BRR>1)", "Pr(BRR_base>1)", .x, fixed = TRUE),
      .cols = -c(Outcome, Setting, `Age group`)
    )

  pr_gt1_wide_serostatus_adj_cov <- pr_gt1_wide_setting_adj_cov %>%
    dplyr::rename_with(
      ~ gsub("Pr(BRR>1)", "Pr(BRR_adj>1)", .x, fixed = TRUE),
      .cols = -c(Outcome, Setting, `Age group`)
    )

  brr_table_wide_serostatus2_cov <- brr_table_wide_serostatus_cov %>%
    dplyr::left_join(pr_gt1_wide_serostatus_base_cov,
                     by = c("Outcome", "Setting", "Age group")) %>%
    dplyr::left_join(pr_gt1_wide_serostatus_adj_cov,
                     by = c("Outcome", "Setting", "Age group")) %>%
    dplyr::relocate(`DB (Averted)`,                                  .after = `Age group`) %>%
    dplyr::relocate(`DB (BRR base)`,                                 .after = `DB (Averted)`) %>%
    dplyr::relocate(`DB (BRR adj)`,                                  .after = `DB (BRR base)`) %>%
    dplyr::relocate(`Disease blocking only Pr(BRR_base>1)`,          .after = `DB (BRR adj)`) %>%
    dplyr::relocate(`Disease blocking only Pr(BRR_adj>1)`,           .after = `Disease blocking only Pr(BRR_base>1)`) %>%
    dplyr::relocate(`DB (SAE caused base)`,                          .after = `Disease blocking only Pr(BRR_adj>1)`) %>%
    dplyr::relocate(`DB (SAE caused adj)`,                           .after = `DB (SAE caused base)`) %>%
    dplyr::relocate(`DB (Death caused base)`,                        .after = `DB (SAE caused adj)`) %>%
    dplyr::relocate(`DB (Death caused adj)`,                         .after = `DB (Death caused base)`) %>%
    dplyr::relocate(`DB (DALY caused base)`,                         .after = `DB (Death caused adj)`) %>%
    dplyr::relocate(`DB (DALY caused adj)`,                          .after = `DB (DALY caused base)`) %>%
    dplyr::relocate(`DIB (Averted)`,                                 .after = `DB (DALY caused adj)`) %>%
    dplyr::relocate(`DIB (BRR base)`,                                .after = `DIB (Averted)`) %>%
    dplyr::relocate(`DIB (BRR adj)`,                                 .after = `DIB (BRR base)`) %>%
    dplyr::relocate(`Disease and infection blocking Pr(BRR_base>1)`, .after = `DIB (BRR adj)`) %>%
    dplyr::relocate(`Disease and infection blocking Pr(BRR_adj>1)`,  .after = `Disease and infection blocking Pr(BRR_base>1)`) %>%
    dplyr::relocate(`DIB (SAE caused base)`,                         .after = `Disease and infection blocking Pr(BRR_adj>1)`) %>%
    dplyr::relocate(`DIB (SAE caused adj)`,                          .after = `DIB (SAE caused base)`) %>%
    dplyr::relocate(`DIB (Death caused base)`,                       .after = `DIB (SAE caused adj)`) %>%
    dplyr::relocate(`DIB (Death caused adj)`,                        .after = `DIB (Death caused base)`) %>%
    dplyr::relocate(`DIB (DALY caused base)`,                        .after = `DIB (Death caused adj)`) %>%
    dplyr::relocate(`DIB (DALY caused adj)`,                         .after = `DIB (DALY caused base)`) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        ~ tidyr::replace_na(as.character(.x), "beneficial")
      )
    )

  # -------------------------------------------------------------------------
  # SECTION 08C. Reporting table assembly
  # -------------------------------------------------------------------------
  part1 <- brr_table_wide_serostatus2_cov %>%
    dplyr::mutate(
      Risk_base = dplyr::case_when(
        Outcome == "SAE"   ~ `DB (SAE caused base)`,
        Outcome == "Death" ~ `DB (Death caused base)`,
        Outcome == "DALY"  ~ `DB (DALY caused base)`,
        TRUE ~ NA_character_
      ),
      Risk_adj = dplyr::case_when(
        Outcome == "SAE"   ~ `DB (SAE caused adj)`,
        Outcome == "Death" ~ `DB (Death caused adj)`,
        Outcome == "DALY"  ~ `DB (DALY caused adj)`,
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(
      Outcome, Setting, `Age group`,
      Benefit = `DB (Averted)`,
      Risk_base, Risk_adj,
      BRR_base  = `DB (BRR base)`,
      BRR_adj   = `DB (BRR adj)`,
      prob_base = `Disease blocking only Pr(BRR_base>1)`,
      prob_adj  = `Disease blocking only Pr(BRR_adj>1)`
    ) %>%
    dplyr::mutate(mechanism = "Disease blocking only")

  part2 <- brr_table_wide_serostatus2_cov %>%
    dplyr::mutate(
      Risk_base = dplyr::case_when(
        Outcome == "SAE"   ~ `DIB (SAE caused base)`,
        Outcome == "Death" ~ `DIB (Death caused base)`,
        Outcome == "DALY"  ~ `DIB (DALY caused base)`,
        TRUE ~ NA_character_
      ),
      Risk_adj = dplyr::case_when(
        Outcome == "SAE"   ~ `DIB (SAE caused adj)`,
        Outcome == "Death" ~ `DIB (Death caused adj)`,
        Outcome == "DALY"  ~ `DIB (DALY caused adj)`,
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(
      Outcome, Setting, `Age group`,
      Benefit = `DIB (Averted)`,
      Risk_base, Risk_adj,
      BRR_base  = `DIB (BRR base)`,
      BRR_adj   = `DIB (BRR adj)`,
      prob_base = `Disease and infection blocking Pr(BRR_base>1)`,
      prob_adj  = `Disease and infection blocking Pr(BRR_adj>1)`
    ) %>%
    dplyr::mutate(mechanism = "Disease and infection blocking")

  brr_table_final_long_cov <- dplyr::bind_rows(part1, part2) %>%
    dplyr::mutate(
      Setting = factor(Setting, levels = c("High", "Moderate", "Low"))
    ) %>%
    dplyr::arrange(Outcome, Setting, `Age group`) %>%
    dplyr::select(
      Outcome, Setting, `Age group`, mechanism,
      Benefit, Risk_base, Risk_adj, BRR_base, BRR_adj, prob_base, prob_adj
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(Benefit, Risk_base, Risk_adj, BRR_base, BRR_adj),
        ~ {
          x <- gsub(" \\(", "\n(", .x, fixed = FALSE)
          dplyr::case_when(
            x %in% c(
              "NA (NA-NA)", "NA (NA\u2013NA)",
              "NA\n(NA-NA)", "NA\n(NA\u2013NA)",
              "NA (NA NA)", "NA\n(NA NA)"
            ) ~ "beneficial",
            TRUE ~ x
          )
        }
      )
    )

  ft_brr_cov <- flextable::flextable(brr_table_final_long_cov) %>%
    flextable::set_header_labels(
      Outcome     = "Outcome",
      Setting     = "Setting",
      `Age group` = "Age group",
      mechanism   = "Vaccine protection\nmechanism",
      Benefit     = "Benefit:\nOutcomes averted\n(per 10,000)",
      Risk_base   = "Risk (base):\nOutcomes attributable\n(per 10,000)",
      Risk_adj    = "Risk (adjusted):\nOutcomes attributable\n(per 10,000)",
      BRR_base    = "BRR (base):\n(Prevented per 1 caused)",
      BRR_adj     = "BRR (adjusted):\n(Prevented per 1 caused)",
      prob_base   = "Probability\n(BRR_base > 1)\n(%)",
      prob_adj    = "Probability\n(BRR_adj > 1)\n(%)"
    ) %>%
    flextable::theme_booktabs() %>%
    flextable::bold(part = "header") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::align(j = 1:4, align = "left", part = "all") %>%
    flextable::merge_v(j = c("Outcome", "Setting", "Age group")) %>%
    flextable::valign(j = c("Outcome", "Setting", "Age group"), valign = "top") %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::autofit()

  if (save_docx) {
    if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
    out_docx <- file.path(
      results_dir,
      paste0("BRR_table_ori_setting_", coverage_code, ".docx")
    )
    doc <- officer::read_docx() %>%
      officer::body_add_par(
        paste0("Benefit-Risk Ratio (BRR) by Outcome, Age group, and VE (",
               coverage_code, ")"),
        style = "heading 2"
      ) %>%
      flextable::body_add_flextable(ft_brr_cov)
    print(doc, target = out_docx)
    message("Saved: ", out_docx)
  }

  # -------------------------------------------------------------------------
  # SECTION 08D. BRR acceptability curves
  # -------------------------------------------------------------------------
  p_daly_ceac_outbreak <- plot_brr_ceac_outbreak_ve(ceac_ob_cov, "DALY") +
    ggplot2::theme(text = ggplot2::element_text(family = "Calibri")) +
    ggplot2::labs(tag = "E") +
    ggplot2::theme(
      plot.tag          = ggplot2::element_text(face = "bold", size = 16),
      plot.tag.position = c(0, 1)
    )

  plot_sae <- plot_brr_ceac_outbreak_ve(ceac_ob_cov, "SAE") +
    ggplot2::theme(text = ggplot2::element_text(family = "Calibri")) +
    ggplot2::labs(tag = "E") +
    ggplot2::theme(
      plot.tag          = ggplot2::element_text(face = "bold", size = 16),
      plot.tag.position = c(0, 1)
    )

  plot_death <- plot_brr_ceac_outbreak_ve(ceac_ob_cov, "Death") +
    ggplot2::theme(text = ggplot2::element_text(family = "Calibri")) +
    ggplot2::labs(tag = "E") +
    ggplot2::theme(
      plot.tag          = ggplot2::element_text(face = "bold", size = 16),
      plot.tag.position = c(0, 1)
    )

  # -------------------------------------------------------------------------
  invisible(list(
    coverage                    = coverage_code,
    tot_vacc_map                = tot_vacc_map_cov,
    benefit_draws               = benefit_draws_cov,
    benefit_base                = benefit_base_cov,
    benefit_draw_df             = benefit_draw_df_cov,
    draw_level_xy_serostatus    = draw_level_xy_serostatus_cov,
    brr_draw_summary            = brr_draw_summary_cov,
    brr_draw_summary_filtered   = brr_draw_summary_filtered_cov,
    brr_table_long_serostatus   = brr_table_long_serostatus_cov,
    brr_long_serostatus         = brr_long_serostatus_cov,
    brr_long_serostatus_long    = brr_long_serostatus_long_cov,
    ceac_ob                     = ceac_ob_cov,
    pr_gt1_wide_setting_base    = pr_gt1_wide_setting_base_cov,
    pr_gt1_wide_setting_adj     = pr_gt1_wide_setting_adj_cov,
    brr_table_wide_serostatus   = brr_table_wide_serostatus_cov,
    brr_table_wide_serostatus2  = brr_table_wide_serostatus2_cov,
    brr_table_final_long        = brr_table_final_long_cov,
    ft_brr                      = ft_brr_cov,
    p_daly_ceac_outbreak        = p_daly_ceac_outbreak,
    plot_sae                    = plot_sae,
    plot_death                  = plot_death
  ))
}

# =============================================================================
# Run the pipeline for cov10 and cov90
# =============================================================================
ori_cov10 <- run_brr_pipeline_for_coverage("cov10")
ori_cov90 <- run_brr_pipeline_for_coverage("cov90")

# Convenience aliases mirroring the cov50 object names, with a _cov10/_cov90 suffix
tot_vacc_map_true_cov10                <- ori_cov10$tot_vacc_map
benefit_draws_true_cov10               <- ori_cov10$benefit_draws
benefit_base_true_cov10                <- ori_cov10$benefit_base
benefit_draw_df_cov10                  <- ori_cov10$benefit_draw_df
draw_level_xy_serostatus_cov10         <- ori_cov10$draw_level_xy_serostatus
brr_draw_summary_true_cov10            <- ori_cov10$brr_draw_summary
brr_draw_summary_true_filtered_cov10   <- ori_cov10$brr_draw_summary_filtered
ceac_ob_cov10                          <- ori_cov10$ceac_ob
brr_table_wide_serostatus2_cov10       <- ori_cov10$brr_table_wide_serostatus2
brr_table_final_long_cov10             <- ori_cov10$brr_table_final_long
ft_brr_cov10                           <- ori_cov10$ft_brr

tot_vacc_map_true_cov90                <- ori_cov90$tot_vacc_map
benefit_draws_true_cov90               <- ori_cov90$benefit_draws
benefit_base_true_cov90                <- ori_cov90$benefit_base
benefit_draw_df_cov90                  <- ori_cov90$benefit_draw_df
draw_level_xy_serostatus_cov90         <- ori_cov90$draw_level_xy_serostatus
brr_draw_summary_true_cov90            <- ori_cov90$brr_draw_summary
brr_draw_summary_true_filtered_cov90   <- ori_cov90$brr_draw_summary_filtered
ceac_ob_cov90                          <- ori_cov90$ceac_ob
brr_table_wide_serostatus2_cov90       <- ori_cov90$brr_table_wide_serostatus2
brr_table_final_long_cov90             <- ori_cov90$brr_table_final_long
ft_brr_cov90                           <- ori_cov90$ft_brr

message(
  "\n[ORI cov_alt] Finished.\n",
  "  cov10 rows (draw-level): ", nrow(draw_level_xy_serostatus_cov10), "\n",
  "  cov90 rows (draw-level): ", nrow(draw_level_xy_serostatus_cov90), "\n",
  "  Tables saved to 06_Results/BRR_table_ori_setting_cov10.docx and ..._cov90.docx"
)
