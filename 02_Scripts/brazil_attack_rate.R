get_total_S0 <- function(sim_results, t0 = 1) {
  out <- data.frame()
  
  for (region in names(sim_results)) {
    for (VE in names(sim_results[[region]])) {
      for (VC in names(sim_results[[region]][[VE]])) {
        
        scenario_list <- sim_results[[region]][[VE]][[VC]]
        
        for (sc in seq_along(scenario_list)) {
          S <- scenario_list[[sc]]$sim_out$S   # matrix: age x time
          total_S0 <- sum(S[, t0], na.rm = TRUE)
          
          out <- rbind(out, data.frame(
            region = region,
            VE = VE,
            VC = VC,
            scenario = paste0("Scenario_", sc),
            t0 = t0,
            total_S0 = total_S0
          ))
        }
      }
    }
  }
  
  out
}

S0_df <- get_total_S0(sim_results_vc_ixchiq_model, t0 = 1)

S0_df <- S0_df %>% distinct(region, total_S0)

p_symp <- 0.5242478

true_inf_region <- combined_nnv_df_region_coverage_model %>%
  distinct(region, pre_vacc, pre_vacc_lo, pre_vacc_hi) %>%
  group_by(region) %>%
  summarise(
    tot_symp_mid = sum(pre_vacc),
    tot_symp_lo  = sum(pre_vacc_lo),
    tot_symp_hi  = sum(pre_vacc_hi),
    .groups = "drop"
  ) %>%
  mutate(
    true_inf_mid = tot_symp_mid / p_symp,
    true_inf_lo  = tot_symp_lo  / p_symp,
    true_inf_hi  = tot_symp_hi  / p_symp
  )


## updates (using all draws)
true_inf_region <- all_draws_ix_true %>%
  group_by(Region, draw_id) %>%
  summarise(
    tot_symp = sum(total_pre, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Region) %>%
  summarise(
    tot_symp_mid = median(tot_symp, na.rm = TRUE),
    tot_symp_lo  = quantile(tot_symp, probs = 0.025, na.rm = TRUE),
    tot_symp_hi  = quantile(tot_symp, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    true_inf_mid = tot_symp_mid / p_symp,
    true_inf_lo  = tot_symp_lo  / p_symp,
    true_inf_hi  = tot_symp_hi  / p_symp
  ) %>% dplyr::rename(region = Region)


ar_target_region <- true_inf_region %>%
  left_join(pop_by_state, by = "region") %>%    # tot_pop 
  left_join(S0_df,    by = "region") %>%    # S0_pop 
  mutate(
   
    AR_mid_totpop = 1 - exp(- true_inf_mid / tot_pop),
    AR_lo_totpop  = 1 - exp(- true_inf_lo  / tot_pop),
    AR_hi_totpop  = 1 - exp(- true_inf_hi  / tot_pop),
    
    AR_mid_S0 = true_inf_mid / total_S0,
    AR_lo_S0  = true_inf_lo  / total_S0,
    AR_hi_S0  = true_inf_hi  / total_S0,
    
    S_prop    = total_S0 / tot_pop
  )

AR_draw_by_region <- all_draws_ix_true %>%
  group_by(Region, draw_id) %>%
  summarise(tot_symp = sum(total_pre, na.rm = TRUE), .groups="drop") %>%
  left_join(S0_df, by = c("Region" = "region")) %>%
  mutate(
    tot_inf = tot_symp / p_symp,
    AR_S0   = tot_inf / total_S0
  )

AR_summary_by_region <- AR_draw_by_region %>%
  filter(!is.na(AR_S0), AR_S0 > 0, AR_S0 < 1) %>%   # support 밖 제거
  group_by(Region) %>%
  summarise(
    AR_mid_S0 = median(AR_S0),
    AR_lo_S0  = quantile(AR_S0, 0.025),
    AR_hi_S0  = quantile(AR_S0, 0.975),
    .groups="drop"
  ) %>%
  dplyr::rename(region = Region)

AR_by_state <- list(
  mid = setNames(AR_summary_by_region$AR_mid_S0, AR_summary_by_region$region),
  lo  = setNames(AR_summary_by_region$AR_lo_S0,  AR_summary_by_region$region),
  hi  = setNames(AR_summary_by_region$AR_hi_S0,  AR_summary_by_region$region)
)

scale_foi_to_target_ar <- function(foi_daily, ar_target) {
  H0 <- sum(foi_daily, na.rm = TRUE)
  Ht <- -log(1 - ar_target)
  m  <- Ht / H0
  foi_daily * m
}

foi_daily_adj_by_state <- lapply(names(foi_daily_by_state), function(st) {
  foi0 <- foi_daily_by_state[[st]]
  list(
    mid = scale_foi_to_target_ar(foi0, AR_by_state$mid[[st]]),
    lo  = scale_foi_to_target_ar(foi0, AR_by_state$lo[[st]]),
    hi  = scale_foi_to_target_ar(foi0, AR_by_state$hi[[st]])
  )
})
names(foi_daily_adj_by_state) <- names(foi_daily_by_state)

foi_daily_by_state_mid <- lapply(foi_daily_adj_by_state, `[[`, "mid")

states_to_run <- names(foi_daily_by_state_mid)

## end ------------------------------------------------------------------------