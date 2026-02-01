
true_inf_region <- combined_nnv_df_region_coverage_model %>%
  distinct(region, pre_inf, pre_inf_lo, pre_inf_hi) %>%   # remove VE/VC duplicates
  group_by(region) %>%
  summarise(
    tot_inf_mid = sum(pre_inf),
    tot_inf_lo  = sum(pre_inf_lo),
    tot_inf_hi  = sum(pre_inf_hi),
    .groups = "drop"
  ) %>%
  left_join(rho_df, by = "region") %>%
  mutate(
    true_inf_mid = tot_inf_mid / rho_p50,
    true_inf_lo  = tot_inf_lo  / rho_p97.5,  # lower true infections: divide by larger rho
    true_inf_hi  = tot_inf_hi  / rho_p2.5    # upper true infections: divide by smaller rho
  )
ar_target_region <- true_inf_region %>%
  left_join(pop_by_state, by = "region") %>%
  mutate(
    AR_mid = 1 - exp(- true_inf_mid / tot_pop),
    AR_lo  = 1 - exp(- true_inf_lo  / tot_pop),
    AR_hi  = 1 - exp(- true_inf_hi  / tot_pop)
  )

AR_by_state <- list(
  mid = setNames(ar_target_region$AR_mid, ar_target_region$region),
  lo  = setNames(ar_target_region$AR_lo,  ar_target_region$region),
  hi  = setNames(ar_target_region$AR_hi,  ar_target_region$region)
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

