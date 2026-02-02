rho_global <- median(rho_df$rho_p50, na.rm = TRUE)
rho_global_lo <- median(rho_df$rho_p97.5, na.rm = TRUE) # conservative lower infections
rho_global_hi <- median(rho_df$rho_p2.5,  na.rm = TRUE)


true_inf_region <- combined_nnv_df_region_coverage_model %>%
  distinct(region, pre_inf, pre_inf_lo, pre_inf_hi) %>%   # remove VE/VC duplicates
  group_by(region) %>%
  summarise(
    tot_inf_mid = sum(pre_inf),
    tot_inf_lo  = sum(pre_inf_lo),
    tot_inf_hi  = sum(pre_inf_hi),
    .groups = "drop"
  ) %>%
  mutate(
    true_inf_mid = tot_inf_mid / rho_global,
    true_inf_lo  = tot_inf_lo  / rho_global_lo,  # lower true infections: divide by larger rho
    true_inf_hi  = tot_inf_hi  / rho_global_hi    # upper true infections: divide by smaller rho
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

get_AR_from_sim_out <- function(sim_out) {
  S0 <- sum(sim_out$S[, 1], na.rm = TRUE)
  ST <- sum(sim_out$S[, ncol(sim_out$S)], na.rm = TRUE)
  1 - ST / S0
}


extract_AR_all <- function(sim_results) {
  out <- list(); k <- 1L
  for (st in names(sim_results)) {
    for (ve in names(sim_results[[st]])) {
      for (cv in names(sim_results[[st]][[ve]])) {
        reps <- sim_results[[st]][[ve]][[cv]]
        for (r in seq_along(reps)) {
          sim_out <- reps[[r]]$sim_out
          out[[k]] <- data.frame(
            state = st, VE = ve, cov = cv, rep = r,
            AR = get_AR_from_sim_out(sim_out)
          )
          k <- k + 1L
        }
      }
    }
  }
  do.call(rbind, out)
}

AR_df <- extract_AR_all(sim_results_vc_ixchiq_model)

# make AR* using global burden data
p_symp <- 0.5242478

# If rho_global was calibrated as (reported symptomatic) / (true infections),
# then infection-hazard scaling should use: FOI_true = FOI_model * p_symp / rho_global
foi_daily_true_by_state <- lapply(foi_daily_by_state, function(x) list(
  mid = x * p_symp / rho_global,
  lo  = x * p_symp / rho_global_lo,  # lower risk: larger rho -> smaller FOI
  hi  = x * p_symp / rho_global_hi   # higher risk: smaller rho -> larger FOI
))

AR_state_mid <- sapply(
  foi_daily_true_by_state,
  function(z) 1 - exp(-sum(z$mid, na.rm = TRUE))
)

foi_daily_by_state_mid <- lapply(foi_daily_true_by_state, `[[`, "mid")


