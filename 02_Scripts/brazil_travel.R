## load simulation data for all 11 states in brazil ----------------------------
load("01_Data/sim_results_vc_ixchiq_model.RData")

# extract weekly FOI for each state
extract_phi_long <- function(sim_results) {
  out <- list()
  k <- 1L
  
  states <- names(sim_results)
  for (st in states) {
    ve_names <- names(sim_results[[st]])
    for (ve in ve_names) {
      cov_names <- names(sim_results[[st]][[ve]])
      for (cv in cov_names) {
        reps <- sim_results[[st]][[ve]][[cv]]
        for (r in seq_along(reps)) {
          phi <- reps[[r]]$sim_out$phi
          
          out[[k]] <- data.frame(
            state = st,
            VE    = ve,
            cov   = cv,
            rep   = r,
            week  = seq_along(phi),
            phi   = as.numeric(phi)
          )
          k <- k + 1L
        }
      }
    }
  }
  do.call(rbind, out)
}

phi_df <- extract_phi_long(sim_results_vc_ixchiq_model)

phi_df <- phi_df %>%
  filter(VE == "VE98.9", cov == "cov50", rep == 1)

phi_attack <- phi_df %>% 
  group_by(state, VE, cov, rep) %>%
  mutate(
    FOI_total = sum(phi),
    AR_total  = 1 - exp(-FOI_total)
  ) %>% 
  ungroup()

weekly_to_daily_foi <- function(phi_weekly) {
  phi_weekly <- as.numeric(phi_weekly)
  rep(phi_weekly / 7, each = 7) 
}

foi_daily_by_state <- split(phi_attack$phi, phi_attack$state)
foi_daily_by_state <- lapply(foi_daily_by_state, weekly_to_daily_foi)


states_to_run <- names(foi_daily_by_state)

### psa loop
# ---------- helpers (once, before the PSA loop) ----------
weekly_to_daily_foi <- function(phi_weekly) {
  phi_weekly <- as.numeric(phi_weekly)
  rep(phi_weekly / 7, each = 7)  # 52*7 = 364 days
}

# foi_daily_by_state: named list, each element is daily FOI vector for a state
# assumes you already created phi_df (state, VE, cov, rep, week, phi)

phi_state_weekly <- phi_df %>%
  filter(VE == "VE98.9", cov == "cov50") %>%
  group_by(state, week) %>%
  summarise(phi = median(phi, na.rm = TRUE), .groups = "drop") %>%
  arrange(state, week)

foi_daily_by_state <- split(phi_state_weekly$phi, phi_state_weekly$state)
foi_daily_by_state <- lapply(foi_daily_by_state, weekly_to_daily_foi)

states_to_run <- names(foi_daily_by_state)

rho_med <- setNames(rho_df$rho_p50, rho_df$region)
rho_lo  <- setNames(rho_df$rho_p2.5,rho_df$region)
rho_hi  <- setNames(rho_df$rho_p97.5, rho_df$region)

# ---------- PSA loop ----------
psa_out_list <- list()
n_entry_samples <- 50  # Number of entry day samples to average timing uncertainty

for (d in seq_len(nrow(lhs_sample))) {
  
  # Sample uncertain parameters from LHS
  draw_pars <- lhs_sample[d, ]
  ve_d      <- lhs_sample$ve[d]
  
  # Vaccine-related SAE / death (u65 vs 65+)
  p_sae_vacc_u65   <- lhs_sample$p_sae_vacc_u65[d]
  p_sae_vacc_65    <- lhs_sample$p_sae_vacc_65[d]
  p_death_vacc_u65 <- lhs_sample$p_death_vacc_u65[d]
  p_death_vacc_65  <- lhs_sample$p_death_vacc_65[d]
  
  # Natural hosp (symptomatic -> hosp) by age group
  p_sae_nat_11 <- lhs_sample$p_sae_nat_11[d]  # for 1-11
  p_sae_nat_17 <- lhs_sample$p_sae_nat_17[d]  # for 12-17
  p_sae_nat_64 <- lhs_sample$p_sae_nat_64[d]  # for 18-64
  p_sae_nat_65 <- lhs_sample$p_sae_nat_65[d]  # for 65+
  
  # Natural death (symptomatic -> death) by age group
  p_death_nat_11 <- lhs_sample$p_death_nat_11[d]
  p_death_nat_17 <- lhs_sample$p_death_nat_17[d]
  p_death_nat_64 <- lhs_sample$p_death_nat_64[d]  # for 18-64
  p_death_nat_65 <- lhs_sample$p_death_nat_65[d]  # for 65+
  
  # Travel durations
  travel_days <- list(
    "7d"  = lhs_sample$trav_7d[d],
    "14d" = lhs_sample$trav_14d[d],
    "30d" = lhs_sample$trav_30d[d],
    "90d" = lhs_sample$trav_90d[d]
  )
  
  # Symptomatic proportion for this draw (choose overall; swap if needed)
  symp_prop_d <- as.numeric(draw_pars$symp_overall)
  
  # ---------- loop over states (real FOI) ----------
  for (st in states_to_run) {
    
    foi_daily <- foi_daily_by_state[[st]]
    L_days    <- length(foi_daily)
    
    # state-level total AR implied by FOI (useful to store)
    AR_total_state <- 1 - exp(-sum(foi_daily, na.rm = TRUE))
    
    # ---------- loop over age groups (from all_risk) ----------
    for (i in seq_len(nrow(all_risk))) {
      
      age <- as.character(all_risk$age_group[i])  # "1-11","12-17","18-64","65+"
      
      # 1) Assign vaccine-related risks (u65 vs 65+)
      if (age == "65+") {
        p_sae_vacc   <- p_sae_vacc_65
        p_death_vacc <- p_death_vacc_65
      } else {
        p_sae_vacc   <- p_sae_vacc_u65
        p_death_vacc <- p_death_vacc_u65
      }
      
      # 2) Assign natural risks (age-specific; with proxies)
      if (age == "1-11") {
        p_sae_nat   <- p_sae_nat_11
        p_death_nat <- p_death_nat_11
      } else if (age == "12-17") {
        p_sae_nat   <- p_sae_nat_17
        p_death_nat <- p_death_nat_17
      } else if (age == "18-64") {
        p_sae_nat   <- p_sae_nat_64   
        p_death_nat <- p_death_nat_64
      } else if (age == "65+") {
        p_sae_nat   <- p_sae_nat_65   
        p_death_nat <- p_death_nat_65
      } else {
        stop("Unknown age group in all_risk: ", age)
      }
      
      # ---------- loop over travel durations ----------
      for (days_label in names(travel_days)) {
        
        D <- max(1L, round(travel_days[[days_label]]))
        
        max_entry  <- max(1L, L_days - D + 1L)
        entry_days <- sample.int(max_entry, size = n_entry_samples, replace = TRUE)
        
        temp_results <- vector("list", length(entry_days))
        
        # ---------- loop over sampled entry days ----------
        for (entry_idx in seq_along(entry_days)) {
          
          entry_day <- entry_days[entry_idx]
          
          # Travel-specific infection probability (AR_travel) from real FOI
          AR_travel <- compute_ar_travel(foi_daily, entry_day, D)
          
          # ---- Use compute_outcome (your version) ----
          res <- compute_outcome(
            AR           = AR_travel,
            p_hosp       = symp_prop_d * p_sae_nat,     # (infection -> symp) * (symp -> hosp)
            p_death      = symp_prop_d * p_death_nat,   # (infection -> symp) * (symp -> death)
            p_sae_vacc   = p_sae_vacc,
            p_death_vacc = p_death_vacc,
            VE_hosp      = ve_d,
            VE_death     = ve_d
          )
          
          # symptomatic per 10k (needed for chronic YLD in compute_daly_one)
          symp_nv_10k <- 1e4 * AR_travel * symp_prop_d
          symp_v_10k  <- symp_nv_10k * (1 - ve_d)
          
          # non-hosp symptomatic per 10k
          nonhosp_symp_nv_10k <- pmax(0, symp_nv_10k - res$risk_nv_hosp)
          nonhosp_symp_v_10k  <- pmax(0, symp_v_10k  - res$risk_v_hosp)
          
          # DALY: disease (nv) / disease (v) / vaccine SAE
          dz_nv <- compute_daly_one(
            age_group = age,
            deaths_10k = res$risk_nv_death,
            hosp_10k = res$risk_nv_hosp,
            nonhosp_symp_10k = nonhosp_symp_nv_10k,
            symp_10k = symp_nv_10k,
            sae_10k = 0, deaths_sae_10k = 0,
            draw_pars = draw_pars
          )
          
          dz_v <- compute_daly_one(
            age_group = age,
            deaths_10k = res$risk_v_death,
            hosp_10k = res$risk_v_hosp,
            nonhosp_symp_10k = nonhosp_symp_v_10k,
            symp_10k = symp_v_10k,
            sae_10k = 0, deaths_sae_10k = 0,
            draw_pars = draw_pars
          )
          
          sae <- compute_daly_one(
            age_group = age,
            deaths_10k = 0, hosp_10k = 0, nonhosp_symp_10k = 0, symp_10k = 0,
            sae_10k = res$excess_10k_sae,
            deaths_sae_10k = res$excess_10k_death,
            draw_pars = draw_pars
          )
          
          daly_averted <- dz_nv$daly_dz - dz_v$daly_dz
          daly_sae     <- sae$daly_sae
          brr_daly     <- ifelse(daly_sae > 0, daly_averted / daly_sae, NA_real_)
          
          temp_results[[entry_idx]] <- list(
            AR_travel = AR_travel,
            
            # event metrics (keep same names)
            risk_nv_10k_sae = res$risk_nv_hosp,
            risk_v_10k_sae  = res$risk_v_hosp + res$excess_10k_sae,
            averted_10k_sae = res$averted_10k_sae,
            excess_10k_sae  = res$excess_10k_sae,
            brr_sae         = res$brr_sae,
            
            risk_nv_10k_death = res$risk_nv_death,
            risk_v_10k_death  = res$risk_v_death + res$excess_10k_death,
            averted_10k_death = res$averted_10k_death,
            excess_10k_death  = res$excess_10k_death,
            brr_death         = res$brr_death,
            
            # DALY metrics
            daly_dz_nv   = dz_nv$daly_dz,
            daly_dz_v    = dz_v$daly_dz,
            daly_averted = daly_averted,
            daly_sae     = daly_sae,
            brr_daly     = brr_daly,
            yll_dz_nv    = dz_nv$yll_dz,
            yll_dz_v     = dz_v$yll_dz,
            yld_dz_nv    = dz_nv$yld_dz,
            yld_dz_v     = dz_v$yld_dz
          )
        }
        
        # ---------- aggregate over entry days (median) ----------
        psa_out_list[[length(psa_out_list) + 1]] <- data.frame(
          draw         = d,
          state        = st,
          AR_total     = AR_total_state,
          AR_total_pct = AR_total_state * 100,
          days         = days_label,
          age_group    = age,
          entry_day    = median(entry_days),
          AR           = median(sapply(temp_results, `[[`, "AR_travel"), na.rm = TRUE),
          
          risk_nv_10k_sae = median(sapply(temp_results, `[[`, "risk_nv_10k_sae"), na.rm = TRUE),
          risk_v_10k_sae  = median(sapply(temp_results, `[[`, "risk_v_10k_sae"),  na.rm = TRUE),
          averted_10k_sae = median(sapply(temp_results, `[[`, "averted_10k_sae"), na.rm = TRUE),
          excess_10k_sae  = median(sapply(temp_results, `[[`, "excess_10k_sae"),  na.rm = TRUE),
          brr_sae         = median(sapply(temp_results, `[[`, "brr_sae"),         na.rm = TRUE),
          
          risk_nv_10k_death = median(sapply(temp_results, `[[`, "risk_nv_10k_death"), na.rm = TRUE),
          risk_v_10k_death  = median(sapply(temp_results, `[[`, "risk_v_10k_death"),  na.rm = TRUE),
          averted_10k_death = median(sapply(temp_results, `[[`, "averted_10k_death"), na.rm = TRUE),
          excess_10k_death  = median(sapply(temp_results, `[[`, "excess_10k_death"),  na.rm = TRUE),
          brr_death         = median(sapply(temp_results, `[[`, "brr_death"),         na.rm = TRUE),
          
          daly_nv      = median(sapply(temp_results, `[[`, "daly_dz_nv"),   na.rm = TRUE),
          daly_v       = median(sapply(temp_results, `[[`, "daly_dz_v"),    na.rm = TRUE),
          daly_averted = median(sapply(temp_results, `[[`, "daly_averted"), na.rm = TRUE),
          daly_sae     = median(sapply(temp_results, `[[`, "daly_sae"),     na.rm = TRUE),
          brr_daly     = median(sapply(temp_results, `[[`, "brr_daly"),     na.rm = TRUE),
          
          yll_nv = median(sapply(temp_results, `[[`, "yll_dz_nv"), na.rm = TRUE),
          yll_v  = median(sapply(temp_results, `[[`, "yll_dz_v"),  na.rm = TRUE),
          yld_nv = median(sapply(temp_results, `[[`, "yld_dz_nv"), na.rm = TRUE),
          yld_v  = median(sapply(temp_results, `[[`, "yld_dz_v"),  na.rm = TRUE),
          
          row.names = NULL
        )
      }
    }
  }
  
  if (d %% 10 == 0) {
    cat("Completed", d, "of", nrow(lhs_sample), "PSA draws\n")
  }
}

psa_df <- bind_rows(psa_out_list)

save(psa_df, file = "01_Data/psa_df_bra_travel.RData")

## loop with under-reporting adjsuted FOI--------------------------------------------------
psa_out_list <- list()
n_entry_samples <- 50  # Number of entry day samples to average timing uncertainty

for (d in seq_len(nrow(lhs_sample))) {
  
  # Sample uncertain parameters from LHS
  draw_pars <- lhs_sample[d, ]
  ve_d      <- lhs_sample$ve[d]
  
  # Vaccine-related SAE / death (u65 vs 65+)
  p_sae_vacc_u65   <- lhs_sample$p_sae_vacc_u65[d]
  p_sae_vacc_65    <- lhs_sample$p_sae_vacc_65[d]
  p_death_vacc_u65 <- lhs_sample$p_death_vacc_u65[d]
  p_death_vacc_65  <- lhs_sample$p_death_vacc_65[d]
  
  # Natural hosp (symptomatic -> hosp) by age group
  p_sae_nat_11 <- lhs_sample$p_sae_nat_11[d]
  p_sae_nat_17 <- lhs_sample$p_sae_nat_17[d]
  p_sae_nat_64 <- lhs_sample$p_sae_nat_64[d]
  p_sae_nat_65 <- lhs_sample$p_sae_nat_65[d]
  
  # Natural death (symptomatic -> death) by age group
  p_death_nat_11 <- lhs_sample$p_death_nat_11[d]
  p_death_nat_17 <- lhs_sample$p_death_nat_17[d]
  p_death_nat_64 <- lhs_sample$p_death_nat_64[d]
  p_death_nat_65 <- lhs_sample$p_death_nat_65[d]
  
  # Travel durations
  travel_days <- list(
    "7d"  = lhs_sample$trav_7d[d],
    "14d" = lhs_sample$trav_14d[d],
    "30d" = lhs_sample$trav_30d[d],
    "90d" = lhs_sample$trav_90d[d]
  )
  
  # Symptomatic proportion for this draw
  symp_prop_d <- as.numeric(draw_pars$symp_overall)
  
  # ---------- loop over states (FOI + underreporting adjustment) ----------
  for (st in states_to_run) {
    
    # reported-calibrated daily FOI (state_key로 접근)
    foi_daily_cal <- foi_daily_by_state[[st]]
    
    # state-specific reporting probability (median)
    rho_st <- as.numeric(rho_med[st])
    
    # underreporting-adjusted FOI (hazard scaling)
    foi_daily <- foi_daily_cal / rho_st
    L_days    <- length(foi_daily)
    
    # state-level total AR implied by adjusted FOI
    AR_total_state <- 1 - exp(-sum(foi_daily, na.rm = TRUE))
    
    # ---------- loop over age groups ----------
    for (i in seq_len(nrow(all_risk))) {
      
      age <- as.character(all_risk$age_group[i])
      
      # 1) Assign vaccine-related risks (u65 vs 65+)
      if (age == "65+") {
        p_sae_vacc   <- p_sae_vacc_65
        p_death_vacc <- p_death_vacc_65
      } else {
        p_sae_vacc   <- p_sae_vacc_u65
        p_death_vacc <- p_death_vacc_u65
      }
      
      # 2) Assign natural risks (age-specific; with proxies)
      if (age == "1-11") {
        p_sae_nat   <- p_sae_nat_11
        p_death_nat <- p_death_nat_11
      } else if (age == "12-17") {
        p_sae_nat   <- p_sae_nat_17
        p_death_nat <- p_death_nat_17
      } else if (age == "18-64") {
        p_sae_nat   <- p_sae_nat_64
        p_death_nat <- p_death_nat_64
      } else if (age == "65+") {
        p_sae_nat   <- p_sae_nat_65
        p_death_nat <- p_death_nat_65
      } else {
        stop("Unknown age group in all_risk: ", age)
      }
      
      # ---------- loop over travel durations ----------
      for (days_label in names(travel_days)) {
        
        D <- max(1L, round(travel_days[[days_label]]))
        
        max_entry  <- max(1L, L_days - D + 1L)
        entry_days <- sample.int(max_entry, size = n_entry_samples, replace = TRUE)
        
        temp_results <- vector("list", length(entry_days))
        
        # ---------- loop over sampled entry days ----------
        for (entry_idx in seq_along(entry_days)) {
          
          entry_day <- entry_days[entry_idx]
          
          # Travel-specific infection probability from adjusted FOI
          AR_travel <- compute_ar_travel(foi_daily, entry_day, D)
          
          # ---- Use compute_outcome ----
          res <- compute_outcome(
            AR           = AR_travel,
            p_hosp       = symp_prop_d * p_sae_nat,
            p_death      = symp_prop_d * p_death_nat,
            p_sae_vacc   = p_sae_vacc,
            p_death_vacc = p_death_vacc,
            VE_hosp      = ve_d,
            VE_death     = ve_d
          )
          
          # symptomatic per 10k
          symp_nv_10k <- 1e4 * AR_travel * symp_prop_d
          symp_v_10k  <- symp_nv_10k * (1 - ve_d)
          
          # non-hosp symptomatic per 10k
          nonhosp_symp_nv_10k <- pmax(0, symp_nv_10k - res$risk_nv_hosp)
          nonhosp_symp_v_10k  <- pmax(0, symp_v_10k  - res$risk_v_hosp)
          
          # DALY: disease (nv) / disease (v) / vaccine SAE
          dz_nv <- compute_daly_one(
            age_group = age,
            deaths_10k = res$risk_nv_death,
            hosp_10k = res$risk_nv_hosp,
            nonhosp_symp_10k = nonhosp_symp_nv_10k,
            symp_10k = symp_nv_10k,
            sae_10k = 0, deaths_sae_10k = 0,
            draw_pars = draw_pars
          )
          
          dz_v <- compute_daly_one(
            age_group = age,
            deaths_10k = res$risk_v_death,
            hosp_10k = res$risk_v_hosp,
            nonhosp_symp_10k = nonhosp_symp_v_10k,
            symp_10k = symp_v_10k,
            sae_10k = 0, deaths_sae_10k = 0,
            draw_pars = draw_pars
          )
          
          sae <- compute_daly_one(
            age_group = age,
            deaths_10k = 0, hosp_10k = 0, nonhosp_symp_10k = 0, symp_10k = 0,
            sae_10k = res$excess_10k_sae,
            deaths_sae_10k = res$excess_10k_death,
            draw_pars = draw_pars
          )
          
          daly_averted <- dz_nv$daly_dz - dz_v$daly_dz
          daly_sae     <- sae$daly_sae
          brr_daly     <- ifelse(daly_sae > 0, daly_averted / daly_sae, NA_real_)
          
          temp_results[[entry_idx]] <- list(
            AR_travel = AR_travel,
            
            risk_nv_10k_sae = res$risk_nv_hosp,
            risk_v_10k_sae  = res$risk_v_hosp + res$excess_10k_sae,
            averted_10k_sae = res$averted_10k_sae,
            excess_10k_sae  = res$excess_10k_sae,
            brr_sae         = res$brr_sae,
            
            risk_nv_10k_death = res$risk_nv_death,
            risk_v_10k_death  = res$risk_v_death + res$excess_10k_death,
            averted_10k_death = res$averted_10k_death,
            excess_10k_death  = res$excess_10k_death,
            brr_death         = res$brr_death,
            
            daly_dz_nv   = dz_nv$daly_dz,
            daly_dz_v    = dz_v$daly_dz,
            daly_averted = daly_averted,
            daly_sae     = daly_sae,
            brr_daly     = brr_daly,
            yll_dz_nv    = dz_nv$yll_dz,
            yll_dz_v     = dz_v$yll_dz,
            yld_dz_nv    = dz_nv$yld_dz,
            yld_dz_v     = dz_v$yld_dz
          )
        }
        
        # ---------- aggregate over entry days (median) ----------
        psa_out_list[[length(psa_out_list) + 1]] <- data.frame(
          draw         = d,
          state        = st,  # (정규화된 key로 저장됨; 필요하면 나중에 NAME_1와 동일한 키로 조인)
          AR_total     = AR_total_state,
          AR_total_pct = AR_total_state * 100,
          days         = days_label,
          age_group    = age,
          entry_day    = median(entry_days),
          AR           = median(sapply(temp_results, `[[`, "AR_travel"), na.rm = TRUE),
          
          risk_nv_10k_sae = median(sapply(temp_results, `[[`, "risk_nv_10k_sae"), na.rm = TRUE),
          risk_v_10k_sae  = median(sapply(temp_results, `[[`, "risk_v_10k_sae"),  na.rm = TRUE),
          averted_10k_sae = median(sapply(temp_results, `[[`, "averted_10k_sae"), na.rm = TRUE),
          excess_10k_sae  = median(sapply(temp_results, `[[`, "excess_10k_sae"),  na.rm = TRUE),
          brr_sae         = median(sapply(temp_results, `[[`, "brr_sae"),         na.rm = TRUE),
          
          risk_nv_10k_death = median(sapply(temp_results, `[[`, "risk_nv_10k_death"), na.rm = TRUE),
          risk_v_10k_death  = median(sapply(temp_results, `[[`, "risk_v_10k_death"),  na.rm = TRUE),
          averted_10k_death = median(sapply(temp_results, `[[`, "averted_10k_death"), na.rm = TRUE),
          excess_10k_death  = median(sapply(temp_results, `[[`, "excess_10k_death"),  na.rm = TRUE),
          brr_death         = median(sapply(temp_results, `[[`, "brr_death"),         na.rm = TRUE),
          
          daly_nv      = median(sapply(temp_results, `[[`, "daly_dz_nv"),   na.rm = TRUE),
          daly_v       = median(sapply(temp_results, `[[`, "daly_dz_v"),    na.rm = TRUE),
          daly_averted = median(sapply(temp_results, `[[`, "daly_averted"), na.rm = TRUE),
          daly_sae     = median(sapply(temp_results, `[[`, "daly_sae"),     na.rm = TRUE),
          brr_daly     = median(sapply(temp_results, `[[`, "brr_daly"),     na.rm = TRUE),
          
          yll_nv = median(sapply(temp_results, `[[`, "yll_dz_nv"), na.rm = TRUE),
          yll_v  = median(sapply(temp_results, `[[`, "yll_dz_v"),  na.rm = TRUE),
          yld_nv = median(sapply(temp_results, `[[`, "yld_dz_nv"), na.rm = TRUE),
          yld_v  = median(sapply(temp_results, `[[`, "yld_dz_v"),  na.rm = TRUE),
          
          row.names = NULL
        )
      }
    }
  }
  
  if (d %% 10 == 0) {
    cat("Completed", d, "of", nrow(lhs_sample), "PSA draws\n")
  }
}


psa_df <- bind_rows(psa_out_list)


save(psa_df, file = "01_Data/psa_df_bra_travel.RData")

ar_summary_all <- psa_df %>%
  pivot_longer(
    cols = c(
      brr_sae, brr_death, brr_daly
    ),
    names_to = c(".value", "outcome"),
    names_pattern = "(brr|net_10k)_(sae|death|daly)"
  ) %>%
  mutate(
    outcome = dplyr::recode(
      outcome,
      "sae"  = "SAE",
      "death" = "Death",
      "daly" = "DALY"
    ),
    days = factor(days, levels = c("7d", "14d", "30d", "90d"))
  ) %>%
  group_by(AR_total_pct, age_group, days, outcome, state) %>%
  summarise(
    brr_med = median(brr, na.rm = TRUE),
    brr_lo  = quantile(brr, 0.025, na.rm = TRUE),
    brr_hi  = quantile(brr, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


################################################################################
## brr space
################################################################################
# --- STEP 1: Standardize Data Types ---
br_space_data_benefit <- psa_df %>%
  mutate(days = factor(days, levels = c("7d", "14d", "30d", "90d"))) %>%
  pivot_longer(
    cols = c(
      excess_10k_sae, excess_10k_death,
      averted_10k_sae, averted_10k_death
    ),
    names_to = c(".value", "outcome"),
    names_pattern = "(excess_10k|averted_10k)_(sae|death)"
  ) %>%
  mutate(outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death")) %>%
  group_by(AR_total_pct, age_group, days, outcome) %>%
  summarise(
    # x-axis: vaccine risk
    x_med = median(excess_10k),
    x_lo  = quantile(excess_10k, 0.025),
    x_hi  = quantile(excess_10k, 0.975),
    
    # y-axis: PURE benefit (infection-related outcomes averted)
    y_med = median(averted_10k),
    y_lo  = quantile(averted_10k, 0.025),
    y_hi  = quantile(averted_10k, 0.975),
    
    .groups = "drop"
  )

br_space_sd <- psa_df %>%
  mutate(days = factor(days, levels = c("7d", "14d", "30d", "90d"))) %>%
  pivot_longer(
    cols = c(
      excess_10k_sae, excess_10k_death,
      averted_10k_sae, averted_10k_death
    ),
    names_to = c(".value", "outcome"),
    names_pattern = "(excess_10k|averted_10k)_(sae|death)"
  ) %>%
  mutate(outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death")) %>%
  group_by(AR_total_pct, age_group, days, outcome) %>%
  summarise(
    x_med = median(excess_10k, na.rm = TRUE),
    x_lo = quantile(excess_10k, 0.025, na.rm = TRUE),
    x_hi = quantile(excess_10k, 0.975, na.rm = TRUE),
    
    
    y_med = median(averted_10k, na.rm = TRUE),
    y_lo = quantile(averted_10k, 0.025, na.rm = TRUE),
    y_hi = quantile(averted_10k, 0.975, na.rm = TRUE),
    
    
    .groups = "drop"
  )

# 2) DALY (daly_sae = risk, daly_averted = benefit)
br_space_daly <- psa_df %>%
  mutate(days = factor(days, levels = c("7d", "14d", "30d", "90d")),
         outcome = "DALY") %>%
  group_by(AR_total_pct, age_group, days, outcome) %>%
  summarise(
    x_med = median(daly_sae, na.rm = TRUE),
    x_lo = quantile(daly_sae, 0.025, na.rm = TRUE),
    x_hi = quantile(daly_sae, 0.975, na.rm = TRUE),
    
    
    y_med = median(daly_averted, na.rm = TRUE),
    y_lo = quantile(daly_averted, 0.025, na.rm = TRUE),
    y_hi = quantile(daly_averted, 0.975, na.rm = TRUE),
    
    
    .groups = "drop"
  )

br_representative_benefit <- bind_rows(br_space_sd, br_space_daly)

br_representative_benefit <- br_representative_benefit %>%
  mutate(
    age_group = ifelse(age_group == "65", "65+", age_group),
    ar_category = case_when(
      AR_total_pct < 1                    ~ "<1%",
      AR_total_pct >= 1 & AR_total_pct < 5 ~ "1–5%",
      AR_total_pct >= 5                   ~ "5%+",
      TRUE                                ~ NA_character_
    ),
    ar_category = factor(ar_category, levels = c("<1%", "1–5%", "5%+"))
  ) %>%
  filter(!is.na(ar_category)) %>%
  dplyr::rename(AgeCat = age_group) 

panel_ranges <- br_representative_benefit %>%
  group_by(ar_category, days) %>%
  summarise(
    x_min = min(c(0, x_lo), na.rm = TRUE),  # 0 또는 실제 최소값 중 작은 값
    x_max = max(x_hi, na.rm = TRUE) * 1.05,  # 5% 여유
    y_min = min(c(0, y_lo), na.rm = TRUE),  # 0 또는 실제 최소값 중 작은 값
    y_max = max(y_hi, na.rm = TRUE) * 1.05,  # 5% 여유
    .groups = "drop"
  )

# 3. background
bg_grid_optimized <- panel_ranges %>%
  rowwise() %>%
  do({
    panel_data <- .
    
    x_seq <- seq(panel_data$x_min, panel_data$x_max, length.out = 200)
    y_seq <- seq(panel_data$y_min, panel_data$y_max, length.out = 200)
    
    grid <- expand.grid(x = x_seq, y = y_seq)
    
    grid$brr <- with(grid, ifelse(x > 0, y / x, NA_real_))
    grid$log10_brr <- log10(grid$brr)
    
    # --- 이 부분이 반드시 추가되어야 합니다 ---
    grid$outcome     <- panel_data$outcome      # 지표 이름 (DALY, Death 등)
    grid$ar_category <- panel_data$ar_category  # AR 카테고리
    grid$days        <- panel_data$days         # 여행 기간
    # ---------------------------------------
    
    grid
  }) %>%
  ungroup()


log_min <- -2
log_max <- 2
log_range <- seq(log_min, log_max, by = 1)
brr_labels <- c("0.01", "0.1", "1", "10", "100")


## graph
plot_final <- 
  ggplot() +
  geom_raster(
    data = bg_grid_optimized,
    aes(x = x, y = y, fill = log10_brr),
    interpolate = TRUE
  ) +
  
  scale_fill_gradient2(
    name = "Benefit–risk ratio",
    low = "#ca0020", 
    mid = "#f7f7f7", 
    high = "#0571b0",
    midpoint = 0,
    limits = c(log_min, log_max),
    breaks = log_range,
    labels = brr_labels,
    oob = scales::squish,
    na.value = "white"
  ) +
  
  geom_errorbar(
    data = br_representative_benefit, 
    aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome), 
    width = 0, linewidth = 0.6
  ) +
  geom_errorbarh(
    data = br_representative_benefit, 
    aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome), 
    height = 0, linewidth = 0.6
  ) +
  geom_point(
    data = br_representative_benefit, 
    aes(x = x_med, y = y_med, shape = AgeCat, color = outcome), 
    size = 2.5, stroke = 0.7, fill = "white"
  ) +
  
  # y=x 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "black", linewidth = 0.9, alpha = 0.8) +
  
  # facet_wrap 
  facet_wrap(
    ar_category ~ days, 
    scales = "free", 
    ncol = 4
  ) +
  
  scale_color_manual(
    values = c("SAE" = "#1B7F1B", "Death" = "#B8860B", "DALY" = "#A23B72"), 
    name = "Outcome type"
  ) +
  scale_shape_manual(
    values = c("1-11" = 21, "12-17" = 22, "18-64" = 23, "65+" = 24), 
    name = "Age group"
  ) +
  
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  
  coord_cartesian(clip = "on") +
  
  labs(
    x = "Excess adverse outcomes (per 10,000 vaccinated)",
    y = "Infection-related adverse outcomes averted (per 10,000 vaccinated)",
    title = "Benefit-risk assessment: Traveller vaccination"
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95"),
    panel.spacing = unit(0.8, "lines"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    legend.position = "right"
  )


library(dplyr)

# 11개의 AR 데이터를 카테고리별로 평균 내어 요약

# 1. 지표별/조건별 최적 축 범위 계산
panel_ranges <- br_representative_benefit %>%
  group_by(outcome, ar_category, days) %>% # outcome 반드시 포함
  summarise(
    x_min = 0,
    x_max = max(x_hi, na.rm = TRUE) * 1.1, # 데이터보다 10% 넓게
    y_min = 0,
    y_max = max(y_hi, na.rm = TRUE) * 1.1,
    .groups = "drop"
  )

# 2. 배경 Raster 데이터 생성
bg_grid_optimized <- panel_ranges %>%
  rowwise() %>%
  do({
    panel_data <- .
    
    # 지표별 범위에 맞춘 200x200 그리드
    x_seq <- seq(panel_data$x_min, panel_data$x_max, length.out = 200)
    y_seq <- seq(panel_data$y_min, panel_data$y_max, length.out = 200)
    
    grid <- expand.grid(x = x_seq, y = y_seq)
    
    # BRR 및 Log 변환
    grid$brr <- with(grid, ifelse(x > 0, y / x, NA_real_))
    grid$log10_brr <- log10(grid$brr)
    
    # 필터링과 패싯을 위한 메타데이터 주입
    grid$outcome     <- panel_data$outcome
    grid$ar_category <- panel_data$ar_category
    grid$days        <- panel_data$days
    
    grid
  }) %>%
  ungroup()

br_summarized <- br_representative_benefit %>%
  group_by(outcome, ar_category, days, AgeCat) %>%
  summarise(
    # 중심점 (중앙값 혹은 평균)
    x_med = mean(x_med, na.rm = TRUE),
    y_med = mean(y_med, na.rm = TRUE),
    
    # 에러바 범위 (카테고리 내 데이터들의 평균적인 범위를 사용)
    x_lo = mean(x_lo, na.rm = TRUE),
    x_hi = mean(x_hi, na.rm = TRUE),
    y_lo = mean(y_lo, na.rm = TRUE),
    y_hi = mean(y_hi, na.rm = TRUE),
    
    .groups = "drop"
  )

library(ggplot2)

plot_brr_outcome <- function(target_outcome, title_text, color_val) {
  
  # 1) 요약된 데이터 필터링: .data 대명사를 사용하여 컬럼명 명시
  plot_data <- br_summarized %>% 
    filter(.data$outcome == target_outcome)
  
  # 2) 배경 그리드 필터링
  plot_bg <- bg_grid_optimized %>% 
    filter(.data$outcome == target_outcome)
  
  # [안전 검사] 데이터가 없으면 에러 메시지 출력
  if(nrow(plot_data) == 0) {
    message(paste0("⚠️ Skipping: No data for '", target_outcome, "'"))
    message(paste0("   Available outcomes: ", paste(unique(br_summarized$outcome), collapse=", ")))
    return(NULL)
  }
  
  # 3) 그래프 그리기
  ggplot() +
    # 배경 Raster (해당 Outcome의 축 범위에 최적화)
    geom_raster(
      data = plot_bg,
      aes(x = x, y = y, fill = log10_brr),
      interpolate = TRUE
    ) +
    scale_fill_gradient2(
      name = "Benefit–risk ratio",
      low = "#ca0020", mid = "#f7f7f7", high = "#0571b0",
      midpoint = 0, limits = c(log_min, log_max),
      breaks = log_range, labels = brr_labels,
      oob = scales::squish, na.value = "white"
    ) +
    
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.4) +
    
    # 에러바 및 포인트
    geom_errorbar(data = plot_data, aes(x = x_med, ymin = y_lo, ymax = y_hi), 
                  color = color_val, width = 0, linewidth = 0.6) +
    geom_errorbarh(data = plot_data, aes(y = y_med, xmin = x_lo, xmax = x_hi), 
                   color = color_val, height = 0, linewidth = 0.6) +
    geom_point(data = plot_data, aes(x = x_med, y = y_med, shape = AgeCat), 
               fill = "white", color = color_val, size = 3, stroke = 1) +
    
    facet_wrap(~ ar_category + days, scales = "free", ncol = 4) +
    
    scale_shape_manual(values = c("1-11"=21, "12-17"=22, "18-64"=23, "65+"=24), name = "Age group") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    
    labs(title = title_text, x = "Excess Risk", y = "Averted Outcome") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "gray95"),
      legend.position = "right"
    )
}
# 데이터에 있는 이름이 "DALY"이므로 "p_daly" 대신 "DALY"로 호출
p_daly  <- plot_brr_outcome("DALY",  "Benefit-Risk Assessment: DALY",  "#A23B72")

# Death와 SAE도 데이터셋에 있는 실제 이름(아마 "Death", "SAE")으로 호출하세요
p_death <- plot_brr_outcome("Death", "Benefit-Risk Assessment: Death", "#B8860B")
p_sae   <- plot_brr_outcome("SAE",   "Benefit-Risk Assessment: SAE",   "#1B7F1B")

p_daly
p_death
p_sae

