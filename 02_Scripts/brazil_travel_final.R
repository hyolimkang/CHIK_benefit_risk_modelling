## load simulation data for all 11 states in brazil ----------------------------
load("01_Data/sim_results_vc_ixchiq_model.RData")

# extract weekly FOI for each state
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

# this is raw foi -- use this FOI to scale up
foi_daily_by_state <- split(phi_attack$phi, phi_attack$state)
foi_daily_by_state <- lapply(foi_daily_by_state, weekly_to_daily_foi)


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

adjust_foi_using_rho_via_ar <- function(foi_daily, rho) {
  H <- sum(foi_daily, na.rm = TRUE)
  AR_model <- 1 - exp(-H)
  AR_star  <- AR_model / rho
  H_star   <- -log(1 - AR_star)
  m <- H_star / H
  foi_daily * m
}

# make adjusted FOI list (same names as foi_daily_by_state)
foi_daily_adj_by_state <- lapply(names(foi_daily_by_state), function(st) {
  foi0 <- foi_daily_by_state[[st]]
  rho  <- rho_med[[st]]
  adjust_foi_using_rho_via_ar(foi0, rho)
})
names(foi_daily_adj_by_state) <- names(foi_daily_by_state)


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

###### full uncertainty in AR
# ------------------------------------------------------------------------------
# PSA loop (full) with FOI uncertainty levels: lo / mid / hi
# - Assumes you already built:
#   foi_daily_adj_by_state[[state]][["lo"|"mid"|"hi"]]
# - Adds column: ar_level (lo/mid/hi)
# - Reuses the SAME entry_days across lo/mid/hi for fair comparison
# ------------------------------------------------------------------------------

psa_out_list <- list()

for (d in seq_len(nrow(lhs_sample))) {
  
  # ---- Sample uncertain parameters from LHS ----
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
  
  # ---------- loop over states ----------
  for (st in states_to_run) {
    
    # Defensive check
    if (is.null(foi_daily_adj_by_state[[st]])) {
      stop("State not found in foi_daily_adj_by_state: ", st)
    }
    
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
        
        # We need L_days to define entry day range; use mid as reference length
        foi_mid_ref <- foi_daily_adj_by_state[[st]][["mid"]]
        if (is.null(foi_mid_ref)) stop("Missing mid FOI for state: ", st)
        
        L_days <- length(foi_mid_ref)
        max_entry  <- max(1L, L_days - D + 1L)
        
        # ★ Sample entry days ONCE, and reuse across lo/mid/hi
        entry_days <- sample.int(max_entry, size = n_entry_samples, replace = TRUE)
        
        # ---------- loop over AR levels (lo/mid/hi) ----------
        for (ar_level in c("lo", "mid", "hi")) {
          
          foi_daily <- foi_daily_adj_by_state[[st]][[ar_level]]
          if (is.null(foi_daily)) stop("Missing ", ar_level, " FOI for state: ", st)
          
          # In case lengths differ unexpectedly, enforce safety:
          if (length(foi_daily) != L_days) {
            stop("FOI length mismatch for state=", st, " ar_level=", ar_level,
                 " (got ", length(foi_daily), " expected ", L_days, ")")
          }
          
          AR_total_state <- 1 - exp(-sum(foi_daily, na.rm = TRUE))
          
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
              sae_10k = 1e4 * p_sae_vacc,
              deaths_sae_10k = 1e4 * p_death_vacc,
              draw_pars = draw_pars
            )
            
            daly_averted <- dz_nv$daly_dz - dz_v$daly_dz
            daly_sae     <- sae$daly_sae
            brr_daly     <- ifelse(daly_sae > 0, daly_averted / daly_sae, NA_real_)
            
            risk_nv_10k_sae <- res$risk_nv_hosp + res$risk_nv_death
            risk_v_10k_sae  <- (res$risk_v_hosp + res$risk_v_death) + res$excess_10k_sae
            
            temp_results[[entry_idx]] <- list(
              AR_travel = AR_travel,
              
              risk_nv_10k_sae = risk_nv_10k_sae,
              risk_v_10k_sae  = risk_v_10k_sae,
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
            state        = st,
            ar_level     = ar_level,                     # ★ added
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
          
        } # end ar_level loop
        
      } # end days_label loop
    } # end age loop
  } # end state loop
  
  if (d %% 10 == 0) {
    cat("Completed", d, "of", nrow(lhs_sample), "PSA draws\n")
  }
}

# Final combined dataframe
psa_df <- dplyr::bind_rows(psa_out_list)


save(psa_df, file = "01_Data/psa_df_bra_travel_allAR.RData")

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
## three data

psa_data_lo  <- psa_df %>% filter(ar_level == "lo")
psa_data_mid <- psa_df %>% filter(ar_level == "mid")
psa_data_hi  <- psa_df %>% filter(ar_level == "hi")

# br data for mid
br_space_data_benefit_mid <- fn_br_space_benefit(psa_data_mid)

br_space_sd_mid <- fn_br_space_event(psa_data_mid)

br_space_daly_mid <- fn_br_space_daly(psa_data_mid)

# br data for lo
br_space_data_benefit_lo <- fn_br_space_benefit(psa_data_lo)

br_space_sd_lo <- fn_br_space_event(psa_data_lo)

br_space_daly_lo <- fn_br_space_daly(psa_data_lo)

# br data for hi
br_space_data_benefit_hi <- fn_br_space_benefit(psa_data_hi)

br_space_sd_hi <- fn_br_space_event(psa_data_hi)

br_space_daly_hi <- fn_br_space_daly(psa_data_hi)

# combining data
br_representative_benefit_mid <- bind_rows(br_space_sd_mid, br_space_daly_mid)
br_representative_benefit_lo  <- bind_rows(br_space_sd_lo, br_space_daly_lo)
br_representative_benefit_hi  <- bind_rows(br_space_sd_hi, br_space_daly_hi)

# add AR category 
br_representative_benefit_mid <- br_representative_benefit_mid %>%
  mutate(
    age_group = ifelse(age_group == "65", "65+", age_group),
    ar_category = case_when(
      AR_total_pct < 1                    ~ "<1%",
      AR_total_pct >= 1  & AR_total_pct < 10 ~ "1–10%",
      AR_total_pct >= 10   ~ "10%+",
      TRUE                                ~ NA_character_
    ),
    ar_category = factor(ar_category, levels = c("<1%", "1–10%", "10%+"))
  ) %>%
  filter(!is.na(ar_category)) %>%
  dplyr::rename(AgeCat = age_group) 


br_representative_benefit_lo <- br_representative_benefit_lo %>%
  mutate(
    age_group = ifelse(age_group == "65", "65+", age_group),
    ar_category = case_when(
      AR_total_pct < 5                    ~ "<5%",
      TRUE                                ~ NA_character_
    ),
    ar_category = factor(ar_category, levels = c("<5%"))
  ) %>%
  filter(!is.na(ar_category)) %>%
  dplyr::rename(AgeCat = age_group) 

br_representative_benefit_hi <- br_representative_benefit_hi %>%
  mutate(
    age_group = ifelse(age_group == "65", "65+", age_group),
    ar_category = case_when(
      AR_total_pct < 5                    ~ "<5%",
      AR_total_pct >= 5  & AR_total_pct < 20 ~ "5–30%",
      AR_total_pct >= 10   ~ "30%+",
      TRUE                                ~ NA_character_
    ),
    ar_category = factor(ar_category, levels = c("<5%", "5–30%", "30%+"))
  ) %>%
  filter(!is.na(ar_category)) %>%
  dplyr::rename(AgeCat = age_group) 

panel_ranges_mid <- fn_panel_range(br_representative_benefit_mid)
panel_ranges_lo <- fn_panel_range(br_representative_benefit_lo)
panel_ranges_hi <- fn_panel_range(br_representative_benefit_hi)

bg_grid_optimized_mid <- fn_br_grid(panel_ranges_mid)
bg_grid_optimized_lo <- fn_br_grid(panel_ranges_lo)
bg_grid_optimized_hi <- fn_br_grid(panel_ranges_hi)

br_summarized_mid <- fn_br_summ(br_representative_benefit_mid)
br_summarized_lo <- fn_br_summ(br_representative_benefit_lo)
br_summarized_hi <- fn_br_summ(br_representative_benefit_hi)

#### make plots ----------------------------------------------------------------
## -----------------------------------------------------------------------------

log_min <- -2
log_max <- 2
log_range <- seq(log_min, log_max, by = 1)
brr_labels <- c("0.01", "0.1", "1", "10", "100")


### plots mid / 95%UIs

p_daly_mid  <- plot_brr_outcome(br_summarized_mid, bg_grid_optimized_mid,
                                "DALY",  "Benefit-Risk Assessment: DALY",  "#A23B72")

p_death_mid <- plot_brr_outcome(br_summarized_mid, bg_grid_optimized_mid,
                                "Death", "Benefit-Risk Assessment: Death", "#B8860B")

p_sae_mid   <- plot_brr_outcome(br_summarized_mid, bg_grid_optimized_mid,
                                "SAE",   "Benefit-Risk Assessment: SAE",   "#1B7F1B")

p_daly_lo  <- plot_brr_outcome(br_summarized_lo, bg_grid_optimized_lo,
                                "DALY",  "Benefit-Risk Assessment: DALY",  "#A23B72")

p_death_lo <- plot_brr_outcome(br_summarized_lo, bg_grid_optimized_lo,
                                "Death", "Benefit-Risk Assessment: Death", "#B8860B")

p_sae_lo   <- plot_brr_outcome(br_summarized_lo, bg_grid_optimized_lo,
                                "SAE",   "Benefit-Risk Assessment: SAE",   "#1B7F1B")

p_daly_hi  <- plot_brr_outcome(br_summarized_hi, bg_grid_optimized_hi,
                               "DALY",  "Benefit-Risk Assessment: DALY",  "#A23B72")

p_death_hi <- plot_brr_outcome(br_summarized_hi, bg_grid_optimized_hi,
                               "Death", "Benefit-Risk Assessment: Death", "#B8860B")

p_sae_hi   <- plot_brr_outcome(br_summarized_hi, bg_grid_optimized_hi,
                               "SAE",   "Benefit-Risk Assessment: SAE",   "#1B7F1B")

ggsave("06_Results/brr_travel_daly_mid.pdf", plot = p_daly_mid, width = 10, height = 8,
       device = cairo_pdf)
ggsave("06_Results/brr_travel_death_mid.pdf", plot = p_death_mid, width = 10, height = 8)
ggsave("06_Results/brr_travel_sae_mid.pdf", plot = p_sae_mid, width = 10, height = 8)

################################################################################
# table 
################################################################################
ar_summary_all <- psa_data_mid %>%
  pivot_longer(
    cols = c(brr_sae, brr_death, brr_daly),
    names_to = c(".value", "outcome"),
    names_pattern = "(brr|net_10k)_(sae|death|daly)"
  ) %>%
  mutate(
    outcome = dplyr::recode(outcome,
                            "sae"   = "SAE",
                            "death" = "Death",
                            "daly"  = "DALY"),
    days = factor(days, levels = c("7d", "14d", "30d", "90d")),
    age_group = ifelse(age_group == "65", "65+", age_group),
    ar_category = case_when(
      AR_total_pct < 1                       ~ "<1%",
      AR_total_pct >= 1  & AR_total_pct < 10 ~ "1–10%",
      AR_total_pct >= 10                     ~ "10%+",
      TRUE                                   ~ NA_character_
    ),
    ar_category = factor(ar_category, levels = c("<1%", "1–10%", "10%+"))
  ) %>%
  group_by(ar_category, age_group, days, outcome) %>%
  summarise(
    brr_med = median(brr, na.rm = TRUE),
    brr_lo  = quantile(brr, 0.025, na.rm = TRUE),
    brr_hi  = quantile(brr, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

ar_table_wide <- ar_summary_all %>%
  mutate(
    brr_formatted = sprintf("%.2f [%.2f–%.2f]", brr_med, brr_lo, brr_hi)
  ) %>%
  dplyr::select(outcome, ar_category, age_group, days, brr_formatted) %>%
  pivot_wider(
    names_from = days,
    values_from = brr_formatted
  ) %>%
  arrange(outcome, ar_category, age_group)

idx_outcome <- table(ar_table_wide$outcome)

kable(
  ar_table_wide,
  format = "html",
  col.names = c("Outcome", "Attack Rate", "Age Group",
                "7 days", "14 days", "30 days", "90 days"),
  caption = "Benefit–Risk Ratio (BRR) by Outcome, Attack Rate Category, Age Group, and Travel Duration",
  align = c("l", "c", "c", "r", "r", "r", "r")
) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    font_size = 12
  ) %>%
  pack_rows(index = idx_outcome) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2, bold = TRUE)