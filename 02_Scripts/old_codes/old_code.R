compute_outcome <- function(AR, p_nat, p_vacc, VE) {
  
  # Infection related risk (per 10,000)
  risk_inf_nv_10k <- 1e4 * AR * p_nat
  risk_inf_v_10k  <- 1e4 * AR * p_nat * (1 - VE)
  
  # Vaccine related AE risk (per 10,000)
  risk_vacc_AE_10k <- 1e4 * p_vacc
  
  # 1. Risk in unvaccinated individuals (per 10,000)
  risk_nv_10k <- risk_inf_nv_10k
  
  # 2. Risk in vaccinated individuals (per 10,000)
  risk_v_10k  <- risk_inf_v_10k + risk_vacc_AE_10k
  
  # 3. Averted outcomes (infection prevented due to vaccine)
  averted_10k <- risk_inf_nv_10k - risk_inf_v_10k   
  
  # 4. Excess risk caused by vaccine
  excess_10k  <- risk_vacc_AE_10k                   
  
  # 5. Benefit-Risk Ratio
  brr         <- ifelse(excess_10k == 0, NA, averted_10k / excess_10k)
  
  # 6. Net benefit (difference in cases prevented vs caused)
  net_10k     <- averted_10k - excess_10k
  
  list(
    risk_nv_10k = risk_nv_10k,
    risk_v_10k  = risk_v_10k,
    averted_10k = averted_10k,
    excess_10k  = excess_10k,
    brr         = brr,
    net_10k     = net_10k
  )
}



# # updated: four age + daly
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
  p_sae_nat_59 <- lhs_sample$p_sae_nat_59[d]
  p_sae_nat_60 <- lhs_sample$p_sae_nat_60[d]
  
  # Natural death (symptomatic -> death) by age group
  p_death_nat_11 <- lhs_sample$p_death_nat_11[d]
  p_death_nat_17 <- lhs_sample$p_death_nat_17[d]
  p_death_nat_59 <- lhs_sample$p_death_nat_59[d]
  p_death_nat_60 <- lhs_sample$p_death_nat_60[d]
  
  # Other parameters
  epi_months <- lhs_sample$epi_months[d]
  trav_7d    <- lhs_sample$trav_7d[d]
  trav_14d   <- lhs_sample$trav_14d[d]
  trav_30d   <- lhs_sample$trav_30d[d]
  trav_90d   <- lhs_sample$trav_90d[d]
  
  travel_days <- list(
    "7d"  = trav_7d,
    "14d" = trav_14d,
    "30d" = trav_30d,
    "90d" = trav_90d
  )
  
  # Convert epidemic duration in months to days
  L_days <- max(1L, round(epi_months) * 30.5)
  
  # Choose symptomatic proportion for this draw 
  symp_prop_d <- as.numeric(draw_pars$symp_overall)
  
  # Loop through age groups
  for (i in seq_len(nrow(all_risk))) {
    group <- all_risk[i, ]
    age   <- as.character(group$age_group)
    
    # 1) Assign vaccine-related risks:
    #    1–11, 12–17, 18–59 -> u65 values
    #    60 -> 65+ values
    if (age == "60") {
      p_sae_vacc   <- p_sae_vacc_65
      p_death_vacc <- p_death_vacc_65
    } else {
      p_sae_vacc   <- p_sae_vacc_u65
      p_death_vacc <- p_death_vacc_u65
    }
    
    # 2) Assign natural risks for each specific age group
    if (age == "1-11") {
      p_sae_nat   <- p_sae_nat_11
      p_death_nat <- p_death_nat_11
    } else if (age == "12-17") {
      p_sae_nat   <- p_sae_nat_17
      p_death_nat <- p_death_nat_17
    } else if (age == "18-59") {
      p_sae_nat   <- p_sae_nat_59
      p_death_nat <- p_death_nat_59
    } else if (age == "60") {
      p_sae_nat   <- p_sae_nat_60
      p_death_nat <- p_death_nat_60
    } else {
      stop("Unknown age group in all_risk: ", age)
    }
    
    # Loop through continuous AR grid
    for (AR_total in ar_grid) {
      
      # Compute daily FOI for this specific AR
      foi_daily <- compute_daily_foi(AR_total = AR_total, L = L_days)
      
      # Loop through travel durations
      for (days_label in names(travel_days)) {
        
        D <- round(travel_days[[days_label]])
        max_entry <- max(1L, L_days - D + 1L)
        
        # Sample multiple entry days to average over entry timing uncertainty
        entry_days <- sample.int(max_entry, size = n_entry_samples, replace = TRUE)
        
        # Storage for entry-day results
        temp_results <- vector("list", length(entry_days))
        
        for (entry_idx in seq_along(entry_days)) {
          entry_day <- entry_days[entry_idx]
          
          # Travel-specific infection probability (AR_travel)
          AR_travel <- compute_ar_travel(foi_daily, entry_day, D)
          
          # infection per 10k (same for nv & v because VE is disease prevention)
          inf_10k <- 1e4 * AR_travel
          
          # symptomatic per 10k
          symp_nv_10k <- inf_10k * symp_prop_d
          symp_v_10k  <- symp_nv_10k * (1 - ve_d)
          
          # natural hosp/death per 10k (conditional on symptomatic)
          hosp_nv_10k  <- symp_nv_10k * p_sae_nat
          hosp_v_10k   <- symp_v_10k  * p_sae_nat
          
          death_nv_10k <- symp_nv_10k * p_death_nat
          death_v_10k  <- symp_v_10k  * p_death_nat
          
          nonhosp_symp_nv_10k <- pmax(0, symp_nv_10k - hosp_nv_10k)
          nonhosp_symp_v_10k  <- pmax(0, symp_v_10k  - hosp_v_10k)
          
          # vaccine SAE/death per 10k
          sae_vacc_10k   <- 1e4 * p_sae_vacc
          death_vacc_10k <- 1e4 * p_death_vacc
          
          # DALY: disease nv / disease v / vaccine SAE
          dz_nv <- compute_daly_one(
            age_group = age,
            deaths_10k = death_nv_10k,
            hosp_10k = hosp_nv_10k,
            nonhosp_symp_10k = nonhosp_symp_nv_10k,
            symp_10k = symp_nv_10k,
            sae_10k = 0, deaths_sae_10k = 0,
            draw_pars = draw_pars
          )
          
          dz_v <- compute_daly_one(
            age_group = age,
            deaths_10k = death_v_10k,
            hosp_10k = hosp_v_10k,
            nonhosp_symp_10k = nonhosp_symp_v_10k,
            symp_10k = symp_v_10k,
            sae_10k = 0, deaths_sae_10k = 0,
            draw_pars = draw_pars
          )
          
          sae <- compute_daly_one(
            age_group = age,
            deaths_10k = 0, hosp_10k = 0, nonhosp_symp_10k = 0, symp_10k = 0,
            sae_10k = sae_vacc_10k,
            deaths_sae_10k = death_vacc_10k,
            draw_pars = draw_pars
          )
          
          daly_averted <- dz_nv$daly_dz - dz_v$daly_dz
          daly_sae     <- sae$daly_sae
          brr_daly     <- ifelse(daly_sae > 0, daly_averted / daly_sae, NA_real_)
          
          averted_10k_sae <- hosp_nv_10k - hosp_v_10k
          excess_10k_sae  <- sae_vacc_10k
          brr_sae         <- ifelse(excess_10k_sae > 0, averted_10k_sae / excess_10k_sae, NA_real_)
          net_10k_sae     <- averted_10k_sae - excess_10k_sae
          
          averted_10k_death <- death_nv_10k - death_v_10k
          excess_10k_death  <- death_vacc_10k
          brr_death         <- ifelse(excess_10k_death > 0, averted_10k_death / excess_10k_death, NA_real_)
          net_10k_death     <- averted_10k_death - excess_10k_death
          
          temp_results[[entry_idx]] <- list(
            AR_travel = AR_travel,
            # event metrics
            risk_nv_10k_sae = hosp_nv_10k,
            risk_v_10k_sae  = hosp_v_10k + sae_vacc_10k,
            averted_10k_sae = averted_10k_sae,
            excess_10k_sae  = excess_10k_sae,
            brr_sae = brr_sae,
            net_10k_sae = net_10k_sae,
            
            risk_nv_10k_death = death_nv_10k,
            risk_v_10k_death  = death_v_10k + death_vacc_10k,
            averted_10k_death = averted_10k_death,
            excess_10k_death  = excess_10k_death,
            brr_death = brr_death,
            net_10k_death = net_10k_death,
            
            # DALY metrics
            daly_dz_nv = dz_nv$daly_dz,
            daly_dz_v  = dz_v$daly_dz,
            daly_averted = daly_averted,
            daly_sae = daly_sae,
            brr_daly = brr_daly,
            yll_dz_nv = dz_nv$yll_dz,
            yll_dz_v  = dz_v$yll_dz,
            yld_dz_nv = dz_nv$yld_dz,
            yld_dz_v  = dz_v$yld_dz
          )
        }
        
        # ---- Aggregate over entry days (median) ----
        AR_travel_vec <- sapply(temp_results, `[[`, "AR_travel")
        
        risk_nv_10k_sae_vec <- sapply(temp_results, `[[`, "risk_nv_10k_sae")
        risk_v_10k_sae_vec  <- sapply(temp_results, `[[`, "risk_v_10k_sae")
        averted_10k_sae_vec <- sapply(temp_results, `[[`, "averted_10k_sae")
        excess_10k_sae_vec  <- sapply(temp_results, `[[`, "excess_10k_sae")
        brr_sae_vec         <- sapply(temp_results, `[[`, "brr_sae")
        net_10k_sae_vec     <- sapply(temp_results, `[[`, "net_10k_sae")
        
        risk_nv_10k_death_vec <- sapply(temp_results, `[[`, "risk_nv_10k_death")
        risk_v_10k_death_vec  <- sapply(temp_results, `[[`, "risk_v_10k_death")
        averted_10k_death_vec <- sapply(temp_results, `[[`, "averted_10k_death")
        excess_10k_death_vec  <- sapply(temp_results, `[[`, "excess_10k_death")
        brr_death_vec         <- sapply(temp_results, `[[`, "brr_death")
        net_10k_death_vec     <- sapply(temp_results, `[[`, "net_10k_death")
        
        daly_dz_nv_vec <- sapply(temp_results, `[[`, "daly_dz_nv")
        daly_dz_v_vec  <- sapply(temp_results, `[[`, "daly_dz_v")
        daly_averted_vec <- sapply(temp_results, `[[`, "daly_averted")
        daly_sae_vec     <- sapply(temp_results, `[[`, "daly_sae")
        brr_daly_vec     <- sapply(temp_results, `[[`, "brr_daly")
        
        yll_dz_nv_vec <- sapply(temp_results, `[[`, "yll_dz_nv")
        yll_dz_v_vec  <- sapply(temp_results, `[[`, "yll_dz_v")
        yld_dz_nv_vec <- sapply(temp_results, `[[`, "yld_dz_nv")
        yld_dz_v_vec  <- sapply(temp_results, `[[`, "yld_dz_v")
        
        # Store one averaged result per combination
        psa_out_list[[length(psa_out_list) + 1]] <- data.frame(
          draw              = d,
          age_group         = age,
          AR_total          = AR_total,
          AR_total_pct      = AR_total * 100,
          days              = days_label,
          entry_day         = median(entry_days),
          AR                = median(AR_travel_vec, na.rm = TRUE),
          
          risk_nv_10k_sae   = median(risk_nv_10k_sae_vec,   na.rm = TRUE),
          risk_v_10k_sae    = median(risk_v_10k_sae_vec,    na.rm = TRUE),
          averted_10k_sae   = median(averted_10k_sae_vec,   na.rm = TRUE),
          excess_10k_sae    = median(excess_10k_sae_vec,    na.rm = TRUE),
          brr_sae           = median(brr_sae_vec,           na.rm = TRUE),
          
          risk_nv_10k_death = median(risk_nv_10k_death_vec, na.rm = TRUE),
          risk_v_10k_death  = median(risk_v_10k_death_vec,  na.rm = TRUE),
          averted_10k_death = median(averted_10k_death_vec, na.rm = TRUE),
          excess_10k_death  = median(excess_10k_death_vec,  na.rm = TRUE),
          brr_death         = median(brr_death_vec,         na.rm = TRUE),
          
          # DALY metrics
          daly_nv      = median(daly_dz_nv_vec,  na.rm = TRUE),
          daly_v       = median(daly_dz_v_vec,   na.rm = TRUE),
          daly_averted = median(daly_averted_vec, na.rm = TRUE),
          daly_sae     = median(daly_sae_vec,     na.rm = TRUE),
          brr_daly     = median(brr_daly_vec,     na.rm = TRUE),
          
          yll_nv = median(yll_dz_nv_vec, na.rm = TRUE),
          yll_v  = median(yll_dz_v_vec,  na.rm = TRUE),
          yld_nv = median(yld_dz_nv_vec, na.rm = TRUE),
          yld_v  = median(yld_dz_v_vec,  na.rm = TRUE)
        )
      }
    }
  }
  
  if (d %% 10 == 0) {
    cat("Completed", d, "of", nrow(lhs_sample), "PSA draws\n")
  }
}

# no dalys
for (d in seq_len(nrow(lhs_sample))) { 
  
  # Sample uncertain parameters from LHS
  draw_pars <- lhs_sample[d, ] 
  ve_d      <- lhs_sample$ve[d]
  
  p_sae_vacc_u65    <- lhs_sample$p_sae_vacc_u65[d]
  p_sae_vacc_65     <- lhs_sample$p_sae_vacc_65[d]
  p_death_vacc_u65  <- lhs_sample$p_death_vacc_u65[d]
  p_death_vacc_65   <- lhs_sample$p_death_vacc_65[d]
  
  p_sae_nat_u65     <- lhs_sample$p_sae_nat_u65[d]
  p_sae_nat_65      <- lhs_sample$p_sae_nat_65[d]
  p_death_nat_u65   <- lhs_sample$p_death_nat_u65[d]
  p_death_nat_65    <- lhs_sample$p_death_nat_65[d]
  
  epi_months        <- lhs_sample$epi_months[d]
  trav_7d           <- lhs_sample$trav_7d[d]
  trav_14d          <- lhs_sample$trav_14d[d]
  trav_30d          <- lhs_sample$trav_30d[d]
  trav_90d          <- lhs_sample$trav_90d[d]
  
  travel_days <- list(
    "7d"  = trav_7d,
    "14d" = trav_14d,
    "30d" = trav_30d,
    "90d" = trav_90d
  )
  
  L_days <- max(1L, round(epi_months) * 30.5)
  
  # Loop through age groups
  for (i in seq_len(nrow(all_risk))) {
    group <- all_risk[i, ]
    age   <- group$age_group
    
    if (grepl("65", age) && !grepl("under", age)) {
      p_sae_vacc   <- p_sae_vacc_65
      p_death_vacc <- p_death_vacc_65
      p_sae_nat    <- p_sae_nat_65
      p_death_nat  <- p_death_nat_65
    } else {
      p_sae_vacc   <- p_sae_vacc_u65
      p_death_vacc <- p_death_vacc_u65
      p_sae_nat    <- p_sae_nat_u65
      p_death_nat  <- p_death_nat_u65
    }
    
    # Loop through continuous AR grid
    for (AR_total in ar_grid) {
      
      # Compute daily FOI for this specific AR
      foi_daily <- compute_daily_foi(AR_total = AR_total, L = L_days)
      
      # Loop through travel durations
      for (days_label in names(travel_days)) {
        
        D <- round(travel_days[[days_label]])
        max_entry <- max(1L, L_days - D + 1L)
        
        # Sample multiple entry days to average over entry timing uncertainty
        entry_days <- sample.int(max_entry, size = n_entry_samples, replace = TRUE)
        
        # Initialize storage for this combination
        temp_results <- list()
        
        # Compute outcomes for each entry day
        for (entry_idx in seq_along(entry_days)) {
          entry_day <- entry_days[entry_idx]
          
          # Travel-specific AR for this entry day
          AR_travel <- compute_ar_travel(foi_daily, entry_day, D)
          
          # Compute outcomes (SAE, Death) only
          sae <- compute_outcome(
            AR_travel, 
            p_nat = p_sae_nat, 
            p_vacc = p_sae_vacc, 
            VE = ve_d
          )
          
          death <- compute_outcome(
            AR_travel, 
            p_nat = p_death_nat, 
            p_vacc = p_death_vacc, 
            VE = ve_d
          )
          
          # Store intermediate results (DALY removed)
          temp_results[[entry_idx]] <- list(
            entry_day  = entry_day,
            AR_travel  = AR_travel,
            sae        = sae,
            death      = death
          )
        }
        
        # Median across entry day samples
        AR_travel_vec <- sapply(temp_results, function(x) x$AR_travel)
        
        risk_nv_10k_sae_vec   <- sapply(temp_results, function(x) x$sae$risk_nv_10k)
        risk_v_10k_sae_vec    <- sapply(temp_results, function(x) x$sae$risk_v_10k)
        averted_10k_sae_vec   <- sapply(temp_results, function(x) x$sae$averted_10k)
        excess_10k_sae_vec    <- sapply(temp_results, function(x) x$sae$excess_10k)
        brr_sae_vec           <- sapply(temp_results, function(x) x$sae$brr)
        net_10k_sae_vec       <- sapply(temp_results, function(x) x$sae$net_10k)
        
        risk_nv_10k_death_vec <- sapply(temp_results, function(x) x$death$risk_nv_10k)
        risk_v_10k_death_vec  <- sapply(temp_results, function(x) x$death$risk_v_10k)
        averted_10k_death_vec <- sapply(temp_results, function(x) x$death$averted_10k)
        excess_10k_death_vec  <- sapply(temp_results, function(x) x$death$excess_10k)
        brr_death_vec         <- sapply(temp_results, function(x) x$death$brr)
        net_10k_death_vec     <- sapply(temp_results, function(x) x$death$net_10k)
        
        # Store single averaged result per combination (DALY columns removed)
        psa_out_list[[length(psa_out_list) + 1]] <- data.frame(
          draw              = d,
          age_group         = age,
          AR_total          = AR_total,
          AR_total_pct      = AR_total * 100,
          days              = days_label,
          entry_day         = median(entry_days),
          AR                = median(AR_travel_vec, na.rm = TRUE),
          
          risk_nv_10k_sae   = median(risk_nv_10k_sae_vec, na.rm = TRUE),
          risk_v_10k_sae    = median(risk_v_10k_sae_vec, na.rm = TRUE),
          averted_10k_sae   = median(averted_10k_sae_vec, na.rm = TRUE),
          excess_10k_sae    = median(excess_10k_sae_vec, na.rm = TRUE),
          brr_sae           = median(brr_sae_vec, na.rm = TRUE),
          net_10k_sae       = median(net_10k_sae_vec, na.rm = TRUE),
          
          risk_nv_10k_death = median(risk_nv_10k_death_vec, na.rm = TRUE),
          risk_v_10k_death  = median(risk_v_10k_death_vec, na.rm = TRUE),
          averted_10k_death = median(averted_10k_death_vec, na.rm = TRUE),
          excess_10k_death  = median(excess_10k_death_vec, na.rm = TRUE),
          brr_death         = median(brr_death_vec, na.rm = TRUE),
          net_10k_death     = median(net_10k_death_vec, na.rm = TRUE)
        )
      }
    }
  }
  
  # Progress indicator (optional)
  if (d %% 10 == 0) {
    cat("Completed", d, "of", nrow(lhs_sample), "PSA draws\n")
  }
}

# four age groups
for (d in seq_len(nrow(lhs_sample))) { 
  
  # Sample uncertain parameters from LHS
  draw_pars <- lhs_sample[d, ] 
  ve_d      <- lhs_sample$ve[d]
  
  # Vaccine-related SAE / death (u65 vs 65+)
  p_sae_vacc_u65   <- lhs_sample$p_sae_vacc_u65[d]
  p_sae_vacc_65    <- lhs_sample$p_sae_vacc_65[d]
  p_death_vacc_u65 <- lhs_sample$p_death_vacc_u65[d]
  p_death_vacc_65  <- lhs_sample$p_death_vacc_65[d]
  
  # Natural SAE by age group
  p_sae_nat_11 <- lhs_sample$p_sae_nat_11[d]
  p_sae_nat_17 <- lhs_sample$p_sae_nat_17[d]
  p_sae_nat_59 <- lhs_sample$p_sae_nat_59[d]
  p_sae_nat_60 <- lhs_sample$p_sae_nat_60[d]
  
  # Natural death by age group
  p_death_nat_11 <- lhs_sample$p_death_nat_11[d]
  p_death_nat_17 <- lhs_sample$p_death_nat_17[d]
  p_death_nat_59 <- lhs_sample$p_death_nat_59[d]
  p_death_nat_60 <- lhs_sample$p_death_nat_60[d]
  
  # Other parameters
  epi_months        <- lhs_sample$epi_months[d]
  trav_7d           <- lhs_sample$trav_7d[d]
  trav_14d          <- lhs_sample$trav_14d[d]
  trav_30d          <- lhs_sample$trav_30d[d]
  trav_90d          <- lhs_sample$trav_90d[d]
  
  travel_days <- list(
    "7d"  = trav_7d,
    "14d" = trav_14d,
    "30d" = trav_30d,
    "90d" = trav_90d
  )
  
  # Convert epidemic duration in months to days
  L_days <- max(1L, round(epi_months) * 30.5)
  
  # Loop through age groups
  for (i in seq_len(nrow(all_risk))) {
    group <- all_risk[i, ]
    age   <- as.character(group$age_group)
    
    # 1) Assign vaccine-related risks:
    #    1–11, 12–17, 18–59 -> u65 values
    #    60 -> 65+ values
    if (age == "60") {
      p_sae_vacc   <- p_sae_vacc_65
      p_death_vacc <- p_death_vacc_65
    } else {
      p_sae_vacc   <- p_sae_vacc_u65
      p_death_vacc <- p_death_vacc_u65
    }
    
    # 2) Assign natural risks for each specific age group
    if (age == "1-11") {
      p_sae_nat   <- p_sae_nat_11
      p_death_nat <- p_death_nat_11
    } else if (age == "12-17") {
      p_sae_nat   <- p_sae_nat_17
      p_death_nat <- p_death_nat_17
    } else if (age == "18-59") {
      p_sae_nat   <- p_sae_nat_59
      p_death_nat <- p_death_nat_59
    } else if (age == "60") {
      p_sae_nat   <- p_sae_nat_60
      p_death_nat <- p_death_nat_60
    } else {
      stop("Unknown age group in all_risk: ", age)
    }
    
    # Loop through continuous AR grid
    for (AR_total in ar_grid) {
      
      # Compute daily FOI for this specific AR
      foi_daily <- compute_daily_foi(AR_total = AR_total, L = L_days)
      
      # Loop through travel durations
      for (days_label in names(travel_days)) {
        
        D <- round(travel_days[[days_label]])
        max_entry <- max(1L, L_days - D + 1L)
        
        # Sample multiple entry days to average over entry timing uncertainty
        entry_days <- sample.int(max_entry, size = n_entry_samples, replace = TRUE)
        
        # Temporary storage for this combination
        temp_results <- list()
        
        # Compute outcomes for each sampled entry day
        for (entry_idx in seq_along(entry_days)) {
          entry_day <- entry_days[entry_idx]
          
          # Travel-specific AR for this entry day
          AR_travel <- compute_ar_travel(foi_daily, entry_day, D)
          
          # SAE outcomes
          sae <- compute_outcome(
            AR_travel, 
            p_nat  = p_sae_nat, 
            p_vacc = p_sae_vacc, 
            VE     = ve_d
          )
          
          # Death outcomes
          death <- compute_outcome(
            AR_travel, 
            p_nat  = p_death_nat, 
            p_vacc = p_death_vacc, 
            VE     = ve_d
          )
          
          # Store intermediate results (without DALYs)
          temp_results[[entry_idx]] <- list(
            entry_day  = entry_day,
            AR_travel  = AR_travel,
            sae        = sae,
            death      = death
          )
        }
        
        # Take median across entry day samples
        AR_travel_vec <- sapply(temp_results, function(x) x$AR_travel)
        
        risk_nv_10k_sae_vec   <- sapply(temp_results, function(x) x$sae$risk_nv_10k)
        risk_v_10k_sae_vec    <- sapply(temp_results, function(x) x$sae$risk_v_10k)
        averted_10k_sae_vec   <- sapply(temp_results, function(x) x$sae$averted_10k)
        excess_10k_sae_vec    <- sapply(temp_results, function(x) x$sae$excess_10k)
        brr_sae_vec           <- sapply(temp_results, function(x) x$sae$brr)
        net_10k_sae_vec       <- sapply(temp_results, function(x) x$sae$net_10k)
        
        risk_nv_10k_death_vec <- sapply(temp_results, function(x) x$death$risk_nv_10k)
        risk_v_10k_death_vec  <- sapply(temp_results, function(x) x$death$risk_v_10k)
        averted_10k_death_vec <- sapply(temp_results, function(x) x$death$averted_10k)
        excess_10k_death_vec  <- sapply(temp_results, function(x) x$death$excess_10k)
        brr_death_vec         <- sapply(temp_results, function(x) x$death$brr)
        net_10k_death_vec     <- sapply(temp_results, function(x) x$death$net_10k)
        
        # Store one averaged result per combination (DALY columns removed)
        psa_out_list[[length(psa_out_list) + 1]] <- data.frame(
          draw              = d,
          age_group         = age,
          AR_total          = AR_total,
          AR_total_pct      = AR_total * 100,
          days              = days_label,
          entry_day         = median(entry_days),
          AR                = median(AR_travel_vec, na.rm = TRUE),
          
          risk_nv_10k_sae   = median(risk_nv_10k_sae_vec,   na.rm = TRUE),
          risk_v_10k_sae    = median(risk_v_10k_sae_vec,    na.rm = TRUE),
          averted_10k_sae   = median(averted_10k_sae_vec,   na.rm = TRUE),
          excess_10k_sae    = median(excess_10k_sae_vec,    na.rm = TRUE),
          brr_sae           = median(brr_sae_vec,           na.rm = TRUE),
          net_10k_sae       = median(net_10k_sae_vec,       na.rm = TRUE),
          
          risk_nv_10k_death = median(risk_nv_10k_death_vec, na.rm = TRUE),
          risk_v_10k_death  = median(risk_v_10k_death_vec,  na.rm = TRUE),
          averted_10k_death = median(averted_10k_death_vec, na.rm = TRUE),
          excess_10k_death  = median(excess_10k_death_vec,  na.rm = TRUE),
          brr_death         = median(brr_death_vec,         na.rm = TRUE),
          net_10k_death     = median(net_10k_death_vec,     na.rm = TRUE)
        )
      }
    }
  }
  
  # Progress indicator (optional)
  if (d %% 10 == 0) {
    cat("Completed", d, "of", nrow(lhs_sample), "PSA draws\n")
  }
}
