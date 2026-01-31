# The clinical attack rates of specific outbreaks ranged from 0.28% to 73.4% in Indonesia in 2001 and Cambodia in 2012
# Define continuous attack rate grid
ar_grid <- seq(0.005, 0.80, by = 0.005)  # 0.5% to 80%, 0.5% increments
# Total: 160 different AR values

# Or more focused range
ar_grid <- seq(0.01, 0.10, by = 0.005)  # 1% to 80%, 0.5% increments
# Total: 79 different AR values

cat("Total AR scenarios:", length(ar_grid), "\n")
cat("AR range:", min(ar_grid) * 100, "% to", max(ar_grid) * 100, "%\n")

###  PSA with continuous AR ----------------------------------------

psa_out_list <- list()

# Note: Remove ar_small, ar_med, ar_large from LHS
# They are now deterministic scenarios, not uncertain parameters

# Number of entry day samples per combination (to average over entry timing uncertainty)
# Higher values = more stable estimates but slower computation
n_entry_samples <- 50  # Sample 50 random entry days per draw/age/AR/duration combination

# updated: four age + daly
# ---------------------------------------------------------
# 2. Main PSA Loop
# ---------------------------------------------------------
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
    
    # 1) Assign vaccine-related risks
    if (age == "65+") {
      p_sae_vacc   <- p_sae_vacc_65
      p_death_vacc <- p_death_vacc_65
    } else {
      p_sae_vacc   <- p_sae_vacc_u65
      p_death_vacc <- p_death_vacc_u65
    }
    
    # 2) Assign natural risks
    if (age == "1-11") {
      p_sae_nat   <- p_sae_nat_11
      p_death_nat <- p_death_nat_11
    } else if (age == "12-17") {
      p_sae_nat   <- p_sae_nat_17
      p_death_nat <- p_death_nat_17
    } else if (age == "18-64") {
      p_sae_nat   <- p_sae_nat_59
      p_death_nat <- p_death_nat_59
    } else if (age == "65+") {
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
          
          # ---- Use compute_outcome  ----
          res <- compute_outcome(
            AR           = AR_travel,
            p_hosp       = symp_prop_d * p_sae_nat,     # symptomatic -> hosp
            p_death      = symp_prop_d * p_death_nat,   # symptomatic -> death
            p_sae_vacc   = p_sae_vacc,
            p_death_vacc = p_death_vacc,
            VE_hosp      = ve_d,
            VE_death     = ve_d
          )
          
          # ---- DALY inputs from res ----
          # symptomatic per 10k
          symp_nv_10k <- 1e4 * AR_travel * symp_prop_d
          symp_v_10k  <- symp_nv_10k * (1 - ve_d)
          
          # non-hosp symptomatic per 10k
          nonhosp_symp_nv_10k <- pmax(0, symp_nv_10k - res$risk_nv_hosp)
          nonhosp_symp_v_10k  <- pmax(0, symp_v_10k  - res$risk_v_hosp)
          
          # DALY: disease nv / disease v / vaccine SAE
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
          
          # ---- Benefit-risk event metrics (match your original outputs) ----
          averted_10k_sae <- res$averted_10k_sae
          excess_10k_sae  <- res$excess_10k_sae
          brr_sae         <- res$brr_sae
          net_10k_sae     <- averted_10k_sae - excess_10k_sae
          
          averted_10k_death <- res$averted_10k_death
          excess_10k_death  <- res$excess_10k_death
          brr_death         <- res$brr_death
          net_10k_death     <- averted_10k_death - excess_10k_death
          
          temp_results[[entry_idx]] <- list(
            AR_travel = AR_travel,
            
            # event metrics (keep same names)
            risk_nv_10k_sae = res$risk_nv_hosp,
            risk_v_10k_sae  = res$risk_v_hosp + res$excess_10k_sae,
            averted_10k_sae = averted_10k_sae,
            excess_10k_sae  = excess_10k_sae,
            brr_sae         = brr_sae,
            net_10k_sae     = net_10k_sae,
            
            risk_nv_10k_death = res$risk_nv_death,
            risk_v_10k_death  = res$risk_v_death + res$excess_10k_death,
            averted_10k_death = averted_10k_death,
            excess_10k_death  = excess_10k_death,
            brr_death         = brr_death,
            net_10k_death     = net_10k_death,
            
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
        
        # ---- Aggregate over entry days (median) ----
        AR_travel_vec <- sapply(temp_results, `[[`, "AR_travel")
        
        risk_nv_10k_sae_vec <- sapply(temp_results, `[[`, "risk_nv_10k_sae")
        risk_v_10k_sae_vec  <- sapply(temp_results, `[[`, "risk_v_10k_sae")
        averted_10k_sae_vec <- sapply(temp_results, `[[`, "averted_10k_sae")
        excess_10k_sae_vec  <- sapply(temp_results, `[[`, "excess_10k_sae")
        brr_sae_vec         <- sapply(temp_results, `[[`, "brr_sae")
        
        risk_nv_10k_death_vec <- sapply(temp_results, `[[`, "risk_nv_10k_death")
        risk_v_10k_death_vec  <- sapply(temp_results, `[[`, "risk_v_10k_death")
        averted_10k_death_vec <- sapply(temp_results, `[[`, "averted_10k_death")
        excess_10k_death_vec  <- sapply(temp_results, `[[`, "excess_10k_death")
        brr_death_vec         <- sapply(temp_results, `[[`, "brr_death")
        
        daly_dz_nv_vec    <- sapply(temp_results, `[[`, "daly_dz_nv")
        daly_dz_v_vec     <- sapply(temp_results, `[[`, "daly_dz_v")
        daly_averted_vec  <- sapply(temp_results, `[[`, "daly_averted")
        daly_sae_vec      <- sapply(temp_results, `[[`, "daly_sae")
        brr_daly_vec      <- sapply(temp_results, `[[`, "brr_daly")
        
        yll_dz_nv_vec <- sapply(temp_results, `[[`, "yll_dz_nv")
        yll_dz_v_vec  <- sapply(temp_results, `[[`, "yll_dz_v")
        yld_dz_nv_vec <- sapply(temp_results, `[[`, "yld_dz_nv")
        yld_dz_v_vec  <- sapply(temp_results, `[[`, "yld_dz_v")
        
        psa_out_list[[length(psa_out_list) + 1]] <- data.frame(
          draw         = d,
          age_group    = age,
          AR_total     = AR_total,
          AR_total_pct = AR_total * 100,
          days         = days_label,
          entry_day    = median(entry_days),
          AR           = median(AR_travel_vec, na.rm = TRUE),
          
          risk_nv_10k_sae = median(risk_nv_10k_sae_vec, na.rm = TRUE),
          risk_v_10k_sae  = median(risk_v_10k_sae_vec,  na.rm = TRUE),
          averted_10k_sae = median(averted_10k_sae_vec, na.rm = TRUE),
          excess_10k_sae  = median(excess_10k_sae_vec,  na.rm = TRUE),
          brr_sae         = median(brr_sae_vec,         na.rm = TRUE),
          
          risk_nv_10k_death = median(risk_nv_10k_death_vec, na.rm = TRUE),
          risk_v_10k_death  = median(risk_v_10k_death_vec,  na.rm = TRUE),
          averted_10k_death = median(averted_10k_death_vec, na.rm = TRUE),
          excess_10k_death  = median(excess_10k_death_vec,  na.rm = TRUE),
          brr_death         = median(brr_death_vec,         na.rm = TRUE),
          
          daly_nv      = median(daly_dz_nv_vec,   na.rm = TRUE),
          daly_v       = median(daly_dz_v_vec,    na.rm = TRUE),
          daly_averted = median(daly_averted_vec, na.rm = TRUE),
          daly_sae     = median(daly_sae_vec,     na.rm = TRUE),
          brr_daly     = median(brr_daly_vec,     na.rm = TRUE),
          
          yll_nv = median(yll_dz_nv_vec, na.rm = TRUE),
          yll_v  = median(yll_dz_v_vec,  na.rm = TRUE),
          yld_nv = median(yld_dz_nv_vec, na.rm = TRUE),
          yld_v  = median(yld_dz_v_vec,  na.rm = TRUE),
          
          row.names = NULL
        )
      }
    }
  }
  
  if (d %% 10 == 0) {
    cat("Completed", d, "of", nrow(lhs_sample), "PSA draws\n")
  }
}

# combine final output
psa_out <- do.call(rbind, psa_out_list)
# Combine results
psa_df <- bind_rows(psa_out_list)

save(psa_df, file = "01_Data/psa_df_conservative_nodaly.RData")
save(psa_df, file = "01_Data/psa_df_four_ages.RData")

cat("\nTotal simulations:", nrow(psa_df), "\n")
cat("Structure:", nrow(lhs_sample), "draws ×", 
    length(ar_grid), "AR values ×",
    nrow(all_risk), "age groups ×",
    length(travel_days), "durations\n")
cat("Note: Each combination averaged over", n_entry_samples, "random entry day samples\n")
cat("Total entry day samples computed:", nrow(psa_df) * n_entry_samples, "\n")

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


# no daly
psa_df <- psa_df |>
  filter(!age_group %in% c("18-64", "65")) |>
  droplevels()

ar_summary_all <- psa_df %>%
  pivot_longer(
    cols = c(
      brr_sae, brr_death,
      net_10k_sae, net_10k_death
    ),
    names_to = c(".value", "outcome"),
    names_pattern = "(brr|net_10k)_(sae|death)"
  ) %>%
  mutate(
    outcome = dplyr::recode(
      outcome,
      "sae"  = "Hospitalisation",
      "death" = "Death"
    ),
    age_group = dplyr::recode(
      age_group,
      "60"    = "60+"
    ),
    days = factor(days, levels = c("7d", "14d", "30d", "90d"))
  ) %>%
  group_by(AR_total_pct, age_group, days, outcome) %>%
  summarise(
    brr_med = median(brr, na.rm = TRUE),
    brr_lo  = quantile(brr, 0.025, na.rm = TRUE),
    brr_hi  = quantile(brr, 0.975, na.rm = TRUE),
    net_med = median(net_10k, na.rm = TRUE),
    net_lo  = quantile(net_10k, 0.025, na.rm = TRUE),
    net_hi  = quantile(net_10k, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Plot with unified Y-axis within each outcome
plot_brr_travel <- 
ggplot(ar_summary_all, 
       aes(x = AR_total_pct, y = brr_med, 
           color = age_group, fill = age_group)) +
  
  geom_ribbon(aes(ymin = brr_lo, ymax = brr_hi), 
              alpha = 0.2, color = NA) +
  
  geom_line(linewidth = 1.2) +
  
  geom_hline(yintercept = 1, linetype = "dashed", 
             color = "gray50", linewidth = 0.7) +
  
  # Fixed Y-axis within rows (outcome), but allow different scales between outcomes
  facet_grid(rows = vars(outcome), 
             cols = vars(days),
             scales = "free_y") +  # This keeps Y consistent within each row
  
  #scale_color_manual(
  #  values = c("<65" = "#E69F00", "65+" = "#56B4E9"),
  #  name = "Age group"
  #) +
  #scale_fill_manual(
  #  values = c("<65" = "#E69F00", "65+" = "#56B4E9"),
  #  name = "Age group"
  #) +
  
  scale_color_manual(
    values = c(
      "1-11"  = "#E69F00", 
      "12-17" = "#009E73", 
      "18-64" = "#56B4E9", 
      "65+"   = "#CC79A7"
    ),
    name = "Age group"
  ) +
  scale_fill_manual(
    values = c(
      "1-11"  = "#E69F00", 
      "12-17" = "#009E73", 
      "18-64" = "#56B4E9", 
      "65+"   = "#CC79A7"
    ),
    name = "Age group"
  ) +
  labs(
    x = "Total outbreak attack rate (%)",
    y = "Benefit-Risk Ratio",
    title = "Benefit-risk ratio: Traveller vaccination"
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    plot.caption = element_text(hjust = 0, size = 11, color = "gray40"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 10, hjust = 0)
  ) + scale_y_continuous(
    trans = pseudo_log_trans(base = 10),
    breaks = c(0, 1, 10, 100, 1000),
    labels = scales::comma_format()
  )

plot_brr_travel

### Benefit Harm trade off -----------------------------------------------------
# Risk space -------------------------------------------------------------------
br_space_data_correct <- psa_df %>%
  mutate(
    days = factor(days, levels = c("7d", "14d", "30d", "90d"))
  ) %>%
  pivot_longer(
    cols = c(risk_v_10k_sae, risk_v_10k_death, risk_v_10k_daly,
             risk_nv_10k_sae, risk_nv_10k_death, risk_nv_10k_daly,
             net_10k_sae, net_10k_death),
    names_to = c(".value", "outcome"),
    names_pattern = "(risk_v_10k|risk_nv_10k|net_10k)_(sae|death|daly)"
  ) %>%
  mutate(
    outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death", "daly" = "DALY")
  ) %>%
  group_by(AR_total_pct, age_group, days, outcome) %>%
  summarise(
    x_med = median(risk_v_10k),
    x_lo = quantile(risk_v_10k, 0.025),
    x_hi = quantile(risk_v_10k, 0.975),
    y_med = median(risk_nv_10k),
    y_lo = quantile(risk_nv_10k, 0.025),
    y_hi = quantile(risk_nv_10k, 0.975),
    net_med = median(net_10k),
    .groups = "drop"
  )


# no daly
br_space_data_correct <- psa_df %>%
  mutate(
    days = factor(days, levels = c("7d", "14d", "30d", "90d"))
  ) %>%
  pivot_longer(
    cols = c(risk_v_10k_sae, risk_v_10k_death,
             risk_nv_10k_sae, risk_nv_10k_death, 
             net_10k_sae, net_10k_death),
    names_to = c(".value", "outcome"),
    names_pattern = "(risk_v_10k|risk_nv_10k|net_10k)_(sae|death)"
  ) %>%
  mutate(
    outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death")
  ) %>%
  group_by(AR_total_pct, age_group, days, outcome) %>%
  summarise(
    x_med = median(risk_v_10k),
    x_lo = quantile(risk_v_10k, 0.025),
    x_hi = quantile(risk_v_10k, 0.975),
    y_med = median(risk_nv_10k),
    y_lo = quantile(risk_nv_10k, 0.025),
    y_hi = quantile(risk_nv_10k, 0.975),
    net_med = median(net_10k),
    .groups = "drop"
  )

representative_ar <- c(10, 30, 50, 80)

br_representative_correct <- br_space_data_correct %>%
  filter(AR_total_pct %in% representative_ar) %>%
  mutate(
    ar_category = factor(
      paste0("AR = ", AR_total_pct, "%"),
      levels = paste0("AR = ", representative_ar, "%")
    )
  ) %>%
  mutate(age_group = ifelse(age_group == "65", "65+", age_group))

br_representative_correct <- br_representative_correct %>%
  dplyr::rename(AgeCat = age_group)

# SIMPLE: Create ONE large background for all facets
# Get global ranges
x_max_global <- max(br_representative_correct$x_hi, na.rm = TRUE) * 1.1
y_max_global <- max(br_representative_correct$y_hi, na.rm = TRUE) * 1.1

cat("Global ranges:\n")
cat("  X: [0,", round(x_max_global, 2), "]\n")
cat("  Y: [0,", round(y_max_global, 2), "]\n")

# Create single background grid
bg_grid_unified <- expand.grid(
  x = seq(0, x_max_global, length.out = 300),
  y = seq(0, y_max_global, length.out = 300)
) %>%
  mutate(outcomes_averted = y - x)

max_averted <- max(abs(bg_grid_unified$outcomes_averted), na.rm = TRUE)

cat("\nBackground grid created:\n")
cat("  Points:", nrow(bg_grid_unified), "\n")
cat("  Outcomes averted range: [", 
    round(min(bg_grid_unified$outcomes_averted), 2), ",",
    round(max(bg_grid_unified$outcomes_averted), 2), "]\n")

# Diagnostic: Check error bar sizes for triangles (65+)
cat("\n=== Error Bar Diagnostic (65+ age group) ===\n")
triangles_data <- br_representative_correct %>% 
  filter(AgeCat == "65+") %>%
  mutate(
    x_range = x_hi - x_lo,
    y_range = y_hi - y_lo,
    x_range_pct = (x_range / x_med) * 100,
    y_range_pct = (y_range / y_med) * 100
  )

cat("X-axis error bar ranges (65+):\n")
print(triangles_data %>% 
  group_by(outcome) %>%
  summarise(
    min_x_range = min(x_range, na.rm = TRUE),
    max_x_range = max(x_range, na.rm = TRUE),
    mean_x_range = mean(x_range, na.rm = TRUE),
    mean_x_range_pct = mean(x_range_pct, na.rm = TRUE),
    .groups = "drop"
  ))

cat("\nY-axis error bar ranges (65+):\n")
print(triangles_data %>% 
  group_by(outcome) %>%
  summarise(
    min_y_range = min(y_range, na.rm = TRUE),
    max_y_range = max(y_range, na.rm = TRUE),
    mean_y_range = mean(y_range, na.rm = TRUE),
    mean_y_range_pct = mean(y_range_pct, na.rm = TRUE),
    .groups = "drop"
  ))

cat("\nY-axis scale vs error bar size:\n")
cat("  Y-axis max:", round(y_max_global, 2), "\n")
cat("  Mean Y error bar range (65+):", 
    round(mean(triangles_data$y_range, na.rm = TRUE), 2), "\n")
cat("  Error bar as % of Y-axis:", 
    round(mean(triangles_data$y_range, na.rm = TRUE) / y_max_global * 100, 3), "%\n")

# Plot
plot_risk_travel <- 
  ggplot() +
  # Unified background (same for all facets)
  geom_raster(
    data = bg_grid_unified,
    aes(x = x, y = y, fill = outcomes_averted),
    alpha = 0.85,
    interpolate = TRUE
  ) +
  
  scale_fill_gradientn(
    name = "Outcomes averted\nper 10,000 (median)",
    colours = c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0"),
    values = scales::rescale(c(-max_averted, -max_averted/4, 0, 
                               max_averted/4, max_averted)),
    limits = c(-max_averted, max_averted),
    na.value = "gray90"
  ) +
  
  # Break-even line
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 0.7, color = "black") +
  
  # Error bars - Y-axis (vertical)
  # Make more visible with thicker lines and relative width
  geom_errorbar(
    data = br_representative_correct,
    aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome),
    width = x_max_global * 0.01,  # 1% of x-axis range (instead of 0)
    linewidth = 0.8,  # Slightly thicker
    alpha = 0.8  # More opaque
  ) +
  # Error bars - X-axis (horizontal)
  geom_errorbarh(
    data = br_representative_correct,
    aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome),
    height = y_max_global * 0.01,  # 1% of y-axis range (instead of 0)
    linewidth = 0.8,  # Slightly thicker
    alpha = 0.8  # More opaque
  ) +
  
  # Points
  geom_point(
    data = br_representative_correct,
    aes(x = x_med, y = y_med, shape = AgeCat, color = outcome),
    size = 2, stroke = 0.7, fill = "white"
  ) +
  
  # CRITICAL: Free X only, Y is fixed!
  facet_grid(
    rows = vars(ar_category),
    cols = vars(days),
    scales = "free_x"  # Only X is free, Y is shared
  ) +
  
  # Scales
  #scale_color_manual(
  #  values = c("Death" = "#F8766D", "SAE" = "#00BFC4", "DALY"  = "#7B3CFF" ),
  #  name = "Chikungunya\nadverse events"
  #) +
  
  scale_color_manual(
    values = c(
      "SAE"   = "#1B7F1B",   
      "Death" = "#B8860B",   
      "DALY"  = "#7E2F3B"    
    ),
    name = "Outcome type"
  ) + 
  scale_shape_manual(
    values = c("18-64" = 22, "65+" = 24),
    name = "Age group"
  ) +
  
  guides(
    fill = guide_colorbar(order = 1, title.position = "top"),
    color = guide_legend(order = 2, title.position = "top"),
    shape = guide_legend(order = 3, title.position = "top")
  ) + 
  
  # Set Y limits explicitly (optional, for cleaner look)
  coord_cartesian(ylim = c(0, y_max_global)) +
  
  # Labels
  labs(
    x = "Risk of vaccine associated adverse events (per 10,000 individuals)",
    y = "Risk of severe outcomes following infection (per 10,000 individuals)",
  ) +
  
  # Theme
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "right",
    legend.box = "vertical",
    panel.spacing = unit(0.8, "lines")
  ) 

plot_risk_travel <- plot_risk_travel + 
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = '',
    tag_suffix = '.',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18, hjust = 0, vjust = 1)
    )
  )

ggsave("06_Results/plot_1A.pdf", plot = plot_risk_travel, width = 10, height = 7)

## no daly
br_no_daly <- 
  br_representative_correct %>%
  filter(outcome != "DALY")

plot_risk_travel <- 
  ggplot() +
  # Unified background (same for all facets)
  geom_raster(
    data = bg_grid_unified,
    aes(x = x, y = y, fill = outcomes_averted),
    alpha = 0.85,
    interpolate = TRUE
  ) +
  
  scale_fill_gradientn(
    name = "Outcomes averted\nper 10,000 (median)",
    colours = c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0"),
    values = scales::rescale(c(-max_averted, -max_averted/4, 0, 
                               max_averted/4, max_averted)),
    limits = c(-max_averted, max_averted),
    na.value = "gray90"
  ) +
  
  # Break-even line
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 0.7, color = "black") +
  
  # Error bars - Y-axis (vertical)
  # Make more visible with thicker lines and relative width
  geom_errorbar(
    data = br_representative_correct,
    aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome),
    width = x_max_global * 0.01,  # 1% of x-axis range (instead of 0)
    linewidth = 0.8,  # Slightly thicker
    alpha = 0.8  # More opaque
  ) +
  # Error bars - X-axis (horizontal)
  geom_errorbarh(
    data = br_representative_correct,
    aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome),
    height = y_max_global * 0.01,  # 1% of y-axis range (instead of 0)
    linewidth = 0.8,  # Slightly thicker
    alpha = 0.8  # More opaque
  ) +
  
  # Points
  geom_point(
    data = br_representative_correct,
    aes(x = x_med, y = y_med, shape = AgeCat, color = outcome),
    size = 2, stroke = 0.7, fill = "white"
  ) +
  
  # CRITICAL: Free X only, Y is fixed!
  facet_grid(
    rows = vars(ar_category),
    cols = vars(days),
    scales = "free_x"  # Only X is free, Y is shared
  ) +
  
  # Scales
  #scale_color_manual(
  #  values = c("Death" = "#F8766D", "SAE" = "#00BFC4", "DALY"  = "#7B3CFF" ),
  #  name = "Chikungunya\nadverse events"
  #) +
  
  scale_color_manual(
    values = c(
      "SAE"   = "#1B7F1B",   
      "Death" = "#B8860B",   
      "DALY"  = "#7E2F3B"    
    ),
    name = "Outcome type"
  ) + 
  scale_shape_manual(
    values = c("18-64" = 22, "65+" = 24),
    name = "Age group"
  ) +
  
  guides(
    fill = guide_colorbar(order = 1, title.position = "top"),
    color = guide_legend(order = 2, title.position = "top"),
    shape = guide_legend(order = 3, title.position = "top")
  ) + 
  
  # Set Y limits explicitly (optional, for cleaner look)
  coord_cartesian(ylim = c(0, y_max_global)) +
  
  # Labels
  labs(
    x = "Risk of vaccine associated adverse events (per 10,000 individuals)",
    y = "Risk of severe outcomes following infection (per 10,000 individuals)",
  ) +
  
  # Theme
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "right",
    legend.box = "vertical",
    panel.spacing = unit(0.8, "lines")
  ) 



### final benefit-risk space----------------------------------------------------
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

representative_ar <- c(10, 30, 50, 80)

br_representative_benefit <- br_space_data_benefit %>%
  filter(AR_total_pct %in% representative_ar) %>%
  mutate(
    ar_category = factor(
      paste0("AR = ", AR_total_pct, "%"),
      levels = paste0("AR = ", representative_ar, "%")
    )
  ) %>%
  mutate(age_group = ifelse(age_group == "65", "65+", age_group))

br_representative_benefit <- br_representative_benefit |>
  filter(!age_group %in% c("18-64", "65+")) |>
  droplevels()

br_representative_benefit <- br_representative_benefit %>%
  dplyr::rename(AgeCat = age_group) 
  

# 2. 각 패널별 적절한 축 범위 계산
panel_ranges <- br_representative_benefit %>%
  group_by(ar_category, days) %>%
  summarise(
    x_min = min(c(0, x_lo), na.rm = TRUE),  # 0 또는 실제 최소값 중 작은 값
    x_max = max(x_hi, na.rm = TRUE) * 1.05,  # 5% 여유
    y_min = min(c(0, y_lo), na.rm = TRUE),  # 0 또는 실제 최소값 중 작은 값
    y_max = max(y_hi, na.rm = TRUE) * 1.05,  # 5% 여유
    .groups = "drop"
  )

# 3. 각 패널별 배경 그리드 생성
bg_grid_optimized <- panel_ranges %>%
  rowwise() %>%
  do({
    panel_data <- .
    
    # 해당 패널에 맞는 그리드 생성
    x_seq <- seq(panel_data$x_min, panel_data$x_max, length.out = 200)
    y_seq <- seq(panel_data$y_min, panel_data$y_max, length.out = 200)
    
    grid <- expand.grid(x = x_seq, y = y_seq)
    
    # BRR 계산
    grid$brr <- with(grid, ifelse(x > 0, y / x, NA_real_))
    grid$log10_brr <- log10(grid$brr)
    
    # 패널 정보 추가
    grid$ar_category <- panel_data$ar_category
    grid$days <- panel_data$days
    
    grid
  }) %>%
  ungroup()

# 4. log10 BRR 범위 설정
log_min <- -2
log_max <- 2
log_range <- seq(log_min, log_max, by = 1)
brr_labels <- c("0.01", "0.1", "1", "10", "100")

# 5. 시각화
plot_final <- ggplot() +
  # 패널별 최적화된 배경 레이어
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
  
  # 에러바 및 포인트
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
  
  # y=x 라인 (데이터 레이어 다음에 추가하여 clipping 적용)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "black", linewidth = 0.9, alpha = 0.8) +
  
  # facet_wrap 사용 (독립적인 축, 4열)
  facet_wrap(
    ar_category ~ days, 
    scales = "free", 
    ncol = 4
  ) +
  
  scale_color_manual(
    values = c("Hospitalisation" = "#1B7F1B", "Death" = "#B8860B"), 
    name = "Outcome type"
  ) +
  scale_shape_manual(
    values = c("1-11" = 21, "12-17" = 22, "18-59" = 23, "60" = 24), 
    name = "Age group"
  ) +
  
  # x축과 y축: 여백 없이 꽉 채우기
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  
  # clip = "on"으로 변경하여 패널 경계 밖으로 나가는 요소 잘라내기
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

print(plot_final)

combined_plot <- plot_final / plot_brr_travel + 
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = '',
    tag_suffix = '.',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18, hjust = 0, vjust = 1)
    )
  )

combined_plot

ggsave("06_Results/brr_final_plot1.pdf", plot = combined_plot, width = 14, height = 18)



# BRR tables -------------------------------------------------------------------
# Prepare summary with representative ARs only
representative_ar <- c(10, 30, 50, 80)

ar_summary_representative <- ar_summary_all %>%
  filter(AR_total_pct %in% representative_ar) %>%
  mutate(
    # Format BRR with glue
    brr_formatted = glue("{format(round(brr_med, 2), nsmall = 2)} ({format(round(brr_lo, 2), nsmall = 2)}, {format(round(brr_hi, 2), nsmall = 2)})"),
    # Order outcome: Death first, then SAE
    outcome = factor(outcome, levels = c("Death", "SAE", "DALY")),
    # Format AR
    ar_label = paste0(AR_total_pct, "%")
  ) %>%
  arrange(outcome, AR_total_pct, age_group, days)

# Create wide format table
summary_table_wide <- ar_summary_representative %>%
  select(outcome, ar_label, age_group, days, brr_formatted) %>%
  pivot_wider(
    names_from = days,
    values_from = brr_formatted
  ) %>%
  arrange(outcome, ar_label, age_group)

# Display table
kable(summary_table_wide, 
      format = "html",
      col.names = c("Outcome", "Attack Rate", "Age Group", 
                    "7 days", "14 days", "30 days", "90 days"),
      caption = "Benefit-Risk Ratio by Outbreak Intensity and Travel Duration",
      align = c("l", "c", "c", "r", "r", "r", "r")) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    font_size = 12
  ) %>%
  pack_rows(index = table(summary_table_wide$outcome)) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2, bold = TRUE)


# no daly
ar_summary_representative <- ar_summary_all %>%
  filter(AR_total_pct %in% representative_ar) %>%
  mutate(
    # Format BRR with glue
    brr_formatted = glue("{format(round(brr_med, 1), nsmall = 1)} ({format(round(brr_lo, 1), nsmall = 1, trim = TRUE)}, {format(round(brr_hi, 1), nsmall = 1, trim = TRUE)})"),
    # Order outcome: Death first, then SAE
    outcome = factor(outcome, levels = c("Death", "SAE")),
    # Format AR
    ar_label = paste0(AR_total_pct, "%")
  ) %>%
  arrange(outcome, AR_total_pct, age_group, days)

# Create wide format table
summary_table_wide <- ar_summary_representative %>%
  select(outcome, ar_label, age_group, days, brr_formatted) %>%
  pivot_wider(
    names_from = days,
    values_from = brr_formatted
  ) %>%
  arrange(outcome, ar_label, age_group)

# Display table
kable(summary_table_wide, 
      format = "html",
      col.names = c("Outcome", "Attack Rate", "Age Group", 
                    "7 days", "14 days", "30 days", "90 days"),
      caption = "Benefit-Risk Ratio by Outbreak Intensity and Travel Duration",
      align = c("l", "c", "c", "r", "r", "r", "r")) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    font_size = 12
  ) %>%
  pack_rows(index = table(summary_table_wide$outcome)) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2, bold = TRUE)

# Save as file
write.csv(summary_table_wide, "06_Results/brr_summary_table.csv", row.names = FALSE)


################################################################################
# plot 3
################################################################################

lhs_draw <- lhs_sample %>%
  mutate(draw = row_number()) %>%
  select(draw, p_sae_vacc_u65, p_sae_vacc_65, p_death_vacc_u65, p_death_vacc_65)

plot_df_net <- psa_df %>%
  mutate(
    draw = as.integer(draw),
    age_group = recode(age_group, `18-64` = "<65", `65` = "65+"),
    days = factor(days, levels = c("7d","14d","30d","90d")),
    AR_pct_10 = round(AR_total_pct/10)*10
  ) %>%
  filter(AR_pct_10 %in% c(10,30,50,80)) %>%
  left_join(lhs_draw, by = "draw") %>%
  mutate(
    # age별 vaccine AE probability 선택
    p_sae_vacc   = if_else(age_group == "65+", p_sae_vacc_65,   p_sae_vacc_u65),
    p_death_vacc = if_else(age_group == "65+", p_death_vacc_65, p_death_vacc_u65),
    
    # averted disease outcomes
    averted_sae_10k   = risk_nv_10k_sae   - risk_v_10k_sae,
    averted_death_10k = risk_nv_10k_death - risk_v_10k_death,
    
    # vaccine adverse outcomes per 10k vaccinated
    vacc_sae_10k   = 1e4 * p_sae_vacc,
    vacc_death_10k = 1e4 * p_death_vacc,
    
    # net averted
    net_sae_10k   = averted_sae_10k   - vacc_sae_10k,
    net_death_10k = averted_death_10k - vacc_death_10k
  ) %>%
  transmute(draw, age_group, days, AR_pct_10,
            SAE = net_sae_10k,
            Death = net_death_10k) %>%
  pivot_longer(cols = c(SAE, Death), names_to = "outcome", values_to = "net_10k") %>%
  group_by(age_group, days, AR_pct_10, outcome) %>%
  summarise(
    med = median(net_10k, na.rm=TRUE),
    lo  = quantile(net_10k, 0.025, na.rm=TRUE),
    hi  = quantile(net_10k, 0.975, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    direction = if_else(med >= 0, "Benefit (net averted)", "Harm (net excess)"),
    AR_label = factor(paste0(AR_pct_10, "%"), levels = paste0(c(10,30,50,80), "%"))
  )

plot_net_travel <- ggplot(plot_df_net, aes(x = med, y = AR_label, fill = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_col() +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, alpha = 0.6) +
  facet_grid(outcome ~ age_group + days, scales = "free_x") +
  labs(
    #x = "Net averted outcomes per 10,000 vaccinated",
    y = "Total outbreak attack rate",
    fill = NULL
  ) +
  xlab("")+
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

write.csv(plot_df_net,
          file = "06_Results/net_benefit_travel.csv",
          row.names = FALSE)
