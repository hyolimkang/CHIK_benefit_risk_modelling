 
# ==============================================================================
# STEP 1: Risk summary (vaccine risks only)
# ==============================================================================
risk_summary <- lhs_sample %>%
  summarise(
    sae_u65_med = median(p_sae_vacc_u65),
    sae_u65_lo  = quantile(p_sae_vacc_u65, 0.025),
    sae_u65_hi  = quantile(p_sae_vacc_u65, 0.975),
    
    sae_65_med  = median(p_sae_vacc_65),
    sae_65_lo   = quantile(p_sae_vacc_65, 0.025),
    sae_65_hi   = quantile(p_sae_vacc_65, 0.975),
    
    death_u65_med = median(p_death_vacc_u65),
    death_u65_lo  = quantile(p_death_vacc_u65, 0.025),
    death_u65_hi  = quantile(p_death_vacc_u65, 0.975),
    
    death_65_med  = median(p_death_vacc_65),
    death_65_lo   = quantile(p_death_vacc_65, 0.025),
    death_65_hi   = quantile(p_death_vacc_65, 0.975),
    
    # daly related params -------------------------
    le_lost_1_11  = median(le_lost_1_11),
    le_lost_12_17 = median(le_lost_12_17),
    le_lost_18_64 = median(le_lost_18_64),
    le_lost_65    = median(le_lost_65),
    
    dw_hosp    = median(dw_hosp),
    dw_nonhosp = median(dw_nonhosp),
    dw_subac   = median(dw_subac),
    dw_chronic = median(dw_chronic),
    
    dur_acute   = median(dur_acute),
    dur_nonhosp = median(dur_nonhosp),
    dur_subac   = median(dur_subac),
    dur_6m      = median(dur_6m),
    dur_12m     = median(dur_12m),
    dur_30m     = median(dur_30m),
    
    acute  = median(acute),
    subac  = median(subac),
    chr6m  = median(chr6m),
    chr12m = median(chr12m),
    chr30m = median(chr30m)
  )
# ==============================================================================
# STEP 2: Prepare df2 
# ==============================================================================
# 1. Aggregate counts by Age Category and Scenario first
combined_nnv_national_vc50 <- combined_nnv_national_age_ixchiq %>% filter(VC == "cov50")

df_aggregated <- combined_nnv_national_vc50 %>%
  mutate(
    # Define analysis categories
    AgeCat = case_when(
      AgeGroup %in% 2:4 ~ "1-11",
      AgeGroup == 5 ~ "12-17",
      AgeGroup %in% 6:15 ~ "18-64",
      AgeGroup %in% 16:20 ~ "65+",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(AgeCat)) %>%
  group_by(scenario, AgeCat, VE) %>%
  summarise(
    # Sum of Averted outcomes (Benefit side)
    averted_death_med = sum(diff_fatal, na.rm = TRUE),
    averted_death_lo  = sum(diff_fatal_low, na.rm = TRUE),
    averted_death_hi  = sum(diff_fatal_hi, na.rm = TRUE),
    
    averted_hosp_med  = sum(diff_hosp, na.rm = TRUE),
    averted_hosp_lo   = sum(diff_hosp_lo, na.rm = TRUE),
    averted_hosp_hi   = sum(diff_hosp_hi, na.rm = TRUE),
    
    averted_daly_med  = sum(diff_daly, na.rm = TRUE),
    averted_daly_lo   = sum(diff_daly_low, na.rm = TRUE),
    averted_daly_hi   = sum(diff_daly_hi, na.rm = TRUE),
    
    # Total vaccinated in this combined group
    tot_vacc_grp = sum(tot_vacc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 2. Map Risk Probabilities and calculate Absolute Risk Events
  mutate(
    # Use risk_summary values based on AgeCat
    risk_sae_med = ifelse(AgeCat == "65+", risk_summary$sae_65_med, risk_summary$sae_u65_med),
    risk_sae_lo  = ifelse(AgeCat == "65+", risk_summary$sae_65_lo,  risk_summary$sae_u65_lo),
    risk_sae_hi  = ifelse(AgeCat == "65+", risk_summary$sae_65_hi,  risk_summary$sae_u65_hi),
    
    risk_death_med = ifelse(AgeCat == "65+", risk_summary$death_65_med, risk_summary$death_u65_med),
    risk_death_lo  = ifelse(AgeCat == "65+", risk_summary$death_65_lo,  risk_summary$death_u65_lo),
    risk_death_hi  = ifelse(AgeCat == "65+", risk_summary$death_65_hi,  risk_summary$death_u65_hi),
    
    # Calculate Total Absolute Events for the whole AgeCat
    vax_sae_med   = tot_vacc_grp * risk_sae_med,
    vax_sae_lo    = tot_vacc_grp * risk_sae_lo,
    vax_sae_hi    = tot_vacc_grp * risk_sae_hi,
    
    vax_death_med = tot_vacc_grp * risk_death_med,
    vax_death_lo  = tot_vacc_grp * risk_death_lo,
    vax_death_hi  = tot_vacc_grp * risk_death_hi
  ) %>%
  mutate(
    # Assign target group based on scenario
    target = case_when(
      scenario == "Scenario_1" & AgeCat == "1-11"  ~ 1,
      scenario == "Scenario_2" & AgeCat == "12-17" ~ 1,
      scenario == "Scenario_3" & AgeCat == "18-64" ~ 1,
      scenario == "Scenario_4" & AgeCat == "65+"   ~ 1,
      TRUE ~ 0 # Non-vaccinated groups (beneficiaries of indirect effects)
    )
  )


df_final_averted <- df_aggregated %>%
  group_by(scenario, VE) %>%
  mutate(
    across(starts_with("averted_"), 
           ~ if(first(VE) > 0) sum(.x, na.rm = TRUE) else .x)
  ) %>%
  ungroup() %>%
  filter(target == 1)

df_final <- df_final_averted %>%
  rowwise() %>%
  mutate(
    # --- Medium Scenario Risk DALY ---
    res_med = list(compute_daly_one(
      age_group      = AgeCat,
      sae_10k        = vax_sae_med,   # Events in target group
      deaths_sae_10k = vax_death_med, # Deaths in target group
      draw_pars      = as.list(risk_summary)
    )),
    vax_daly_med = res_med$daly_sae,
    
    # --- Low Scenario Risk DALY ---
    res_lo = list(compute_daly_one(
      age_group      = AgeCat,
      sae_10k        = vax_sae_lo,
      deaths_sae_10k = vax_death_lo,
      draw_pars      = as.list(risk_summary)
    )),
    vax_daly_lo = res_lo$daly_sae,
    
    # --- High Scenario Risk DALY ---
    res_hi = list(compute_daly_one(
      age_group      = AgeCat,
      sae_10k        = vax_sae_hi,
      deaths_sae_10k = vax_death_hi,
      draw_pars      = as.list(risk_summary)
    )),
    vax_daly_hi = res_hi$daly_sae
  ) %>%
  ungroup() %>%
  select(-starts_with("res_"))

summary_df <- df_final %>%
  transmute(
    scenario, AgeCat, VE, tot_vacc_grp,
    
    # --- DALY ---
    x_DALY_med = (vax_daly_med / tot_vacc_grp) * 1e4,
    x_DALY_lo  = (vax_daly_lo  / tot_vacc_grp) * 1e4,
    x_DALY_hi  = (vax_daly_hi  / tot_vacc_grp) * 1e4,
    y_DALY_med = (averted_daly_med / tot_vacc_grp) * 1e4,
    y_DALY_lo  = (averted_daly_lo  / tot_vacc_grp) * 1e4,
    y_DALY_hi  = (averted_daly_hi  / tot_vacc_grp) * 1e4,
    
    # --- Death ---
    x_Death_med = (vax_death_med / tot_vacc_grp) * 1e4,
    x_Death_lo  = (vax_death_lo  / tot_vacc_grp) * 1e4,
    x_Death_hi  = (vax_death_hi  / tot_vacc_grp) * 1e4,
    y_Death_med = (averted_death_med / tot_vacc_grp) * 1e4,
    y_Death_lo  = (averted_death_lo  / tot_vacc_grp) * 1e4,
    y_Death_hi  = (averted_death_hi  / tot_vacc_grp) * 1e4,
    
    # --- SAE  ---
    x_SAE_med = (vax_sae_med / tot_vacc_grp) * 1e4,
    x_SAE_lo  = (vax_sae_lo  / tot_vacc_grp) * 1e4,
    x_SAE_hi  = (vax_sae_hi  / tot_vacc_grp) * 1e4,
    y_SAE_med = (averted_hosp_med / tot_vacc_grp) * 1e4,
    y_SAE_lo  = (averted_hosp_lo  / tot_vacc_grp) * 1e4,
    y_SAE_hi  = (averted_hosp_hi  / tot_vacc_grp) * 1e4
  ) %>%
  pivot_longer(
    cols = starts_with(c("x_", "y_")),
    names_to = c(".value", "outcome", "stat"),
    names_sep = "_"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = c(x, y)
  )


# plot--------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 1. Data Preparation and Labeling
# ------------------------------------------------------------------------------
log_min <- -2
log_max <- 2
log_range <- seq(log_min, log_max, by = 1)
brr_labels <- c("0.01", "0.1", "1", "10", "100")

# Convert VE to factors with the requested descriptive labels
summary_df <- summary_df %>%
  mutate(VE_label = factor(VE, 
                           levels = c("VE0", "VE98.9"),
                           labels = c("Disease blocking only", "Disease and infection blocking")))

# Define scenario labels for facet headers
scenario_labels <- c(
  "Scenario_1" = "1–11 years",
  "Scenario_2" = "12–17 years",
  "Scenario_3" = "18–64 years",
  "Scenario_4" = "65+ years"
)

# ------------------------------------------------------------------------------
# 2. Global Limit Calculation with Buffer (To Prevent Clipping)
# ------------------------------------------------------------------------------

# Calculate global maximums and add a 15% buffer to prevent shapes from being cut off
x_max_buffered <- max(summary_df$x_hi, na.rm = TRUE) * 1.15
y_max_buffered <- max(summary_df$y_hi, na.rm = TRUE) * 1.15

# Set a small negative buffer for the origin to ensure points at 0 are fully visible
x_min_buffered <- -x_max_buffered * 0.01
y_min_buffered <- -y_max_buffered * 0.01


# ------------------------------------------------------------------------------
# 4. Final Visualization
# ------------------------------------------------------------------------------

create_br_plot <- function(data, target_outcome, log_min = -1, log_max = 3) {
  
  # 1. outcome filtering
  plot_data <- data %>% filter(outcome == target_outcome)
  
  # 2. outcome specific limits
  panel_limits <- plot_data %>%
    group_by(AgeCat) %>%
    summarise(
      x_max = max(x_hi, na.rm = TRUE) * 1.1,
      y_max = max(y_hi, na.rm = TRUE) * 1.1,
      .groups = "drop"
    )
  
  bg_grid_specific <- panel_limits %>%
    group_by(AgeCat) %>%
    reframe({
      x_seq <- seq(0, x_max, length.out = 100)
      y_seq <- seq(0, y_max, length.out = 100)
      grid <- expand.grid(x = x_seq, y = y_seq)
      grid$brr <- grid$y / (grid$x + 1e-9)
      grid$log10_brr <- pmax(pmin(log10(grid$brr), log_max), log_min)
      grid
    })
  
  # 4. ggplot 
  ggplot() +
    geom_raster(
      data = bg_grid_specific,
      aes(x = x, y = y, fill = log10_brr),
      alpha = 0.8, interpolate = TRUE
    ) +
    scale_fill_gradient2(
      name = "Benefit–risk ratio",
      low = "#ca0020", mid = "#f7f7f7", high = "#0571b0", midpoint = 0,
      limits = c(log_min, log_max),
      oob = scales::squish,
      breaks = log_range, 
      labels = brr_labels,
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.4) +
    
    # error bar and points
    geom_errorbar(data = plot_data, aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome), width = 0) +
    geom_errorbarh(data = plot_data, aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome), height = 0) +
    geom_point(data = plot_data, aes(x = x_med, y = y_med, shape = VE_label, color = outcome), 
               size = 3, fill = "white", stroke = 1) +
    
    # facet by age group
    facet_wrap(~ AgeCat, scales = "free", ncol = 2) +
    
    scale_color_manual(values = c("SAE" = "#1B7F1B", "Death" = "#B8860B", "DALY" = "#A23B72")) +
    scale_shape_manual(values = c("Disease blocking only" = 21, "Disease and infection blocking" = 24)) +
    
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    
    labs(
      title = paste("Benefit-risk Assessment:", target_outcome),
      x = "Excess adverse outcomes (per 10,000 vaccinated)",
      y = "Averted outcomes (per 10,000 vaccinated)"
    ) +
    theme_bw() +
    theme(legend.position = "right", strip.text = element_text(face = "bold"))
}

plot_sae   <- create_br_plot(summary_df, "SAE")
plot_death <- create_br_plot(summary_df, "Death")
plot_daly  <- create_br_plot(summary_df, "DALY")

plot_daly
plot_sae
plot_death

combined_plot <- plot_vaccine_mechanism_final / plot_brr_ori + 
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = '',
    tag_suffix = '.',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18, hjust = 0, vjust = 1)
    )
  )

ggsave("06_Results/brr_final_plot2.pdf", plot = combined_plot, width = 10, height = 11)
