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
  # 1. 각 Outcome별로 x(Risk)와 y(Benefit)를 1만 명당 수치로 계산
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
    
    # --- SAE (Benefit은 averted_hosp 사용) ---
    x_SAE_med = (vax_sae_med / tot_vacc_grp) * 1e4,
    x_SAE_lo  = (vax_sae_lo  / tot_vacc_grp) * 1e4,
    x_SAE_hi  = (vax_sae_hi  / tot_vacc_grp) * 1e4,
    y_SAE_med = (averted_hosp_med / tot_vacc_grp) * 1e4,
    y_SAE_lo  = (averted_hosp_lo  / tot_vacc_grp) * 1e4,
    y_SAE_hi  = (averted_hosp_hi  / tot_vacc_grp) * 1e4
  ) %>%
  # 2. Outcome별로 데이터 재구조화 (Tidy format)
  pivot_longer(
    cols = starts_with(c("x_", "y_")),
    names_to = c(".value", "outcome", "stat"),
    names_sep = "_"
  ) %>%
  # 3. 에러바를 그리기 위해 다시 가로로 살짝 넓히기 (med, lo, hi를 컬럼으로)
  pivot_wider(
    names_from = stat,
    values_from = c(x, y)
  )
# ==============================================================================
# STEP 3: Create risk space data (CORRECTED!)
# ==============================================================================

death_df <- df_final %>%
  transmute(
    scenario, VE, AgeCat, outcome = "Death",
    tot_vacc,
    
    # X = Vaccine adverse deaths
    x_med = vax_death_med,
    x_lo  = vax_death_lo,
    x_hi  = vax_death_hi,
    
    # Y = Deaths averted
    y_med = averted_death_med,
    y_lo  = averted_death_lo,
    y_hi  = averted_death_hi
  )

sae_df <- df2_vc50 %>%
  transmute(
    scenario, VE, VC, AgeCat, outcome = "SAE",
    tot_vacc,
    
    # X = Vaccine adverse events
    x_med = vax_sae_med,
    x_lo  = vax_sae_lo,
    x_hi  = vax_sae_hi,
    
    # Y = Hospitalizations without vaccine (CORRECTED: use pre_hosp!)
    y_med = pre_hosp,
    y_lo  = pre_hosp_lo,
    y_hi  = pre_hosp_hi
  )

daly_df <- df2_vc50 %>%
  transmute(
    scenario, VE, VC, AgeCat, outcome = "DALY",
    tot_vacc,
    
    # X = Vaccine adverse events
    x_med = vax_daly_med,
    x_lo  = vax_daly_lo,
    x_hi  = vax_daly_hi,
    
    # Y = Hospitalizations without vaccine (CORRECTED: use pre_hosp!)
    y_med = pre_daly,
    y_lo  = pre_daly_lo,
    y_hi  = pre_daly_hi
  )

summary_df <- bind_rows(death_df, sae_df, daly_df)



# ==============================================================================
# STEP 5: Convert to per 10k
# ==============================================================================
summary_df2 <- summary_df %>%
  group_by(scenario, VE, VC, outcome, AgeCat) %>%
  summarise(
    x_med = sum(x_med, na.rm = TRUE),
    x_lo  = sum(x_lo,  na.rm = TRUE),
    x_hi  = sum(x_hi,  na.rm = TRUE),
    y_med = sum(y_med, na.rm = TRUE),
    y_lo  = sum(y_lo,  na.rm = TRUE),
    y_hi  = sum(y_hi,  na.rm = TRUE),
    tot_vacc_grp = sum(tot_vacc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Per 10k
    x_med = (x_med / tot_vacc_grp) * 1e4,
    x_lo  = (x_lo  / tot_vacc_grp) * 1e4,
    x_hi  = (x_hi  / tot_vacc_grp) * 1e4,
    
    y_med = (y_med / tot_vacc_grp) * 1e4,
    y_lo  = (y_lo  / tot_vacc_grp) * 1e4,
    y_hi  = (y_hi  / tot_vacc_grp) * 1e4
      ) %>%
  filter(is.finite(x_med) & is.finite(y_med))

# ==============================================================================
# STEP 6: Split by VE
# ==============================================================================
summary_by_ve <- summary_df2 %>%
  #filter(VC == "cov90") %>%
  split(.$VE) 

# ==============================================================================
# STEP 7: Create unified background
# ==============================================================================
x_vals <- c(summary_df$x_lo, summary_df$x_med, summary_df$x_hi)
y_vals <- c(summary_df$y_lo, summary_df$y_med, summary_df$y_hi)

x_max_global <- max(x_vals[is.finite(x_vals)], na.rm = TRUE) * 1.1
y_max_global <- max(y_vals[is.finite(y_vals)], na.rm = TRUE) * 1.1

bg_grid_unified <- expand.grid(
  x = seq(0, x_max_global, length.out = 300),
  y = seq(0, y_max_global, length.out = 300)
) %>%
  mutate(outcomes_averted = y - x)

max_averted <- max(abs(bg_grid_unified$outcomes_averted), na.rm = TRUE)
cat("Max averted:", round(max_averted, 2), "\n")

# ==============================================================================
# STEP 8: Create plots
# ==============================================================================
scenario_labels <- c(
  "Scenario_1" = "1–11 years",
  "Scenario_2" = "12–17 years",
  "Scenario_3" = "18–64 years",
  "Scenario_4" = "65+ years"
)


df <- summary_by_ve[["VE0"]]

df <- df %>%
  mutate(
    AgeCat = case_when(
      AgeCat == "<65" ~ "18-64",
      TRUE ~ AgeCat
    )
  ) 

plot_risk_ori <- ggplot() +
  # Background gradient
  geom_raster(
    data = bg_grid_unified,
    aes(x = x, y = y, fill = outcomes_averted),
    alpha = 0.85, interpolate = TRUE
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
  
  # Error bars
  geom_errorbar(
    data = df,
    aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome),
    width = 0, linewidth = 0.6, alpha = 0.7
  ) +
  geom_errorbarh(
    data = df,
    aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome),
    height = 0, linewidth = 0.6, alpha = 0.7
  ) +
  
  # Points
  geom_point(
    data = df,
    aes(x = x_med, y = y_med, shape = AgeCat, color = outcome),
    size = 2, stroke = 0.7, fill = "white"
  ) +
  
  # Facets
  facet_wrap(~scenario, labeller = as_labeller(scenario_labels)) +
  
  # Scales
  scale_color_manual(
    values = c("Death" = "#F8766D", "SAE" = "#00BFC4"),
    name = "Chikungunya\nadverse events"
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
  
  coord_cartesian(xlim = c(0, x_max_global), ylim = c(0, y_max_global)) +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_y_continuous(labels = scales::comma_format()) +
  
  labs(
    x = "Risk of vaccine associated adverse events (per 10,000 individuals)",
    y = "Risk of severe outcomes following infection (per 10,000 individuals)",
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "right",
    legend.box = "vertical",
    panel.spacing = unit(0.8, "lines")
  )



plot_risk_ori <- plot_risk_ori +
  labs(tag = "B.") +
  theme(
    plot.tag = element_text(size = 18),
    plot.tag.position = c(0, 1) 
  )

# BRR ==============================================================================

brr_summary <- df2_vc50 %>%
  group_by(scenario, VE, AgeCat) %>%
  summarise(
    tot_vacc_grp = sum(tot_vacc, na.rm = TRUE),
    
    excess_death_med_sum = sum(vax_death_med, na.rm = TRUE),
    excess_death_lo_sum  = sum(vax_death_lo,  na.rm = TRUE),
    excess_death_hi_sum  = sum(vax_death_hi,  na.rm = TRUE),
    
    averted_death_med_sum = sum(averted_death_med, na.rm = TRUE),
    averted_death_lo_sum  = sum(averted_death_lo,  na.rm = TRUE),
    averted_death_hi_sum  = sum(averted_death_hi,  na.rm = TRUE),
    
    excess_hosp_med_sum = sum(vax_sae_med, na.rm = TRUE),
    excess_hosp_lo_sum  = sum(vax_sae_lo,  na.rm = TRUE),
    excess_hosp_hi_sum  = sum(vax_sae_hi,  na.rm = TRUE),
    
    averted_hosp_med_sum = sum(averted_hosp_med, na.rm = TRUE),
    averted_hosp_lo_sum  = sum(averted_hosp_lo,  na.rm = TRUE),
    averted_hosp_hi_sum  = sum(averted_hosp_hi,  na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    # per 10k (Death)
    excess_death_10k_med = ifelse(tot_vacc_grp > 0, excess_death_med_sum / tot_vacc_grp * 1e4, NA_real_),
    excess_death_10k_lo  = ifelse(tot_vacc_grp > 0, excess_death_lo_sum  / tot_vacc_grp * 1e4, NA_real_),
    excess_death_10k_hi  = ifelse(tot_vacc_grp > 0, excess_death_hi_sum  / tot_vacc_grp * 1e4, NA_real_),
    
    averted_death_10k_med = ifelse(tot_vacc_grp > 0, averted_death_med_sum / tot_vacc_grp * 1e4, NA_real_),
    averted_death_10k_lo  = ifelse(tot_vacc_grp > 0, averted_death_lo_sum  / tot_vacc_grp * 1e4, NA_real_),
    averted_death_10k_hi  = ifelse(tot_vacc_grp > 0, averted_death_hi_sum  / tot_vacc_grp * 1e4, NA_real_),
    
    # per 10k (Hospitalisation)
    excess_hosp_10k_med = ifelse(tot_vacc_grp > 0, excess_hosp_med_sum / tot_vacc_grp * 1e4, NA_real_),
    excess_hosp_10k_lo  = ifelse(tot_vacc_grp > 0, excess_hosp_lo_sum  / tot_vacc_grp * 1e4, NA_real_),
    excess_hosp_10k_hi  = ifelse(tot_vacc_grp > 0, excess_hosp_hi_sum  / tot_vacc_grp * 1e4, NA_real_),
    
    averted_hosp_10k_med = ifelse(tot_vacc_grp > 0, averted_hosp_med_sum / tot_vacc_grp * 1e4, NA_real_),
    averted_hosp_10k_lo  = ifelse(tot_vacc_grp > 0, averted_hosp_lo_sum  / tot_vacc_grp * 1e4, NA_real_),
    averted_hosp_10k_hi  = ifelse(tot_vacc_grp > 0, averted_hosp_hi_sum  / tot_vacc_grp * 1e4, NA_real_),
    
    brr_death_med = ifelse(excess_death_med_sum > 0, averted_death_med_sum / excess_death_med_sum, NA_real_),
    brr_death_lo  = ifelse(excess_death_hi_sum  > 0, averted_death_lo_sum  / excess_death_hi_sum,  NA_real_),
    brr_death_hi  = ifelse(excess_death_lo_sum  > 0, averted_death_hi_sum  / excess_death_lo_sum,  NA_real_),
    
    brr_hosp_med  = ifelse(excess_hosp_med_sum  > 0, averted_hosp_med_sum  / excess_hosp_med_sum,  NA_real_),
    brr_hosp_lo   = ifelse(excess_hosp_hi_sum   > 0, averted_hosp_lo_sum   / excess_hosp_hi_sum,   NA_real_),
    brr_hosp_hi   = ifelse(excess_hosp_lo_sum   > 0, averted_hosp_hi_sum   / excess_hosp_lo_sum,   NA_real_)
  )

# ==============================================================================
# STEP 4: Pivot to long format for plotting
# ==============================================================================
brr_long <- brr_summary %>%
  pivot_longer(
    cols = c(starts_with("brr_death"), starts_with("brr_hosp"),
             starts_with("averted_death"), starts_with("averted_hosp"),
             starts_with("excess_death"), starts_with("excess_hosp")),
    names_to = c("metric", "outcome", ".value"),
    names_pattern = "(.+)_(death|hosp)_(.+)"
  ) %>%
  mutate(
    outcome = dplyr::recode(outcome, "death" = "Death", "hosp" = "Hospitalisation"),
    AgeCat = factor(AgeCat, levels = c("<65", "65+"))
  )

# ==============================================================================
# STEP 6: BRR Line Plot
# ==============================================================================
brr_plot_data <- brr_long %>%
  filter(metric == "brr") %>%
  filter(is.finite(med) & med > 0)%>%
  mutate(
    scenario = dplyr::recode(
      scenario,
      "Scenario_1" = "1–11",
      "Scenario_2" = "12–17",
      "Scenario_3" = "18–59",
      "Scenario_4" = "60+"
    ),
    outcome = dplyr::recode(
      outcome,
      "Hospitalisation"  = "SAE"
    )
  )

plot_brr_ori <- 
  ggplot(brr_plot_data, 
         aes(x = VE, y = med, color = VE, group = VE)) + # color와 group을 VE로 변경
  
  # Break-even line
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "gray30") +
  
  # Error bars
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.15, alpha = 0.7,
                position = position_dodge(width = 0.3)) +
  
  # Points
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  
  # Facets
  facet_grid(outcome ~ scenario, scales = "free_y") +
  
  # Scales (색상 기준을 VE 값으로 변경)
  scale_color_manual(
    values = c("VE0" = "#E69F00", "VE98.9" = "#56B4E9"), # 데이터의 실제 VE 값에 색상 매핑
    name = "Vaccine protection mechanism", # 레전드 제목 변경
    labels = c(
      "VE0"    = "Disease blocking only",
      "VE98.9" = "Disease and infection blocking"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "VE0"    = "Disease blocking only",
      "VE98.9" = "Disease and infection blocking"
    )
  ) +
  scale_y_log10(labels = scales::comma) +
  
  # Labels
  labs(
    x = "Vaccine protection mechanism",
    y = "Benefit–Risk Ratio",
    title = "Benefit-Risk ratio: Outbreak response immunisation"
  ) +
  
  # Theme
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(hjust = 0, size = 11, lineheight = 1.2, color = "gray40"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

print(plot_brr_ori)
ggsave("06_Results/plot2.pdf", plot = combined_plot, width = 10, height = 8)


##### 
brr_table_input <- brr_plot_data %>%
  mutate(
    outcome_disp = dplyr::recode(outcome,
                                 "Hospitalisation" = "SAE",
                                 "Death"           = "Death"),
    outcome_disp = factor(outcome_disp, levels = c("Death", "SAE")),
    AgeCat       = factor(AgeCat, levels = c("<65", "65+")),
    brr_formatted = sprintf("%.2f (%.2f, %.2f)", med, lo, hi)
  )

brr_summary_table_wide <- brr_table_input %>%
  select(outcome = outcome_disp, AgeCat, scenario, VE, brr_formatted) %>%
  arrange(outcome, scenario, AgeCat, VE) %>%
  pivot_wider(
    names_from  = VE,
    values_from = brr_formatted
  )

print(brr_summary_table_wide)

write.csv(brr_summary_table_wide, "06_Results/brr_summary_table_ori.csv", row.names = FALSE)


brr_tables_by_scenario <- brr_table_input %>%
  select(outcome = outcome_disp, AgeCat, scenario, VE, brr_formatted) %>%
  arrange(outcome, AgeCat, VE) %>%
  split(.$scenario) %>%
  lapply(function(d) {
    d %>%
      select(-scenario) %>%
      pivot_wider(
        names_from  = VE,
        values_from = brr_formatted
      )
  })


################################################################################

death_df <- df2_vc50 %>%
  transmute(
    scenario, VE, VC, AgeCat, outcome = "Death",
    tot_vacc,
    x_med = vax_death_med,
    x_lo  = vax_death_lo,
    x_hi  = vax_death_hi,
    y_med = diff_fatal,
    y_lo  = diff_fatal_low,
    y_hi  = diff_fatal_hi
  )

hosp_df <- df2_vc50 %>%
  transmute(
    scenario, VE, VC, AgeCat, outcome = "SAE",  
    tot_vacc,
    x_med = vax_sae_med,
    x_lo  = vax_sae_lo,
    x_hi  = vax_sae_hi,
    y_med = diff_hosp,
    y_lo  = diff_hosp_lo,
    y_hi  = diff_hosp_hi
  )

summary_df <- bind_rows(death_df, hosp_df)

summary_df2 <- summary_df %>%
  group_by(scenario, VE, VC, outcome, AgeCat) %>%
  summarise(
    x_med = sum(x_med, na.rm = TRUE),
    x_lo  = sum(x_lo,  na.rm = TRUE),
    x_hi  = sum(x_hi,  na.rm = TRUE),
    y_med = sum(y_med, na.rm = TRUE),
    y_lo  = sum(y_lo,  na.rm = TRUE),
    y_hi  = sum(y_hi,  na.rm = TRUE),
    tot_vacc_grp = sum(tot_vacc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_med = (x_med / tot_vacc_grp) * 1e4,
    x_lo  = (x_lo  / tot_vacc_grp) * 1e4,
    x_hi  = (x_hi  / tot_vacc_grp) * 1e4,
    y_med = (y_med / tot_vacc_grp) * 1e4,
    y_lo  = (y_lo  / tot_vacc_grp) * 1e4,
    y_hi  = (y_hi  / tot_vacc_grp) * 1e4,
    net_med = y_med - x_med
  ) %>%
  filter(is.finite(x_med) & is.finite(y_med))

net_df <- summary_df2 %>%
  mutate(
    outcome = recode(outcome, "Hospitalisation" = "SAE"),
    VE_label = recode(as.character(VE),
                      "VE0" = "Disease blocking only",
                      "VE98.9" = "Disease and infection blocking",
                      .default = as.character(VE)),
    
    scenario = recode(as.character(scenario),
                      "Scenario_1" = "1–11 years",
                      "Scenario_2" = "12–17 years",
                      "Scenario_3" = "18–59 years",
                      "Scenario_4" = "60+ years",
                      .default = as.character(scenario)),
    scenario = factor(scenario, levels = c("1–11 years","12–17 years","18–59 years","60+ years")),
    
    net_lo = y_lo - x_hi,
    net_hi = y_hi - x_lo,
    
    direction = if_else(net_med >= 0, "Benefit (net averted)", "Harm (net excess)")
  )

# 3) Net averted plot
plot_net_outbreak <- ggplot(
  net_df,
  aes(x = net_med, y = VE_label, fill = direction)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_col(width = 0.7) +
  geom_errorbarh(aes(xmin = net_lo, xmax = net_hi), height = 0.18, alpha = 0.6) +
  
  facet_grid(outcome ~ scenario, scales = "free_x") +
  
  scale_fill_manual(
    values = c("Benefit (net averted)" = "#F8766D",
               "Harm (net excess)"     = "#00BFC4"),
    name = NULL
  ) +
  scale_x_continuous(labels = scales::comma) +
  
  labs(
    x = "Net averted outcomes per 10,000 vaccinated",
    y = "Vaccine protection mechanism"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


plot_net_outbreak

combined_plot <- plot_net_travel / plot_net_outbreak + 
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = '',
    tag_suffix = '.',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18, hjust = 0, vjust = 1)
    )
  )

ggsave("06_Results/plot3.pdf", plot = combined_plot, width = 14, height = 7)

write.csv(net_df,
          file = "06_Results/net_benefit_ori.csv",
          row.names = FALSE)
