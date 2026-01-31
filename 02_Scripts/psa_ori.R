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
    death_65_hi   = quantile(p_death_vacc_65, 0.975)
    
  )

df2 <- combined_nnv_national_age_ixchiq %>%
  mutate(
    risk_sae_med   = ifelse(AgeGroup < 16, risk_summary$sae_u65_med, risk_summary$sae_65_med),
    risk_sae_lo    = ifelse(AgeGroup < 16, risk_summary$sae_u65_lo,  risk_summary$sae_65_lo),
    risk_sae_hi    = ifelse(AgeGroup < 16, risk_summary$sae_u65_hi,  risk_summary$sae_65_hi),
    
    risk_death_med = ifelse(AgeGroup < 16, risk_summary$death_u65_med, risk_summary$death_65_med),
    risk_death_lo  = ifelse(AgeGroup < 16, risk_summary$death_u65_lo,  risk_summary$death_65_lo),
    risk_death_hi  = ifelse(AgeGroup < 16, risk_summary$death_u65_hi,  risk_summary$death_65_hi),
    
    vax_sae_med    = tot_vacc * risk_sae_med,
    vax_sae_lo     = tot_vacc * risk_sae_lo,
    vax_sae_hi     = tot_vacc * risk_sae_hi,
    
    vax_death_med  = tot_vacc * risk_death_med,
    vax_death_lo   = tot_vacc * risk_death_lo,
    vax_death_hi   = tot_vacc * risk_death_hi
  )

summary_df <- df2 %>%
  transmute(
    scenario, AgeGroup, VE, VC,
    # vaccine attributed death (95%UIs)
    x_death_med = vax_death_med,
    x_death_lo  = vax_death_lo,
    x_death_hi  = vax_death_hi,
    # vaccine averted deaths (95%UIs)
    y_death_med = diff_fatal,
    y_death_lo  = diff_fatal_low,
    y_death_hi  = diff_fatal_hi,
    # vaccine attributed non-fatal severe cases (95%UIs)
    x_hosp_med  = vax_sae_med,
    x_hosp_lo   = vax_sae_lo,
    x_hosp_hi   = vax_sae_hi,
    # vaccine averted non-fatal severe cases (95%UIs)
    y_hosp_med  = diff_hosp,
    y_hosp_lo   = diff_hosp_lo,
    y_hosp_hi   = diff_hosp_hi
  ) %>%
  pivot_longer(
    cols = starts_with("x_") | starts_with("y_"),
    names_to = c(".value","outcome","stat"),
    names_pattern = "([xy])_(death|hosp)_(med|lo|hi)"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = c(x, y)
  ) %>%
  mutate(
    outcome = ifelse(outcome == "death", "Death", "Hospitalisation"),
    net_med = y_med - x_med
  )


summary_df2 <- summary_df %>%
  mutate(
    AgeCat = ifelse(AgeGroup < 16, "<65", "65+"),   
    AgeCat = factor(AgeCat, levels = c("<65", "65+"))
  ) %>%
  group_by(scenario, VE, VC, outcome, AgeCat) %>%
  summarise(
    x_med = sum(x_med, na.rm = TRUE),
    x_lo  = sum(x_lo, na.rm = TRUE),
    x_hi  = sum(x_hi, na.rm = TRUE),
    y_med = sum(y_med, na.rm = TRUE),
    y_lo  = sum(y_lo, na.rm = TRUE),
    y_hi  = sum(y_hi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(net_med = y_med - x_med)

x_max <- max(summary_df2$x_hi, na.rm=TRUE) * 1.1
y_max <- max(summary_df2$y_hi, na.rm=TRUE) * 1.1

bg_grid <- tidyr::expand_grid(
  x = seq(0, x_max, length.out=200),
  y = seq(0, y_max, length.out=200)
) %>%
  mutate(net_med = y - x) %>%
  crossing(
    summary_df2 %>% distinct(scenario, VE, VC)   
  )

ggplot() +
  geom_raster(data = bg_grid,
              aes(x = x, y = y, fill = net_med),
              alpha = 0.7) +
  geom_contour(data = bg_grid,
               aes(x = x, y = y, z = net_med),
               breaks = 0, color = "black",
               linetype = "dashed", linewidth = 0.6) +
  geom_errorbar(data = summary_df2,
                aes(x = x_med, ymin = y_lo, ymax = y_hi,
                    color = outcome), width = 0.05) +
  geom_errorbarh(data = summary_df2,
                 aes(y = y_med, xmin = x_lo, xmax = x_hi,
                     color = outcome), height = 0.05) +
  geom_point(data = summary_df2,
             aes(x = x_med, y = y_med,
                 color = outcome, shape = AgeCat), size = 3) +
  facet_grid(rows = vars(scenario), cols = vars(VE, VC), scales = "free") +
  scale_fill_gradient2(
    name = "Net benefit",
    low = "red", mid = "white", high = "skyblue",
    midpoint = 0
  ) +
  labs(
    x = "Vaccine-attributable",
    y = "Averted",
    color = "Outcome", shape = "Age group"
  ) +
  theme_bw(base_size = 13)
