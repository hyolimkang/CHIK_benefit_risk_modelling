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
  filter(!age_group %in% c("18-64", "65")) |>  # remove 18-64 and 65
  droplevels() 

br_representative_benefit <- br_representative_benefit %>%
  dplyr::rename(AgeCat = age_group)

x_max_global <- max(br_representative_benefit$x_hi, na.rm = TRUE) * 1.1
y_max_global <- max(br_representative_benefit$y_hi, na.rm = TRUE) * 1.1

cat("Global ranges:\n")
cat("  X: [0,", round(x_max_global, 2), "]\n")
cat("  Y: [0,", round(y_max_global, 2), "]\n")

# Create single background grid
bg_grid_unified <- expand.grid(
  x = seq(0, x_max_global, length.out = 300),
  y = seq(0, y_max_global, length.out = 300)
) %>%
  mutate(outcomes_averted = y)

## no daly
br_no_daly <- 
  br_representative_benefit %>%
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
    data = br_representative_benefit,
    aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome),
    width = x_max_global * 0.01,  # 1% of x-axis range (instead of 0)
    linewidth = 0.8,  # Slightly thicker
    alpha = 0.8  # More opaque
  ) +
  # Error bars - X-axis (horizontal)
  geom_errorbarh(
    data = br_representative_benefit,
    aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome),
    height = y_max_global * 0.01,  # 1% of y-axis range (instead of 0)
    linewidth = 0.8,  # Slightly thicker
    alpha = 0.8  # More opaque
  ) +
  
  # Points
  geom_point(
    data = br_representative_benefit,
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
    x = "Adverse events associated with vaccination (per 10,000 vaccinated individuals)",
    y = "Averted adverse events due to vaccination (per 10,000 vaccinated individuals)",
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

