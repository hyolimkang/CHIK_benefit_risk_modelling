library(tidyverse)

# --- PRE-STEP: Synchronize Scales & Breaks ---
actual_brr_values <- br_representative_benefit %>%
  mutate(brr = y_med / (x_med + 1e-6)) %>%
  pull(brr)

log_min <- floor(log10(min(actual_brr_values, na.rm = TRUE)))
log_max <- ceiling(log10(max(actual_brr_values, na.rm = TRUE)))
log_range <- seq(log_min, log_max, by = 1)
brr_labels <- paste0("BRR = ", 10^log_range)

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

# --- STEP 2: Calculate Smart Limits (With Negative Buffer) ---
row_limits <- br_representative_benefit %>%
  group_by(ar_category) %>%
  summarise(
    x_max = max(x_hi, na.rm = TRUE) * 1.5,
    y_max = max(y_hi, na.rm = TRUE) * 1.5,
    
    x_min = -max(x_hi, na.rm = TRUE) * 0.03,
    y_min = -max(y_hi, na.rm = TRUE) * 0.03,
    .groups = "drop"
  )

# --- STEP 3: Generate Grid for geom_raster ---
bg_grid_raster <- br_representative_benefit %>%
  distinct(ar_category, days) %>%
  left_join(row_limits, by = "ar_category") %>%
  rowwise() %>%
  mutate(grid = list(
    expand.grid(
      x = seq(x_min, x_max, length.out = 150),
      y = seq(y_min, y_max, length.out = 150)
    ) %>%
      mutate(
        calc_x = pmax(x, 1e-6),
        calc_y = pmax(y, 1e-6),
        brr = calc_y / calc_x,
        log10_brr = log10(brr)
      )
  )) %>%
  ungroup() %>%
  tidyr::unnest(cols = grid)

# --- STEP 4: Final Plotting ---
plot_final <- ggplot() +
  geom_raster(
    data = bg_grid_raster,
    aes(x = x, y = y, fill = log10_brr),
    interpolate = TRUE,
    hjust = 0, vjust = 0 
  ) +
  
  scale_fill_gradient2(
    name = "Benefit–risk ratio",
    low = "#ca0020", mid = "#f7f7f7", high = "#0571b0",
    midpoint = 0,
    limits = c(log_min, log_max),
    breaks = log_range,
    labels = brr_labels,
    oob = scales::squish
  ) +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.9, alpha = 0.8) +
  
  geom_errorbar(data = br_representative_benefit, 
                aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome), 
                width = 0, linewidth = 0.6) +
  geom_errorbarh(data = br_representative_benefit, 
                 aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome), 
                 height = 0, linewidth = 0.6) +
  geom_point(data = br_representative_benefit, 
             aes(x = x_med, y = y_med, shape = AgeCat, color = outcome), 
             size = 2.5, stroke = 0.7, fill = "white") +
  
  facet_grid(
    rows = vars(ar_category), 
    cols = vars(days), 
    scales = "free_y"
  ) +
  
  scale_color_manual(values = c("SAE" = "#1B7F1B", "Death" = "#B8860B"), name = "Outcome type") +
  scale_shape_manual(values = c("1-11"=21, "12-17"=22, "18-59"=23, "60"=24), name = "Age group") +
  
  coord_cartesian(expand = FALSE, clip = "on") +
  
  labs(
    x = "Excess adverse outcomes (per 10,000 vaccinated)",
    y = "Infection-related adverse outcomes averted (per 10,000 vaccinated)",
    title = "Benefit-Risk Space: Traveller vaccination"
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95"),
    panel.spacing = unit(0.8, "lines"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 10, hjust = 0)
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
ggsave("06_Results/brr_final_plot1.pdf", plot = combined_plot, width = 10, height = 12)
