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

