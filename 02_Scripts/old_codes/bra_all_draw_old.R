## old legacy codes
# Legacy BRRAC curve helper for the all-age final block
plot_brr_ceac_outbreak_ve_legacy <- function(ceac_df,
                                             target_outcome = "DALY",
                                             setting_levels = c("Low","Moderate","High"),
                                             age_levels = c("1-11","12-17","18-64","65+"),
                                             ve_levels = NULL) {
  
  df <- ceac_df %>%
    filter(outcome == target_outcome) %>%
    mutate(
      setting = factor(setting, levels = setting_levels),
      AgeCat  = factor(AgeCat,  levels = age_levels)
    )
  
  if (!is.null(ve_levels)) {
    df <- df %>% mutate(VE_label = factor(VE_label, levels = ve_levels))
  }
  
  ggplot(df, aes(x = threshold, y = p_accept, colour = AgeCat)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.6) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    facet_grid(setting ~ VE_label) +
    labs(
      x = "Benefit-risk ratio (BRR)",
      y = "Probability (BRR > t)",
      title = paste0("BRR acceptability curve: ", target_outcome),
      colour = "Age group"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

benefit_draws_true <- dplyr::bind_rows(
  make_averted_draws_true(all_draws_daly_true,  "DALY"),
  make_averted_draws_true(all_draws_sae_true,   "SAE"),
  make_averted_draws_true(all_draws_fatal_true, "Death")
)


daly_pars_true <- lhs_sample %>%
  dplyr::mutate(draw_id = dplyr::row_number()) %>%
  dplyr::select(
    draw_id,
    le_lost_1_11, le_lost_12_17, le_lost_18_64, le_lost_65,
    dw_hosp, dw_nonhosp, dw_subac, dw_chronic,
    dur_acute, dur_nonhosp, dur_subac, dur_6m, dur_12m, dur_30m,
    acute, subac, chr6m, chr12m, chr30m
  )


# Assessment plot helper (heatmap + uncertainty bars)
create_br_plot <- function(data, target_outcome, log_min = -1, log_max = 3,
                           grid_n = 150, pad = 1.10,
                           show_prop = TRUE) {
  log_range <- seq(log_min, log_max, by = 1)
  brr_labels <- as.character(10^log_range)
  
  x_label <- dplyr::case_when(
    target_outcome == "DALY"  ~ "DALYs attributable to vaccination (per 10,000 vaccinated individuals)",
    target_outcome == "SAE"   ~ "SAEs attributable to vaccination (per 10,000 vaccinated individuals)",
    target_outcome == "Death" ~ "Deaths attributable to vaccination (per 10,000 vaccinated individuals)",
    TRUE ~ "Vaccine attributable adverse outcome (per 10,000 vaccinated individuals)"
  )
  
  y_label <- dplyr::case_when(
    target_outcome == "DALY"  ~ "DALYs averted by vaccination (per 10,000 vaccinated individuals)",
    target_outcome == "SAE"   ~ "SAEs averted by vaccination (per 10,000 vaccinated individuals)",
    target_outcome == "Death" ~ "Deaths averted by vaccination (per 10,000 vaccinated individuals)",
    TRUE ~ "Vaccine averted adverse outcome (per 10,000 vaccinated individuals)"
  )
  
  plot_data <- data %>%
    filter(outcome == target_outcome)
  
  panel_limits <- plot_data %>%
    group_by(setting, AgeCat) %>%
    summarise(
      x_max_raw = max(x_hi, na.rm = TRUE),
      y_max_raw = max(y_hi, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      x_max_raw = ifelse(is.finite(x_max_raw), x_max_raw, NA_real_),
      y_max_raw = ifelse(is.finite(y_max_raw), y_max_raw, NA_real_),
      x_min = 1e-6,
      y_min = 0,
      x_max = x_max_raw * pad,
      y_max = y_max_raw * pad
    ) %>%
    filter(!is.na(x_max), !is.na(y_max))
  
  bg_grid_specific <- panel_limits %>%
    rowwise() %>%
    do({
      panel <- .
      x_seq <- seq(panel$x_min, panel$x_max, length.out = grid_n)
      y_seq <- seq(panel$y_min, panel$y_max, length.out = grid_n)
      
      grid_df <- expand.grid(x = x_seq, y = y_seq)
      grid_df$brr <- grid_df$y / grid_df$x
      grid_df$log10_brr <- log10(grid_df$brr)
      grid_df$log10_brr[!is.finite(grid_df$log10_brr)] <- NA_real_
      grid_df$log10_brr <- pmax(pmin(grid_df$log10_brr, log_max), log_min)
      grid_df$is_fav <- !is.na(grid_df$log10_brr) & grid_df$log10_brr > 0
      grid_df$setting <- panel$setting
      grid_df$AgeCat  <- panel$AgeCat
      grid_df
    }) %>%
    ungroup() %>%
    filter(!is.na(log10_brr))
  
  panel_prop <- bg_grid_specific %>%
    group_by(setting, AgeCat) %>%
    summarise(
      prop_fav = mean(is_fav, na.rm = TRUE),
      x_min = min(x, na.rm = TRUE),
      x_max = max(x, na.rm = TRUE),
      y_min = min(y, na.rm = TRUE),
      y_max = max(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = ifelse(prop_fav < 0.005, "BRR>1: <1%",
                     sprintf("BRR>1: %.0f%%", 100 * prop_fav)),
      x_lab = x_min + 0.02 * (x_max - x_min),
      y_lab = y_max - 0.02 * (y_max - y_min)
    )
  
  p <- ggplot() +
    geom_raster(
      data = bg_grid_specific,
      aes(x = x, y = y, fill = log10_brr),
      alpha = 0.88,
      interpolate = TRUE
    ) +
    scale_fill_gradient2(
      name = "Benefit-risk ratio",
      low = "#ca0020", mid = "#f7f7f7", high = "#0571b0",
      midpoint = 0,
      limits = c(log_min, log_max),
      oob = scales::squish,
      breaks = log_range,
      labels = brr_labels,
      na.value = "white"
    ) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed",
      alpha = 0.55,
      linewidth = 0.9,
      colour = "grey35"
    ) +
    geom_errorbar(
      data = plot_data,
      aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome),
      width = 0, linewidth = 0.3
    ) +
    geom_errorbarh(
      data = plot_data,
      aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome),
      height = 0, linewidth = 0.3
    ) +
    geom_point(
      data = plot_data,
      aes(x = x_med, y = y_med, shape = VE_label, color = outcome),
      size = 2.9,
      fill = "white",
      stroke = 1.1,
      alpha = 0.95
    )
  
  if (show_prop) {
    p <- p + geom_label(
      data = panel_prop,
      aes(x = x_lab, y = y_lab, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 3,
      fill = "white",
      alpha = 0.85,
      colour = "black",
      label.size = 0.25,
      label.padding = unit(0.10, "lines")
    )
  }
  
  p +
    facet_wrap(~ setting + AgeCat, scales = "free", ncol = 4) +
    scale_color_manual(
      name = "Outcome",
      values = c("SAE" = "#1B7F1B", "Death" = "#B8860B", "DALY" = "#A23B72")
    ) +
    scale_shape_manual(
      name = "Vaccine protection mechanism",
      values = c("Disease blocking only" = 21,
                 "Disease and infection blocking" = 24)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = paste("Benefit-risk Assessment:", target_outcome),
      x = x_label,
      y = y_label
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.6),
      strip.background = element_rect(fill = "gray95", linewidth = 0.6),
      axis.title = element_text(size = 11),
      axis.text  = element_text(size = 9),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9),
      panel.spacing = unit(0.35, "lines")
    )
}


## old ver. drawing risk from travel scenario
# -------------------------------
# Step 8.4) Build traveler risk lookup
# -------------------------------
risk_lookup <- psa_df %>%
  dplyr::filter(days == "7d", state == states_to_run[1]) %>%
  dplyr::select(draw, age_group, excess_10k_sae, excess_10k_death) %>%
  dplyr::distinct()

# -------------------------------
# Step 8.5) Collapse risks into vaccine risk bands
# - u65 is shared across 1-11 / 12-17 / 18-64
# - 65+ is separate
# -------------------------------
risk_draws_true <- risk_lookup %>%
  dplyr::group_by(draw) %>%
  dplyr::summarise(
    sae_10k_u65   = max(excess_10k_sae[age_group %in% c("1-11", "12-17", "18-64")], na.rm = TRUE),
    death_10k_u65 = max(excess_10k_death[age_group %in% c("1-11", "12-17", "18-64")], na.rm = TRUE),
    
    sae_10k_65    = max(excess_10k_sae[age_group == "65+"], na.rm = TRUE),
    death_10k_65  = max(excess_10k_death[age_group == "65+"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(risk_id = dplyr::row_number())


# -------------------------------
# Step 8.7) Assign one risk draw per draw_id x risk_band
# - all under-65 age groups share the same risk draw
# -------------------------------
set.seed(1)

risk_assign_true <- benefit_base_true %>%
  dplyr::distinct(draw_id, risk_band) %>%
  dplyr::mutate(
    risk_id = sample(risk_draws_true$risk_id, size = dplyr::n(), replace = TRUE)
  ) %>%
  dplyr::left_join(risk_draws_true, by = "risk_id")

# -------------------------------
# Step 8.8) Join risk draws back to benefit rows
# -------------------------------
joint_true <- benefit_base_true %>%
  dplyr::left_join(
    risk_assign_true,
    by = c("draw_id", "risk_band")
  )

# -------------------------------
# Step 8.9) Build draw-level x/y values (benefit vs risk)
# -------------------------------
draw_level_xy_true <- joint_true %>%
  dplyr::left_join(
    tot_vacc_map_true,
    by = c("Region", "Scenario", "VE", "AgeCat")
  ) %>%
  {
    if (any(is.na(.$tot_vacc_grp))) {
      bad <- . %>%
        dplyr::filter(is.na(tot_vacc_grp)) %>%
        dplyr::distinct(Region, Scenario, AgeCat, VE) %>%
        head(20)
      stop(
        "tot_vacc_map_true join failed for some keys. Examples:\n",
        paste(capture.output(print(bad)), collapse = "\n")
      )
    }
    .
  } %>%
  dplyr::mutate(
    # disease burden per 10k vaccinated
    baseline_10k = (baseline / tot_vacc_grp) * 1e4,
    post_10k     = (post / tot_vacc_grp) * 1e4,
    y_10k        = (averted / tot_vacc_grp) * 1e4,
    
    # % reduction in disease burden
    pct_reduction = dplyr::if_else(
      baseline > 0,
      100 * averted / baseline,
      NA_real_
    ),
    
    # vaccine-induced excess risk per 10k
    sae_10k   = dplyr::if_else(risk_band == "65+", sae_10k_65, sae_10k_u65),
    death_10k = dplyr::if_else(risk_band == "65+", death_10k_65, death_10k_u65)
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    x_daly_10k = if (outcome == "DALY") {
      compute_daly_one(
        age_group      = AgeCat,
        sae_10k        = sae_10k,
        deaths_sae_10k = death_10k,
        draw_pars = list(
          le_lost_1_11  = le_lost_1_11,
          le_lost_12_17 = le_lost_12_17,
          le_lost_18_64 = le_lost_18_64,
          le_lost_65    = le_lost_65,
          dw_hosp       = dw_hosp,
          dw_nonhosp    = dw_nonhosp,
          dw_subac      = dw_subac,
          dw_chronic    = dw_chronic,
          dur_acute     = dur_acute,
          dur_nonhosp   = dur_nonhosp,
          dur_subac     = dur_subac,
          dur_6m        = dur_6m,
          dur_12m       = dur_12m,
          dur_30m       = dur_30m,
          acute         = acute,
          subac         = subac,
          chr6m         = chr6m,
          chr12m        = chr12m,
          chr30m        = chr30m
        )
      )$daly_sae
    } else {
      NA_real_
    },
    
    x_10k = dplyr::case_when(
      outcome == "SAE"   ~ sae_10k,
      outcome == "Death" ~ death_10k,
      outcome == "DALY"  ~ x_daly_10k,
      TRUE ~ NA_real_
    ),
    
    # net disease burden including vaccine harm
    net_post_10k = post_10k + x_10k,
    net_averted_10k = baseline_10k - net_post_10k,
    
    net_pct_reduction = dplyr::if_else(
      baseline_10k > 0,
      100 * net_averted_10k / baseline_10k,
      NA_real_
    )
  ) %>%
  dplyr::ungroup()


# -------------------------------
# 8) VE label (plot-friendly)
# -------------------------------
draw_level_xy_true <- draw_level_xy_true %>%
  mutate(
    VE_label = factor(VE,
                      levels = c("VE0", "VE98.9"),
                      labels = c("Disease blocking only", "Disease and infection blocking")
    )
  )%>%
  mutate(setting = unname(setting_key[Region]))

## summary
summary_pct_setting <- draw_level_xy_true %>%
  dplyr::mutate(
    setting = factor(setting, levels = c("Low", "Moderate", "High"))
  ) %>%
  dplyr::group_by(outcome, Scenario, AgeCat, VE_label, setting) %>%
  dplyr::summarise(
    baseline_med = quantile(baseline_10k, 0.50, na.rm = TRUE),
    baseline_lo  = quantile(baseline_10k, 0.025, na.rm = TRUE),
    baseline_hi  = quantile(baseline_10k, 0.975, na.rm = TRUE),
    
    post_med = quantile(post_10k, 0.50, na.rm = TRUE),
    post_lo  = quantile(post_10k, 0.025, na.rm = TRUE),
    post_hi  = quantile(post_10k, 0.975, na.rm = TRUE),
    
    pctred_med = quantile(pct_reduction, 0.50, na.rm = TRUE),
    pctred_lo  = quantile(pct_reduction, 0.025, na.rm = TRUE),
    pctred_hi  = quantile(pct_reduction, 0.975, na.rm = TRUE),
    
    netpct_med = quantile(net_pct_reduction, 0.50, na.rm = TRUE),
    netpct_lo  = quantile(net_pct_reduction, 0.025, na.rm = TRUE),
    netpct_hi  = quantile(net_pct_reduction, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


# -------------------------------
# 9) Summary for plotting (x/y quantiles)
#    Output matches your create_br_plot() expectation:
#    columns: outcome, Scenario, AgeCat, VE_label, x_lo/x_med/x_hi, y_lo/y_med/y_hi
# -------------------------------
summary_long_true <- draw_level_xy_true %>%
  group_by(outcome, Scenario, AgeCat, VE_label) %>%
  summarise(
    x_lo  = quantile(x_10k, 0.025, na.rm = TRUE),
    x_med = quantile(x_10k, 0.50,  na.rm = TRUE),
    x_hi  = quantile(x_10k, 0.975, na.rm = TRUE),
    y_lo  = quantile(y_10k, 0.025, na.rm = TRUE),
    y_med = quantile(y_10k, 0.50,  na.rm = TRUE),
    y_hi  = quantile(y_10k, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------
# 10) BRR summary + table (med + 95% UI)
# -------------------------------
brr_draw_summary_true <- draw_level_xy_true %>%
  mutate(brr = y_10k / pmax(x_10k, 1e-12)) %>%
  group_by(outcome, Scenario, AgeCat, VE_label) %>%
  summarise(
    brr_med = quantile(brr, 0.50,  na.rm = TRUE),
    brr_lo  = quantile(brr, 0.025, na.rm = TRUE),
    brr_hi  = quantile(brr, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario  = Scenario,
    age_group = AgeCat
  )

brr_table_long_true <- brr_draw_summary_true %>%
  mutate(
    brr_formatted = sprintf("%.2f [%.2f–%.2f]", brr_med, brr_lo, brr_hi),
    VE_col = as.character(VE_label)
  ) %>%
  dplyr::select(outcome, scenario, age_group, VE_col, brr_formatted)

brr_table_wide_true <- brr_table_long_true %>%
  pivot_wider(names_from = VE_col, values_from = brr_formatted) %>%
  arrange(outcome, scenario, age_group) %>%
  dplyr::select(-scenario) %>%
  dplyr::rename(`Outcome` = outcome, `Age group` = age_group)

idx_outcome_true <- table(brr_table_wide_true$`Outcome`)

kable(
  brr_table_wide_true,
  format  = "html",
  caption = "Benefit–Risk Ratio (BRR) by Outcome, Scenario, Age Group, and VE",
  align   = "l"
) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    font_size  = 12
  ) %>%
  pack_rows(index = idx_outcome_true) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2, bold = TRUE)

# -------------------------------
# 11) by setting
# -------------------------------
summary_long_setting <- draw_level_xy_true %>%
  mutate(
    setting = factor(setting, levels = c("Low", "Moderate", "High"))
  ) %>%
  group_by(outcome, Scenario, AgeCat, VE_label, setting) %>%
  summarise(
    x_lo  = quantile(x_10k, 0.025, na.rm = TRUE),
    x_med = quantile(x_10k, 0.50,  na.rm = TRUE),
    x_hi  = quantile(x_10k, 0.975, na.rm = TRUE),
    y_lo  = quantile(y_10k, 0.025, na.rm = TRUE),
    y_med = quantile(y_10k, 0.50,  na.rm = TRUE),
    y_hi  = quantile(y_10k, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(outcome, Scenario, setting, AgeCat, VE_label)
# -------------------------------
# 12) BRR summary + table (med + 95% UI)
# -------------------------------
brr_summary_setting <- draw_level_xy_true %>%
  mutate(
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    
    # Benefit-Risk Ratio
    brr = y_10k / pmax(x_10k, 1e-12),
    
    # Harm metric corresponding to each outcome
    val_caused = case_when(
      outcome == "DALY"  ~ x_daly_10k,
      outcome == "SAE"   ~ sae_10k,
      outcome == "Death" ~ death_10k,
      TRUE ~ x_10k
    ),
    
    # Benefit metric
    val_averted = y_10k
  ) %>%
  group_by(outcome, Scenario, AgeCat, VE_label, setting) %>%
  summarise(
    # BRR
    brr_med = quantile(brr, 0.50,  na.rm = TRUE),
    brr_lo  = quantile(brr, 0.025, na.rm = TRUE),
    brr_hi  = quantile(brr, 0.975, na.rm = TRUE),
    
    # Benefit
    av_med = quantile(val_averted, 0.50,  na.rm = TRUE),
    av_lo  = quantile(val_averted, 0.025, na.rm = TRUE),
    av_hi  = quantile(val_averted, 0.975, na.rm = TRUE),
    
    # Risk
    ca_med = quantile(val_caused, 0.50,  na.rm = TRUE),
    ca_lo  = quantile(val_caused, 0.025, na.rm = TRUE),
    ca_hi  = quantile(val_caused, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(outcome, Scenario, setting, AgeCat, VE_label)

brr_table_long_setting <- brr_summary_setting %>%
  mutate(
    brr_formatted = sprintf("%.2f [%.2f–%.2f]", brr_med, brr_lo, brr_hi),
    av_formatted  = sprintf("%.2f [%.2f–%.2f]", av_med,  av_lo,  av_hi),
    ca_formatted  = sprintf("%.2f [%.2f–%.2f]", ca_med,  ca_lo,  ca_hi),
    VE_col   = as.character(VE_label),
    scenario = Scenario,
    age_group = AgeCat
  ) %>%
  dplyr::select(
    outcome, scenario, setting, age_group, VE_col,
    brr_formatted, av_formatted, ca_formatted
  )

brr_table_wide_setting <- brr_table_long_setting %>%
  mutate(
    setting = factor(setting, levels = c("Low", "Moderate", "High"))
  ) %>%
  pivot_wider(
    names_from  = VE_col,
    values_from = c(av_formatted, ca_formatted, brr_formatted),
    names_glue  = "{VE_col}_{.value}"
  ) %>%
  arrange(outcome, scenario, setting, age_group) %>%
  dplyr::select(-scenario) %>%
  dplyr::rename(
    Outcome   = outcome,
    `Age group` = age_group,
    Setting   = setting
  )

idx_outcome_true <- table(brr_table_wide_setting$Outcome)

kable(
  brr_table_wide_setting,
  format  = "html",
  caption = "Integrated Benefit-Risk Assessment Table",
  align   = "l"
) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    font_size  = 11
  ) %>%
  pack_rows(index = idx_outcome_true) %>%
  column_spec(1:3, bold = TRUE) %>%
  scroll_box(width = "100%")

# ------------------------------------------------------------------------------
# Final Visualization for assessment plot
# ------------------------------------------------------------------------------

# by setting
plot_sae   <- create_br_plot(summary_long_setting, "SAE")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "D",
       caption = "Note: Background colour indicates BRR = (SAEs averted by vaccination)/(SAEs attributable to vaccination) = y/x; dashed line indicates BRR = 1 (y = x).") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1),
    plot.caption = element_text(hjust = 0, margin = ggplot2::margin(l = -8)),
    plot.caption.position = "plot",
    plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 5.5, l = 12)  
  )

plot_death <- create_br_plot(summary_long_setting, "Death")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "C",
       caption = "Note: Background colour indicates BRR = (Deaths averted by vaccination)/(Deaths attributable to vaccination) = y/x; dashed line indicates BRR = 1 (y = x).") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1),
    plot.caption = element_text(hjust = 0, margin = margin(l = -8)),
    plot.caption.position = "plot",
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 12)  
  ) 

plot_daly  <- create_br_plot(summary_long_setting, "DALY")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "B",
       caption = "Note: Background colour indicates BRR = (DALYs averted by vaccination)/(DALYs attributable to vaccination) = y/x; dashed line indicates BRR = 1 (y = x).") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1),
    plot.caption = element_text(hjust = 0, margin = margin(l = -8)),
    plot.caption.position = "plot",
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 12)  
  )
# national
#plot_sae   <- create_br_plot(summary_long_true, "SAE")
#plot_death <- create_br_plot(summary_long_true, "Death")
#plot_daly  <- create_br_plot(summary_long_true, "DALY")

plot_daly
plot_sae
plot_death


ggsave("06_Results/brr_daly_ori_setting.pdf", plot = plot_daly, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brr_death_ori_setting.pdf", plot = plot_death, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brr_sae_ori_setting.pdf", plot = plot_sae, width = 10, height = 8, device = cairo_pdf)

# ------------------------------------------------------------------------------
# Final Visualization for BRRAC
# ------------------------------------------------------------------------------

brr_long_ob <- draw_level_xy_true %>%
  mutate(
    brr = y_10k / pmax(x_10k, 1e-12),
    outcome = recode(outcome, "sae"="SAE","death"="Death","daly"="DALY"),
    setting = factor(setting, levels = c("Low","Moderate","High")),
    #AgeCat  = factor(AgeCat, levels = c("1-11","12-17","18-64","65+")),
    AgeCat  = factor(AgeCat, levels = c("18-64","65+"))
  ) %>%
  filter(is.finite(brr), brr > 0) %>%
  transmute(Region, setting, VE_label, AgeCat, outcome, brr)

brr_max <- quantile(brr_long_ob$brr, 0.999, na.rm = TRUE)  # 99.9%
brr_min <- quantile(brr_long_ob$brr, 0.001, na.rm = TRUE)  # 0.1%

lo_exp <- floor(log10(max(0.01, brr_min)))
hi_exp <- ceiling(log10(min(1e3, brr_max)))

thresholds_auto <- 10^seq(lo_exp, hi_exp, by = 0.02)
thresholds_auto

ceac_ob <- make_brr_ceac_outbreak(brr_long_ob,
                                  thresholds = thresholds_auto)


p_daly_ceac_outbreak <- plot_brr_ceac_outbreak_ve_legacy(ceac_ob, "DALY")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "E") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  )

p_death_ceac_outbreak <- plot_brr_ceac_outbreak_ve_legacy(ceac_ob, "Death")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "F") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  )

p_sae_ceac_outbreak <- plot_brr_ceac_outbreak_ve_legacy(ceac_ob, "SAE")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "G") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  )

ggsave("06_Results/brrac_daly_ori.pdf", plot = p_daly_ceac_outbreak, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brrac_death_ori.pdf", plot = p_death_ceac_outbreak, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brrac_sae_ori.pdf", plot = p_sae_ceac_outbreak, width = 10, height = 8, device = cairo_pdf)


write.csv(ceac_ob, file = "06_Results/ceac_ob.csv")


ceac_t1 <- ceac_ob %>%
  filter(abs(log10(threshold)) < 1e-12) %>%   # == threshold=1
  dplyr::select(outcome, setting, AgeCat, p_accept, threshold, VE_label)

ceac_t1_wide <- ceac_t1 %>%
  mutate(p_fmt = sprintf("%.0f%%", 100 * p_accept)) %>%
  dplyr::select(outcome, setting, AgeCat, threshold, VE_label, p_fmt) %>%
  tidyr::pivot_wider(names_from = VE_label, values_from = p_fmt)

pr_gt1 <- ceac_ob %>%
  filter(abs(log10(threshold)) < 1e-12) %>%   # threshold=1
  transmute(outcome, setting, AgeCat, VE_label,
            pr_brr_gt_1 = p_accept)

pr_gt1_wide_setting <- ceac_ob %>%
  filter(abs(log10(threshold)) < 1e-12) %>%  # threshold == 1
  transmute(
    Outcome   = outcome,
    Setting   = setting,
    `Age group` = AgeCat,
    VE_col    = VE_label,
    pr_fmt    = sprintf("%.0f%%", 100 * p_accept)  # Pr(BRR>1)라고 가정
  ) %>%
  pivot_wider(
    names_from  = VE_col,
    values_from = pr_fmt,
    names_glue  = "{VE_col} Pr(BRR>1)"
  )



# -------------------------------
# 14) Join Pr(BRR > 1) table
# -------------------------------
brr_table_wide_setting <- brr_table_wide_setting %>%
  dplyr::rename(
    `Disease blocking only (Benefit)` = `Disease blocking only_av_formatted`,
    `Disease blocking only (Risk)`    = `Disease blocking only_ca_formatted`,
    `Disease blocking only (BRR)`     = `Disease blocking only_brr_formatted`,
    
    `Disease and infection blocking (Benefit)` = `Disease and infection blocking_av_formatted`,
    `Disease and infection blocking (Risk)`    = `Disease and infection blocking_ca_formatted`,
    `Disease and infection blocking (BRR)`     = `Disease and infection blocking_brr_formatted`
  )

brr_table_wide_setting2 <- brr_table_wide_setting %>%
  left_join(
    pr_gt1_wide_setting,
    by = c("Outcome", "Setting", "Age group")
  )

# -------------------------------
# 15) Reorder columns
# -------------------------------
brr_table_wide_setting2 <- brr_table_wide_setting2 %>%
  relocate(`Disease blocking only (Benefit)`, .after = `Age group`) %>%
  relocate(`Disease blocking only (Risk)`, .after = `Disease blocking only (Benefit)`) %>%
  relocate(`Disease blocking only (BRR)`, .after = `Disease blocking only (Risk)`) %>%
  relocate(`Disease blocking only Pr(BRR>1)`, .after = `Disease blocking only (BRR)`) %>%
  relocate(`Disease and infection blocking (Benefit)`, .after = `Disease blocking only Pr(BRR>1)`) %>%
  relocate(`Disease and infection blocking (Risk)`, .after = `Disease and infection blocking (Benefit)`) %>%
  relocate(`Disease and infection blocking (BRR)`, .after = `Disease and infection blocking (Risk)`) %>%
  relocate(`Disease and infection blocking Pr(BRR>1)`, .after = `Disease and infection blocking (BRR)`)

