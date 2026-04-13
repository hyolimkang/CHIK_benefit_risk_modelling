psa_long <- dplyr::bind_rows(
  draw_level_xy_serostatus_filtered %>%
    filter(outcome == "DALY") %>%
    dplyr::select(draw_id, AgeCat, daly_10k_base, daly_10k_adj) %>%
    pivot_longer(c(daly_10k_base, daly_10k_adj),
                 names_to = "col", values_to = "val_10k") %>%
    mutate(
      outcome   = "DALY",
      mechanism = ifelse(grepl("_base", col), "Base", "Adjusted")
    ),
  
  draw_level_xy_serostatus_filtered %>%
    filter(outcome == "SAE") %>%
    dplyr::select(draw_id, AgeCat, sae_10k_base, sae_10k_adj) %>%
    pivot_longer(c(sae_10k_base, sae_10k_adj),
                 names_to = "col", values_to = "val_10k") %>%
    mutate(
      outcome   = "SAE",
      mechanism = ifelse(grepl("_base", col), "Base", "Adjusted")
    ),
  
  draw_level_xy_serostatus_filtered %>%
    filter(outcome == "Death") %>%
    dplyr::select(draw_id, AgeCat, death_10k_base, death_10k_adj) %>%
    pivot_longer(c(death_10k_base, death_10k_adj),
                 names_to = "col", values_to = "val_10k") %>%
    mutate(
      outcome   = "Death",
      mechanism = ifelse(grepl("_base", col), "Base", "Adjusted")
    )
) %>%
  mutate(
    outcome   = factor(outcome,   levels = c("DALY", "Death", "SAE")),
    mechanism = factor(mechanism, levels = c("Base", "Adjusted")),
    doses_per_case = ifelse(val_10k == 0, Inf, 10000 / val_10k)
  ) %>%
  dplyr::select(-col)

# 2. X-axis range — outcome × mechanism 패널별, 유한값만
outcome_ranges <- psa_long %>%
  mutate(outcome = as.character(outcome),
         mechanism = as.character(mechanism)) %>%
  group_by(outcome, mechanism) %>%
  summarise(
    max_range = quantile(doses_per_case[is.finite(doses_per_case)],
                         0.999, na.rm = TRUE),
    has_data  = any(is.finite(doses_per_case)),
    .groups   = "drop"
  )

# 3. CDF curves
accept_results <- psa_long %>%
  mutate(outcome   = as.character(outcome),
         mechanism = as.character(mechanism)) %>%
  group_by(outcome, mechanism, AgeCat) %>%
  group_modify(~ {
    key <- outcome_ranges %>%
      filter(outcome   == unique(.y$outcome),
             mechanism == unique(.y$mechanism))
    
    if (nrow(key) == 0 || !key$has_data || !any(is.finite(.x$doses_per_case))) {
      return(data.frame(denom = numeric(0), prob_risk = numeric(0)))
    }
    
    denom_seq <- seq(0, key$max_range, length.out = 300)
    dpc       <- .x$doses_per_case
    
    data.frame(
      denom     = denom_seq,
      prob_risk = sapply(denom_seq, function(d) mean(dpc < d) * 100)
    )
  }) %>%
  ungroup() %>%
  mutate(
    outcome   = factor(outcome,   levels = c("DALY", "Death", "SAE")),
    mechanism = factor(mechanism, levels = c("Base", "Adjusted"))
  )

# val_10k 원본값과 doses_per_case 직접 비교
psa_long %>%
  mutate(outcome = as.character(outcome),
         mechanism = as.character(mechanism)) %>%
  filter(is.finite(doses_per_case)) %>%
  group_by(outcome, mechanism, AgeCat) %>%
  summarise(
    # val_10k: 10,000명당 이상반응 건수
    val_med  = median(val_10k, na.rm = TRUE),
    val_p001 = quantile(val_10k, 0.001, na.rm = TRUE),  # 가장 작은 값
    val_p999 = quantile(val_10k, 0.999, na.rm = TRUE),  # 가장 큰 값
    # doses_per_case = 10000 / val_10k
    dpc_med  = median(doses_per_case, na.rm = TRUE),
    dpc_p001 = quantile(doses_per_case, 0.001, na.rm = TRUE),
    dpc_p999 = quantile(doses_per_case, 0.999, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 수동 검증: dpc_med ≈ 10000 / val_med 인지 확인
  mutate(check = round(10000 / val_med, 0),
         match = abs(dpc_med - check) < 1) %>%
  print(n = Inf)



### old codes
plot_combined_risk_probability_serostatus <- function(draw_level_xy_serostatus) {
  
  # 1. Data Transformation
  psa_long <- draw_level_xy_serostatus %>%
    dplyr::select(draw_id, AgeCat, 
                  sae_10k_base, death_10k_base, daly_10k_base,
                  sae_10k_adj,  death_10k_adj,  daly_10k_adj) %>%
    pivot_longer(
      cols = c(sae_10k_base, death_10k_base, daly_10k_base,
               sae_10k_adj,  death_10k_adj,  daly_10k_adj),
      names_to  = "outcome_raw",
      values_to = "val_10k"
    ) %>%
    mutate(
      outcome = case_when(
        grepl("death", outcome_raw) ~ "Death",
        grepl("sae",   outcome_raw) ~ "SAE",
        grepl("daly",  outcome_raw) ~ "DALY"
      ),
      mechanism = case_when(
        grepl("_base", outcome_raw) ~ "Base",
        grepl("_adj",  outcome_raw) ~ "Adjusted"
      ),
      doses_per_case = 10000 / pmax(1e-9, val_10k)
    )
  
  # 2. X-axis range by outcome
  outcome_ranges <- psa_long %>%
    group_by(outcome) %>%
    summarise(max_range = quantile(doses_per_case, 0.99, na.rm = TRUE))
  
  # 3. Risk probability curves
  accept_results <- psa_long %>%
    group_by(outcome, AgeCat, mechanism) %>%
    group_modify(~ {
      current_max <- outcome_ranges$max_range[outcome_ranges$outcome == unique(.y$outcome)]
      denom_seq   <- seq(0, current_max, length.out = 200)
      
      data.frame(denom = denom_seq) %>%
        rowwise() %>%
        mutate(prob_risk = mean(.x$doses_per_case < denom, na.rm = TRUE) * 100)
    }) %>%
    ungroup()
  
  # 4. Formatter
  ks_formatter <- function(x) {
    case_when(
      x >= 1e6 ~ paste0(round(x/1e6, 1), "M"),
      x >= 1e3 ~ paste0(round(x/1e3, 0), "K"),
      TRUE     ~ as.character(round(x, 0))
    )
  }
  
  # 5. Threshold points (100%)
  threshold_points <- accept_results %>%
    group_by(outcome, AgeCat, mechanism) %>%
    filter(prob_risk >= (max(prob_risk) * 0.999)) %>%
    filter(denom == min(denom)) %>%
    ungroup() %>%
    mutate(label_text = ks_formatter(denom)) %>%
    group_by(outcome, label_text) %>%
    mutate(is_duplicate = row_number() > 1) %>%
    ungroup() %>%
    group_by(outcome) %>%
    arrange(denom) %>%
    mutate(
      y_pos     = seq(4, by = 10, length.out = n()),
      hjust_val = rep(c(-0.2, 0.5, 1.2), length.out = n())
    ) %>%
    ungroup()
  
  # 6. Plot
  # AgeCat → color, mechanism → linetype
  p <- ggplot(accept_results,
              aes(x = denom, y = prob_risk,
                  color    = AgeCat,
                  linetype = mechanism)) +
    
    geom_hline(yintercept = c(5, 50, 95),
               linetype = "dashed", color = "grey92", linewidth = 0.3) +
    
    geom_vline(data = threshold_points,
               aes(xintercept = denom, color = AgeCat, linetype = mechanism),
               linewidth = 0.5, alpha = 0.6) +
    
    geom_line(linewidth = 0.7, alpha = 0.85) +
    
    geom_point(data = threshold_points,
               aes(x = denom, y = prob_risk),
               size = 1.8, shape = 21, fill = "white", stroke = 0.8) +
    
    geom_text(data = filter(threshold_points, !is_duplicate),
              aes(x = denom, y = y_pos, label = label_text, color = AgeCat),
              vjust = 0, size = 2.8, fontface = "bold",
              family = "Calibri", show.legend = FALSE) +
    
    facet_wrap(~outcome, scales = "free_x", ncol = 3) +
    
    coord_cartesian(ylim = c(0, 105), clip = "off") +
    scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
    scale_x_continuous(labels = ks_formatter, expand = c(0.05, 0)) +
    scale_linetype_manual(values = c("Base" = "solid", "Adjusted" = "dashed")) +
    
    labs(
      x        = "Vaccinated individuals per adverse outcome",
      y        = "Probability of adverse outcome (%)",
      color    = "Age group",
      linetype = "Risk mechanism",
      caption  = "Base: unadjusted vaccine risk; Adjusted: serostatus-weighted risk"
    ) +
    theme_bw() +
    theme(
      text              = element_text(family = "Calibri", size = 12),
      plot.margin       = margin(10, 35, 10, 10),
      panel.grid.major  = element_line(color = "grey98"),
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(face = "bold", size = 11),
      legend.position   = "bottom",
      axis.title        = element_text(face = "bold"),
      axis.text         = element_text(size = 12),
      plot.caption      = element_text(size = 10, face = "plain", margin = margin(t = 8)),
      plot.caption.position = "plot"
    )
  
  return(p)
}

p_serostatus <- plot_combined_risk_probability_serostatus(draw_level_xy_serostatus)
print(p_serostatus)

###
## Travel (seronegative)
plot_combined_risk_probability_ascending <- function(psa_df) {
  
  # 1. Data Transformation and Calculation of 'Doses per Case'
  psa_long <- psa_df %>%
    dplyr::select(draw, age_group, excess_10k_death, excess_10k_sae, daly_sae) %>%
    pivot_longer(
      cols = c(excess_10k_death, excess_10k_sae, daly_sae),
      names_to = "outcome_raw",
      values_to = "val_10k"
    ) %>%
    mutate(outcome = case_when(
      outcome_raw == "excess_10k_death" ~ "Death",
      outcome_raw == "excess_10k_sae"   ~ "SAE",
      outcome_raw == "daly_sae"         ~ "DALY",
      TRUE ~ outcome_raw
    )) %>%
    # Calculate actual safety performance: number of doses per 1 event
    mutate(doses_per_case = 10000 / pmax(1e-9, val_10k))
  
  # 2. Define Data-driven X-axis Range (Based on 99th percentile to avoid outliers)
  outcome_ranges <- psa_long %>%
    group_by(outcome) %>%
    summarise(max_range = quantile(doses_per_case, 0.99, na.rm = TRUE))
  
  # 3. Calculate Risk Probability
  accept_results <- psa_long %>%
    group_by(outcome, age_group) %>%
    group_modify(~ {
      current_max <- outcome_ranges$max_range[outcome_ranges$outcome == unique(.y$outcome)]
      denom_seq <- seq(0, current_max, length.out = 200)
      
      data.frame(denom = denom_seq) %>%
        rowwise() %>%
        mutate(prob_risk = mean(.x$doses_per_case < denom, na.rm = TRUE) * 100)
    }) %>%
    ungroup()
  
  # 4. Custom Formatter for K (Thousands) and M (Millions)
  ks_formatter <- function(x) {
    case_when(
      x >= 1e6 ~ paste0(round(x/1e6, 1), "M"),
      x >= 1e3 ~ paste0(round(x/1e3, 0), "K"),
      TRUE ~ as.character(round(x, 0))
    )
  }
  
  # 5. Extract 100% Points and Handle Label Overlapping/Duplicates
  threshold_points <- accept_results %>%
    group_by(outcome, age_group) %>%
    filter(prob_risk >= (max(prob_risk) * 0.999)) %>% 
    filter(denom == min(denom)) %>% 
    ungroup() %>%
    mutate(label_text = ks_formatter(denom)) %>%
    # Remove duplicate labels if multiple groups have the same threshold value
    group_by(outcome, label_text) %>%
    mutate(is_duplicate = row_number() > 1) %>%
    ungroup() %>%
    # Tiered vertical positioning (y_pos) for labels to prevent overlapping
    group_by(outcome) %>%
    arrange(denom) %>%
    mutate(
      y_pos = seq(4, by = 10, length.out = n()), 
      hjust_val = rep(c(-0.2, 0.5, 1.2), length.out = n())
    ) %>%
    ungroup()
  
  # 6. Visualization
  p <- ggplot(accept_results, aes(x = denom, y = prob_risk, color = age_group)) +
    # Background reference lines
    geom_hline(yintercept = c(5, 50, 95), linetype = "dashed", color = "grey92", linewidth = 0.3) +
    
    # Vertical reference lines (Dotted) for all groups
    geom_vline(data = threshold_points, aes(xintercept = denom, color = age_group), 
               linetype = "dotted", linewidth = 0.5, alpha = 0.8) +
    
    # Main curves (Slightly thickened as requested: 0.95)
    geom_line(linewidth = 0.7, alpha = 0.85) +
    
    # Points at 100% threshold
    geom_point(data = threshold_points, aes(x = denom, y = prob_risk), 
               size = 1.8, shape = 21, fill = "white", stroke = 0.8) +
    
    # Threshold text labels (Only non-duplicates, color-coded)
    geom_text(data = filter(threshold_points, !is_duplicate), 
              aes(x = denom, y = y_pos, label = label_text, color = age_group),
              vjust = 0, size = 2.8, fontface = "bold", 
              family = "Calibri", show.legend = FALSE) +
    
    # Faceting by outcome
    facet_wrap(~outcome, scales = "free_x", ncol = 3) +
    
    # Formatting Axes
    coord_cartesian(ylim = c(0, 105), clip = "off") + 
    scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
    scale_x_continuous(labels = ks_formatter, expand = c(0.05, 0)) +
    
    # Labels and Theme
    labs(
      x = "Vaccinated individuals per adverse outcome",
      y = "Probability of adverse outcome (%)",
      color = "Age Group",
      caption = "Note: Threshold values represent the point where risk reaches 100%"
    ) +
    theme_bw() + # Restoring default theme_bw (Grey strip backgrounds)
    theme(
      text = element_text(family = "Calibri", size = 12),
      plot.margin = margin(10, 35, 10, 10),
      panel.grid.major = element_line(color = "grey98"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(size = 12),
    ) + theme(
      plot.caption = element_text(size = 12, face = "plain", margin = margin(t = 8)),
      plot.caption.position = "plot"
    )
  
  return(p)
}

p_combined <- plot_combined_risk_probability_ascending(psa_df)
print(p_combined)

ggsave("06_Results/brr_risk_curves.pdf", plot = p_combined, width = 12, height = 6, device = cairo_pdf)

