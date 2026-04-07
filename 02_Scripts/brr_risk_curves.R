
plot_combined_risk_acceptability_linear <- function(psa_df) {
  
  psa_long <- psa_df %>%
    dplyr::select(draw, age_group, excess_10k_death, excess_10k_sae, daly_sae) %>%
    pivot_longer(
      cols = c(excess_10k_death, excess_10k_sae, daly_sae),
      names_to = "outcome_raw",
      values_to = "val_to_plot"
    ) %>%
    mutate(outcome = case_when(
      outcome_raw == "excess_10k_death" ~ "Death",
      outcome_raw == "excess_10k_sae"   ~ "SAE",
      outcome_raw == "daly_sae"         ~ "DALY",
      TRUE ~ outcome_raw
    )) %>%
    mutate(val_to_plot = pmax(0, val_to_plot, na.rm = TRUE))
  
  accept_results <- psa_long %>%
    group_by(outcome, age_group) %>%
    nest() %>%
    mutate(results = map2(data, outcome, function(df_sub, out_name) {
      
      max_val <- max(df_sub$val_to_plot, na.rm = TRUE)
      q99_val <- quantile(df_sub$val_to_plot, 0.99, na.rm = TRUE)
      
      target_upper <- if(q99_val > 0) (10000 / (q99_val / 5)) else 100000
      
      upper_limit <- case_when(
        out_name == "Death" ~ max(200000, target_upper), 
        out_name == "SAE"   ~ max(20000, target_upper),  
        TRUE                ~ max(10000, target_upper)   
      )
      
      denom_seq <- seq(100, upper_limit, length.out = 120)
      
      map_dfr(denom_seq, function(d) {
        threshold_val <- 10000 / d
        data.frame(
          denom = d,
          prob_acceptable = mean(df_sub$val_to_plot <= threshold_val) * 100
        )
      })
    })) %>%
    dplyr::select(-data) %>%
    unnest(results)
  
  p <- ggplot(accept_results, aes(x = denom, y = prob_acceptable, color = age_group)) +
    geom_line(linewidth = 0.75, na.rm = TRUE) +
    facet_wrap(~outcome, scales = "free_x", ncol = 3) + 
    scale_x_continuous(labels = label_comma(), expand = c(0.02, 0)) + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 2)) +
    labs(
      x = "Acceptable Risk Threshold (1 in X doses)",
      y = "% of probabilistic simulations satisfying threshold",
      color = "Age Group"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold"),
      text = element_text(family = "Calibri"),
      axis.text  = element_text(size = 11)
    )
  
  return(p)
}

p_combined <- plot_combined_risk_acceptability_linear(psa_df)
print(p_combined)




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


## summ
risk_table_summary <- psa_df %>%
  dplyr::select(draw, age_group, excess_10k_death, excess_10k_sae, daly_sae) %>%
  tidyr::pivot_longer(
    cols = c(excess_10k_death, excess_10k_sae, daly_sae),
    names_to = "outcome_raw",
    values_to = "val_10k"
  ) %>%
  dplyr::mutate(
    outcome = dplyr::case_when(
      outcome_raw == "excess_10k_death" ~ "Death",
      outcome_raw == "excess_10k_sae"   ~ "SAE",
      outcome_raw == "daly_sae"         ~ "DALY",
      TRUE ~ outcome_raw
    ),
    vaccinated_per_1_excess = 10000 / pmax(val_10k, 1e-9)
  ) %>%
  dplyr::group_by(outcome, age_group) %>%
  dplyr::summarise(
    n_draws = dplyr::n(),
    excess_per_10k_med = quantile(val_10k, 0.50, na.rm = TRUE),
    excess_per_10k_lo  = quantile(val_10k, 0.025, na.rm = TRUE),
    excess_per_10k_hi  = quantile(val_10k, 0.975, na.rm = TRUE),
    
    vacc_per_1_med = quantile(vaccinated_per_1_excess, 0.50, na.rm = TRUE),
    vacc_per_1_lo  = quantile(vaccinated_per_1_excess, 0.025, na.rm = TRUE),
    vacc_per_1_hi  = quantile(vaccinated_per_1_excess, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


fmt_num <- function(x) {
  format(round(x, 0), big.mark = ",", scientific = FALSE)
}

risk_table_pretty <- risk_table_summary %>%
  dplyr::mutate(
    excess_per_10k = sprintf(
      "%.2f [%.2f–%.2f]",
      excess_per_10k_med, excess_per_10k_lo, excess_per_10k_hi
    ),
    one_excess_per_vaccinated = paste0(
      fmt_num(vacc_per_1_med), " [",
      fmt_num(vacc_per_1_lo), "–",
      fmt_num(vacc_per_1_hi), "]"
    )
  ) %>%
  dplyr::select(
    outcome,
    age_group,
    excess_per_10k,
    one_excess_per_vaccinated
  )


write_xlsx(risk_table_pretty, "06_Results/risk_tbl.xlsx")