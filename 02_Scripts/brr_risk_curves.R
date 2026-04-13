

## ORI

plot_combined_risk_probability_serostatus <- function(data) {
  
  psa_long <- dplyr::bind_rows(
    data %>%
      filter(outcome == "DALY") %>%
      dplyr::select(draw_id, AgeCat, daly_10k_base, daly_10k_adj) %>%
      pivot_longer(c(daly_10k_base, daly_10k_adj),
                   names_to = "col", values_to = "val_10k") %>%
      mutate(
        outcome   = "DALY",
        mechanism = ifelse(grepl("_base", col), "Base", "Adjusted")
      ),
    
    data %>%
      filter(outcome == "SAE") %>%
      dplyr::select(draw_id, AgeCat, sae_10k_base, sae_10k_adj) %>%
      pivot_longer(c(sae_10k_base, sae_10k_adj),
                   names_to = "col", values_to = "val_10k") %>%
      mutate(
        outcome   = "SAE",
        mechanism = ifelse(grepl("_base", col), "Base", "Adjusted")
      ),
    
    data %>%
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
  
  # 2. X-axis range 
  outcome_ranges <- psa_long %>%
    mutate(outcome = as.character(outcome)) %>%
    group_by(outcome) %>%
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
        filter(outcome == unique(.y$outcome))
      
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
  
  # 4. Formatter
  ks_formatter <- function(x) {
    dplyr::case_when(
      x >= 1e9 ~ paste0(round(x / 1e9, 1), "B"),
      x >= 1e6 ~ paste0(round(x / 1e6, 1), "M"),
      x >= 1e3 ~ paste0(round(x / 1e3, 0), "K"),
      TRUE     ~ as.character(round(x, 0))
    )
  }
  
  # 5. Annotation — Death 18-64 양쪽 패널
  annot_df <- data.frame(
    outcome   = factor(c("Death","Death"), levels = c("DALY","Death","SAE")),
    mechanism = factor(c("Base","Adjusted"), levels = c("Base","Adjusted"))
    #label     = "No observed risk\nin ages 18–64"
  )
  
  # 6. Plot
  p <- ggplot(accept_results,
              aes(x = denom, y = prob_risk, color = AgeCat)) +
    
    geom_hline(yintercept = c(5, 50, 95),
               linetype = "dashed", color = "grey92", linewidth = 0.3) +
    
    geom_line(linewidth = 0.7, alpha = 0.85) +
    
    #geom_text(
    #  data        = annot_df,
     # aes(x = Inf, y = 50, label = label),
    #  inherit.aes = FALSE,
    #  hjust = 1.1, vjust = 0.5,
    #  size = 3, color = "grey40", fontface = "italic"
   # ) +
    
    facet_grid(mechanism ~ outcome, scales = "free_x") +
    
    coord_cartesian(ylim = c(0, 105), clip = "off") +
    scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
    scale_x_continuous(
      labels = ks_formatter,
      breaks = scales::breaks_pretty(n = 4),
      expand = c(0.05, 0)
    ) +
    
    labs(
      x       = "Vaccinated individuals per adverse outcome",
      y       = "Probability of adverse outcome (%)",
      color   = "Age group",
      caption = "Base: unadjusted vaccine risk; Adjusted: serostatus-weighted risk.\nDeath: no observed risk in ages 18–64."
    ) +
    theme_bw() +
    theme(
      text                  = element_text(family = "Calibri", size = 12),
      plot.margin           = margin(10, 20, 10, 10),
      panel.grid.major      = element_line(color = "grey98"),
      panel.grid.minor      = element_blank(),
      strip.text            = element_text(face = "bold", size = 11),
      strip.text.y          = element_text(angle = 0),
      legend.position       = "bottom",
      axis.title            = element_text(face = "bold"),
      axis.text             = element_text(size = 11),
      axis.text.x           = element_text(angle = 30, hjust = 1),
      plot.caption          = element_text(size = 9, face = "plain",
                                           margin = margin(t = 8)),
      plot.caption.position = "plot"
    )
  
  return(p)
}

p_serostatus <- plot_combined_risk_probability_serostatus(draw_level_xy_serostatus_filtered)
print(p_serostatus)

ggsave("06_Results/brr_risk_curves.pdf", plot = p_serostatus, width = 12, height = 6, device = cairo_pdf)
