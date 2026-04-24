

## ORI

plot_combined_risk_probability_serostatus <- function(data) {
  
  psa_long <- dplyr::bind_rows(
    data %>%
      filter(outcome == "DALY") %>%
      dplyr::select(draw_id, setting, AgeCat, daly_10k_base, daly_10k_adj) %>%
      pivot_longer(c(daly_10k_base, daly_10k_adj),
                   names_to = "col", values_to = "val_10k") %>%
      mutate(
        outcome   = "DALY",
        mechanism = ifelse(grepl("_base", col), "Base", "Adjusted")
      ),
    
    data %>%
      filter(outcome == "SAE") %>%
      dplyr::select(draw_id, setting, AgeCat, sae_10k_base, sae_10k_adj) %>%
      pivot_longer(c(sae_10k_base, sae_10k_adj),
                   names_to = "col", values_to = "val_10k") %>%
      mutate(
        outcome   = "SAE",
        mechanism = ifelse(grepl("_base", col), "Base", "Adjusted")
      ),
    
    data %>%
      filter(outcome == "Death") %>%
      dplyr::select(draw_id, setting, AgeCat, death_10k_base, death_10k_adj) %>%
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
      setting   = factor(setting,   levels = c("Low", "Moderate", "High")),
      doses_per_case = ifelse(val_10k == 0, Inf, 10000 / val_10k)
    ) %>%
    dplyr::select(-col)
  
  # 2. X-axis range 
  outcome_ranges <- psa_long %>%
    mutate(outcome = as.character(outcome),
           setting = as.character(setting)) %>%
    group_by(outcome, setting) %>%
    summarise(
      max_range = quantile(doses_per_case[is.finite(doses_per_case)],
                           0.999, na.rm = TRUE),
      has_data  = any(is.finite(doses_per_case)),
      .groups   = "drop"
    )
  
  # 3. CDF curves
  accept_results <- psa_long %>%
    mutate(outcome   = as.character(outcome),
           mechanism = as.character(mechanism),
           setting   = as.character(setting)) %>%
    group_by(outcome, setting, mechanism, AgeCat) %>%
    group_modify(~ {
      key <- outcome_ranges %>%
        filter(outcome == unique(.y$outcome),
               setting == unique(.y$setting))
      
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
      mechanism = factor(mechanism, levels = c("Base", "Adjusted")),
      setting   = factor(setting,   levels = c("Low", "Moderate", "High")),
      # Keep 2-panel layout (Base/Adjusted) but color by setting only in Adjusted
      color_group = dplyr::if_else(mechanism == "Adjusted", as.character(setting), "Base")
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
  p <- ggplot(
    accept_results,
    aes(
      x = denom,
      y = prob_risk,
      color = color_group,
      linetype = AgeCat,
      group = interaction(AgeCat, setting, mechanism)
    )
  ) +
    
    geom_hline(yintercept = c(5, 50, 95),
               linetype = "dashed", color = "grey92", linewidth = 0.3) +
    
    geom_line(linewidth = 0.55, alpha = 0.85) +
    
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
    scale_color_manual(
      values = c(
        "Base" = "#4C78A8",
        "Low" = "#1b9e77",
        "Moderate" = "#d95f02",
        "High" = "#7570b3"
      ),
      breaks = c("Base", "Low", "Moderate", "High"),
      name = "Setting"
    ) +
    
    labs(
      x       = "Vaccinated individuals per adverse outcome",
      y       = "Probability of adverse outcome (%)",
      linetype = "Age group",
      caption = "Base panel uses one reference color. Adjusted panel colors indicate setting (Low/Moderate/High).\nDeath: no observed risk in ages 18–64."
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
  
  list(
    plot            = p,
    accept_results  = accept_results,
    psa_long        = psa_long
  )
}

brr_summarise_doses_per_case <- function(psa_long) {
  psa_long %>%
    dplyr::mutate(
      outcome   = as.character(.data$outcome),
      mechanism = as.character(.data$mechanism),
      setting   = as.character(.data$setting)
    ) %>%
    dplyr::group_by(.data$outcome, .data$setting, .data$mechanism, .data$AgeCat) %>%
    dplyr::summarise(
      n_draws                          = dplyr::n(),
      n_draws_zero_outcome_per_10k     = sum(.data$val_10k == 0, na.rm = TRUE),
      median_outcome_rate_per_10k      = stats::median(.data$val_10k, na.rm = TRUE),
      median_vaccinated_per_1_outcome  = stats::median(
        .data$doses_per_case[is.finite(.data$doses_per_case)]
      ),
      q2_5_vaccinated_per_1_outcome    = stats::quantile(
        .data$doses_per_case[is.finite(.data$doses_per_case)], 0.025, na.rm = TRUE, names = FALSE
      ),
      q50_vaccinated_per_1_outcome     = stats::quantile(
        .data$doses_per_case[is.finite(.data$doses_per_case)], 0.50, na.rm = TRUE, names = FALSE
      ),
      q97_5_vaccinated_per_1_outcome   = stats::quantile(
        .data$doses_per_case[is.finite(.data$doses_per_case)], 0.975, na.rm = TRUE, names = FALSE
      ),
      mean_vaccinated_per_1_outcome    = mean(
        .data$doses_per_case[is.finite(.data$doses_per_case)], na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      outcome   = factor(.data$outcome, levels = c("DALY", "Death", "SAE")),
      mechanism = factor(.data$mechanism, levels = c("Base", "Adjusted")),
      setting   = factor(.data$setting, levels = c("Low", "Moderate", "High"))
    )
}

brr_cdf_denom_at_prob <- function(accept_results, prob_pct = c(5, 25, 50, 75, 95)) {
  accept_results %>%
    dplyr::mutate(
      outcome   = as.character(.data$outcome),
      mechanism = as.character(.data$mechanism),
      setting   = as.character(.data$setting)
    ) %>%
    dplyr::group_by(.data$outcome, .data$setting, .data$mechanism, .data$AgeCat) %>%
    dplyr::group_modify(~ {
      ord <- order(.x$prob_risk, .x$denom)
      pr  <- as.numeric(.x$prob_risk[ord])
      dn  <- as.numeric(.x$denom[ord])
      vals <- vapply(prob_pct, function(target) {
        if (length(pr) == 0L || all(!is.finite(pr))) {
          return(NA_real_)
        }
        stats::approx(pr, dn, xout = target, rule = 2)$y
      }, numeric(1))
      dplyr::tibble(prob_psa_pct = prob_pct, vaccinated_per_1_outcome = vals)
    }) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(
      id_cols      = c(outcome, setting, mechanism, AgeCat),
      names_from   = prob_psa_pct,
      values_from  = vaccinated_per_1_outcome,
      names_prefix = "vaccinated_per_case_cdf_p"
    ) %>%
    dplyr::mutate(
      outcome   = factor(.data$outcome, levels = c("DALY", "Death", "SAE")),
      mechanism = factor(.data$mechanism, levels = c("Base", "Adjusted")),
      setting   = factor(.data$setting, levels = c("Low", "Moderate", "High"))
    )
}

# VE-separated version (additional plot)
plot_combined_risk_probability_serostatus_by_ve <- function(data) {
  psa_long <- dplyr::bind_rows(
    data %>%
      dplyr::filter(outcome == "DALY") %>%
      dplyr::select(draw_id, setting, AgeCat, VE_label, daly_10k_base, daly_10k_adj) %>%
      tidyr::pivot_longer(c(daly_10k_base, daly_10k_adj),
                          names_to = "col", values_to = "val_10k") %>%
      dplyr::mutate(
        outcome   = "DALY",
        mechanism = ifelse(grepl("_base", col), "Base", "Adjusted")
      ),
    data %>%
      dplyr::filter(outcome == "SAE") %>%
      dplyr::select(draw_id, setting, AgeCat, VE_label, sae_10k_base, sae_10k_adj) %>%
      tidyr::pivot_longer(c(sae_10k_base, sae_10k_adj),
                          names_to = "col", values_to = "val_10k") %>%
      dplyr::mutate(
        outcome   = "SAE",
        mechanism = ifelse(grepl("_base", col), "Base", "Adjusted")
      ),
    data %>%
      dplyr::filter(outcome == "Death") %>%
      dplyr::select(draw_id, setting, AgeCat, VE_label, death_10k_base, death_10k_adj) %>%
      tidyr::pivot_longer(c(death_10k_base, death_10k_adj),
                          names_to = "col", values_to = "val_10k") %>%
      dplyr::mutate(
        outcome   = "Death",
        mechanism = ifelse(grepl("_base", col), "Base", "Adjusted")
      )
  ) %>%
    dplyr::mutate(
      outcome   = factor(outcome, levels = c("DALY", "Death", "SAE")),
      mechanism = factor(mechanism, levels = c("Base", "Adjusted")),
      setting   = factor(setting, levels = c("Low", "Moderate", "High")),
      doses_per_case = ifelse(val_10k == 0, Inf, 10000 / val_10k)
    ) %>%
    dplyr::select(-col)

  outcome_ranges <- psa_long %>%
    dplyr::mutate(
      outcome = as.character(outcome),
      setting = as.character(setting),
      VE_label = as.character(VE_label)
    ) %>%
    dplyr::group_by(outcome, setting, VE_label) %>%
    dplyr::summarise(
      max_range = quantile(doses_per_case[is.finite(doses_per_case)], 0.999, na.rm = TRUE),
      has_data  = any(is.finite(doses_per_case)),
      .groups   = "drop"
    )

  accept_results <- psa_long %>%
    dplyr::mutate(
      outcome = as.character(outcome),
      mechanism = as.character(mechanism),
      setting = as.character(setting),
      VE_label = as.character(VE_label)
    ) %>%
    dplyr::group_by(outcome, setting, mechanism, AgeCat, VE_label) %>%
    dplyr::group_modify(~ {
      key <- outcome_ranges %>%
        dplyr::filter(
          outcome == unique(.y$outcome),
          setting == unique(.y$setting),
          VE_label == unique(.y$VE_label)
        )
      if (nrow(key) == 0 || !key$has_data || !any(is.finite(.x$doses_per_case))) {
        return(data.frame(denom = numeric(0), prob_risk = numeric(0)))
      }
      denom_seq <- seq(0, key$max_range, length.out = 300)
      dpc <- .x$doses_per_case
      data.frame(
        denom = denom_seq,
        prob_risk = sapply(denom_seq, function(d) mean(dpc < d) * 100)
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      outcome = factor(outcome, levels = c("DALY", "Death", "SAE")),
      mechanism = factor(mechanism, levels = c("Base", "Adjusted")),
      setting = factor(setting, levels = c("Low", "Moderate", "High")),
      color_group = dplyr::if_else(mechanism == "Adjusted", as.character(setting), "Base")
    )

  ks_formatter <- function(x) {
    dplyr::case_when(
      x >= 1e9 ~ paste0(round(x / 1e9, 1), "B"),
      x >= 1e6 ~ paste0(round(x / 1e6, 1), "M"),
      x >= 1e3 ~ paste0(round(x / 1e3, 0), "K"),
      TRUE ~ as.character(round(x, 0))
    )
  }

  p <- ggplot(
    accept_results,
    aes(
      x = denom,
      y = prob_risk,
      color = color_group,
      linetype = AgeCat,
      group = interaction(AgeCat, setting, mechanism, VE_label)
    )
  ) +
    geom_hline(yintercept = c(5, 50, 95), linetype = "dashed", color = "grey92", linewidth = 0.3) +
    geom_line(linewidth = 0.55, alpha = 0.85) +
    facet_grid(VE_label + mechanism ~ outcome, scales = "free_x") +
    coord_cartesian(ylim = c(0, 105), clip = "off") +
    scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
    scale_x_continuous(
      labels = ks_formatter,
      breaks = scales::breaks_pretty(n = 4),
      expand = c(0.05, 0)
    ) +
    scale_color_manual(
      values = c("Base" = "#4C78A8", "Low" = "#1b9e77", "Moderate" = "#d95f02", "High" = "#7570b3"),
      breaks = c("Base", "Low", "Moderate", "High"),
      name = "Setting"
    ) +
    labs(
      x = "Vaccinated individuals per adverse outcome",
      y = "Probability of adverse outcome (%)",
      linetype = "Age group",
      caption = "VE-separated panels. Base panel uses one reference color; Adjusted panel colors indicate setting."
    ) +
    theme_bw() +
    theme(
      text = element_text(family = "Calibri", size = 12),
      plot.margin = margin(10, 20, 10, 10),
      panel.grid.major = element_line(color = "grey98"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      strip.text.y = element_text(angle = 0),
      legend.position = "bottom",
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 11),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.caption = element_text(size = 9, face = "plain", margin = margin(t = 8)),
      plot.caption.position = "plot"
    )

  list(plot = p, accept_results = accept_results, psa_long = psa_long)
}

# Backward compatibility:
# - prefer pre-filtered object if it exists in the session
# - otherwise derive it from draw_level_xy_serostatus
if (!exists("draw_level_xy_serostatus_filtered", inherits = FALSE)) {
  if (!exists("draw_level_xy_serostatus", inherits = FALSE)) {
    stop("Neither 'draw_level_xy_serostatus_filtered' nor 'draw_level_xy_serostatus' exists. Run brazil_all_draws_ori_v2.R first.")
  }
  draw_level_xy_serostatus_filtered <- draw_level_xy_serostatus %>%
    dplyr::filter(
      AgeCat %in% c("18-64", "65+"),
      RR_seropos == 0.0
    )
}

out_risk_curves <- plot_combined_risk_probability_serostatus(draw_level_xy_serostatus_filtered)
p_serostatus <- out_risk_curves$plot
print(p_serostatus)

ggsave("06_Results/brr_risk_curves.pdf", plot = p_serostatus, width = 12, height = 6, device = cairo_pdf)

out_risk_curves_by_ve <- plot_combined_risk_probability_serostatus_by_ve(draw_level_xy_serostatus_filtered)
p_serostatus_by_ve <- out_risk_curves_by_ve$plot
print(p_serostatus_by_ve)
ggsave("06_Results/brr_risk_curves_by_ve.pdf", plot = p_serostatus_by_ve, width = 12, height = 8, device = cairo_pdf)

summary_doses_per_case <- brr_summarise_doses_per_case(out_risk_curves$psa_long)
cdf_denom_at_prob      <- brr_cdf_denom_at_prob(out_risk_curves$accept_results)

# Diagnostics: q_seroneg mix and adjusted risk by setting/age/VE
q_seroneg_summary <- draw_level_xy_serostatus_filtered %>%
  dplyr::select(draw_id, setting, AgeCat, VE_label, q_seroneg_vacc) %>%
  dplyr::distinct() %>%
  dplyr::group_by(setting, AgeCat, VE_label) %>%
  dplyr::summarise(
    n = dplyr::n(),
    q_seroneg_med = stats::median(q_seroneg_vacc, na.rm = TRUE),
    q_seroneg_p2_5 = stats::quantile(q_seroneg_vacc, 0.025, na.rm = TRUE, names = FALSE),
    q_seroneg_p97_5 = stats::quantile(q_seroneg_vacc, 0.975, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(VE_label, AgeCat, setting)

adjusted_risk_summary <- draw_level_xy_serostatus_filtered %>%
  dplyr::transmute(
    setting, AgeCat, VE_label, outcome,
    risk_adj_10k = x_10k_adj
  ) %>%
  dplyr::filter(is.finite(risk_adj_10k)) %>%
  dplyr::group_by(outcome, setting, AgeCat, VE_label) %>%
  dplyr::summarise(
    n = dplyr::n(),
    risk_adj_med = stats::median(risk_adj_10k, na.rm = TRUE),
    risk_adj_p2_5 = stats::quantile(risk_adj_10k, 0.025, na.rm = TRUE, names = FALSE),
    risk_adj_p97_5 = stats::quantile(risk_adj_10k, 0.975, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(outcome, VE_label, AgeCat, setting)

xlsx_path <- "06_Results/brr_risk_curves_data.xlsx"
writexl::write_xlsx(
  list(
    README = tibble::tibble(
      description = c(
        "CDF_curves: x = vaccinated individuals per one outcome (denom); y = prob_risk = % of PSA draws with doses_per_case < denom (see script).",
        "PSA_draws: one row per draw_id x AgeCat x outcome x mechanism; val_10k = outcome count per 10,000 vaccinated; doses_per_case = 10000/val_10k (Inf if val_10k==0).",
        "Summary_vaccinated_per_case: PSA quantiles of doses_per_case (vaccinated persons per one outcome event).",
        "Summary_cdf_thresholds: vaccinated-per-case on the plotted CDF at which prob_risk reaches 5,25,...,95 (% of draws).",
        "Q_seroneg_summary: distribution of vaccinated-cohort seronegative fraction (q_seroneg_vacc) by setting x age group x VE.",
        "Adjusted_risk_summary: distribution of adjusted risk per 10,000 (x_10k_adj) by outcome x setting x age group x VE."
      )
    ),
    CDF_curves                 = out_risk_curves$accept_results,
    PSA_draws                  = out_risk_curves$psa_long,
    Summary_vaccinated_per_case = summary_doses_per_case,
    Summary_cdf_thresholds     = cdf_denom_at_prob,
    Q_seroneg_summary          = q_seroneg_summary,
    Adjusted_risk_summary      = adjusted_risk_summary
  ),
  path = xlsx_path
)
cat("Saved:", xlsx_path, "\n")
