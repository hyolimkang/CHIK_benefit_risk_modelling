group_map <- c(
  "p_sae_vacc_u65"    = "Vaccine Risk",
  "p_sae_vacc_65"     = "Vaccine Risk",
  "p_death_vacc_u65"  = "Vaccine Risk",
  "p_death_vacc_65"   = "Vaccine Risk",
  "p_sae_nat_11"      = "Natural Disease Risk",
  "p_sae_nat_17"      = "Natural Disease Risk",
  "p_sae_nat_64"      = "Natural Disease Risk",
  "p_sae_nat_65"      = "Natural Disease Risk",
  "p_death_nat_11"    = "Natural Disease Risk",
  "p_death_nat_17"    = "Natural Disease Risk",
  "p_death_nat_64"    = "Natural Disease Risk",
  "p_death_nat_65"    = "Natural Disease Risk",
  "ve"                = "Vaccine Effectiveness",
  "ar_small"          = "Attack Rate / Exposure",
  "ar_med"            = "Attack Rate / Exposure",
  "ar_large"          = "Attack Rate / Exposure",
  "epi_months"        = "Attack Rate / Exposure",
  "trav_7d"           = "Attack Rate / Exposure",
  "trav_14d"          = "Attack Rate / Exposure",
  "trav_30d"          = "Attack Rate / Exposure",
  "trav_90d"          = "Attack Rate / Exposure",
  "symp_asia"         = "Symptomatic / Severity",
  "symp_africa"       = "Symptomatic / Severity",
  "symp_america"      = "Symptomatic / Severity",
  "symp_overall"      = "Symptomatic / Severity",
  "fatal_hosp"        = "Symptomatic / Severity",
  "hosp"              = "Symptomatic / Severity",
  "lt"                = "Symptomatic / Severity",
  "le_lost_1_11"      = "Life Expectancy",
  "le_lost_12_17"     = "Life Expectancy",
  "le_lost_18_64"     = "Life Expectancy",
  "le_lost_65"        = "Life Expectancy",
  "dw_chronic"        = "DALY Components",
  "dur_chronic"       = "DALY Components",
  "dw_hosp"           = "DALY Components",
  "dur_acute"         = "DALY Components",
  "dw_nonhosp"        = "DALY Components",
  "dur_nonhosp"       = "DALY Components",
  "acute"             = "DALY Components",
  "subac"             = "DALY Components",
  "chr6m"             = "DALY Components",
  "chr12m"            = "DALY Components",
  "chr30m"            = "DALY Components",
  "dw_chronic_mild"   = "DALY Components",
  "dw_chronic_severe" = "DALY Components",
  "dur_subac"         = "DALY Components",
  "dw_subac"          = "DALY Components",
  "dur_6m"            = "DALY Components",
  "dur_12m"           = "DALY Components",
  "dur_30m"           = "DALY Components",
  "fatal_nonhosp"     = "DALY Components"
)

group_map_ar <- setNames(
  rep("Regional Attack Rate", length(cols_ar)),
  cols_ar
)

# Merge regional attack-rate columns into group_map
group_map <- c(group_map, group_map_ar)

group_colors <- c(
  "Vaccine Risk"           = "#E63946",
  "Natural Disease Risk"   = "#457B9D",
  "Vaccine Effectiveness"  = "#2A9D8F",
  "Attack Rate / Exposure" = "#E9C46A",
  "Symptomatic / Severity" = "#F4A261",
  "Life Expectancy"        = "#A8DADC",
  "DALY Components"        = "#9B72CF",
  "Regional Attack Rate"   = "#6DBF9E"
)

# Human-readable strip labels; only these (+ regional ar_*) are included in figures
param_display_labels_core <- c(
  p_sae_vacc_u65    = "Vaccine SAE rate (18-64 y)",
  p_sae_vacc_65     = "Vaccine SAE rate (65+ y)",
  p_death_vacc_u65  = "Vaccine death rate (18-64 y)",
  p_death_vacc_65   = "Vaccine death rate (65+ y)",
  p_sae_nat_64      = "Hospitalisation rate among symptomatic (18-64 y)",
  p_sae_nat_65      = "Hospitalisation rate among symptomatic (65+ y)",
  p_death_nat_64    = "Death rate among symptomatic (18-64 y)",
  p_death_nat_65    = "Death rate among symptomatic (65+ y)",
  ve                = "Vaccine efficacy",
  epi_months        = "Epidemic season length (months)",
  trav_7d           = "Travel duration (7 d scenario, days)",
  trav_14d          = "Travel duration (14 d scenario, days)",
  trav_30d          = "Travel duration (30 d scenario, days)",
  trav_90d          = "Travel duration (90 d scenario, days)",
  symp_overall      = "Symptomatic fraction (pooled)",
  fatal_hosp        = "Case fatality, hospitalised acute",
  hosp              = "Hospitalisation given symptomatic infection",
  lt                = "Long-term outcome weight / proportion (lt)",
  le_lost_18_64     = "Life-years lost, ages 18-64",
  le_lost_65        = "Life-years lost, ages 65+",
  dw_chronic        = "Disability weight, chronic sequelae",
  dur_chronic       = "Duration, chronic sequelae (generic)",
  dw_hosp           = "Disability weight, hospitalised acute",
  dur_acute         = "Duration, hospitalised acute phase",
  dw_nonhosp        = "Disability weight, non-hospitalised acute",
  dur_nonhosp       = "Duration, non-hospitalised acute",
  acute             = "Proportion, acute-only health state",
  subac             = "Proportion, subacute health state",
  chr6m             = "Proportion, chronic up to 6 months",
  chr12m            = "Proportion, chronic up to 12 months",
  chr30m            = "Proportion, chronic up to 30 months",
  dur_subac         = "Duration, subacute phase",
  dw_subac          = "Disability weight, subacute",
  dur_6m            = "Duration, chronic up to 6-month",
  dur_12m           = "Duration, chronic up to 12-month",
  dur_30m           = "Duration, chronic up to 30-month",
  fatal_nonhosp     = "Case fatality, non-hospitalised acute"
)

param_display_ar <- stats::setNames(
  paste0("State attack-rate draw (scaled), ", names(region_key)),
  cols_ar
)
param_display_labels <- c(param_display_labels_core, param_display_ar)

# Plot titles by thematic group (clearer than "LHS ...")
group_plot_titles <- c(
  "Vaccine Risk"           = "Probabilistic inputs: vaccine-attributable risks",
  "Natural Disease Risk"   = "Probabilistic inputs: natural infection risks",
  "Vaccine Effectiveness"  = "Probabilistic inputs: vaccine efficacy",
  "Attack Rate / Exposure" = "Probabilistic inputs: attack rate, season length, travel duration",
  "Symptomatic / Severity" = "Probabilistic inputs: symptomatic fraction (pooled)",
  "Life Expectancy"        = "Probabilistic inputs: life-years lost by age group",
  "DALY Components"        = "Probabilistic inputs: disability weights, durations, and fatality",
  "Regional Attack Rate"   = "Probabilistic inputs: state-specific attack-rate draws"
)

# Parameters held fixed (excluded from plots)
zero_fixed_params <- c("p_death_vacc_u65")
# Exclude child / adolescent strata from figures (adult-focused outputs)
exclude_age_params <- c(
  "p_sae_nat_11", "p_sae_nat_17",
  "p_death_nat_11", "p_death_nat_17",
  "le_lost_1_11", "le_lost_12_17"
)
exclude_plot_params <- unique(c(zero_fixed_params, exclude_age_params))

# Proportion-scale attack rates: plot on % (×100); epi_months and trav_* stay in natural units
params_attack_rate_prop <- c("ar_small", "ar_med", "ar_large", cols_ar)

plot_param_core_ordered <- c(names(param_display_labels_core), cols_ar)
# lhs_sample colnames follow cols2 in setup.R (cols + regional ar_*); `cols` omits ar_*
plot_param_levels <- intersect(plot_param_core_ordered, names(lhs_sample))
plot_param_levels <- setdiff(plot_param_levels, exclude_plot_params)
param_label_levels <- vapply(plot_param_levels, function(p) {
  dplyr::coalesce(unname(param_display_labels[p]), p)
}, character(1))

param_label_df <- data.frame(
  param = plot_param_levels,
  param_label = factor(param_label_levels, levels = param_label_levels),
  stringsAsFactors = FALSE
)

# --- 3. Long format: filtered parameters + display labels for facets ---
df_long <- as.data.frame(lhs_sample) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  dplyr::filter(.data$param %in% plot_param_levels) %>%
  dplyr::mutate(param = factor(.data$param, levels = plot_param_levels)) %>%
  dplyr::left_join(param_label_df, by = "param") %>%
  dplyr::mutate(group = unname(group_map[as.character(.data$param)])) %>%
  dplyr::filter(!is.na(.data$group)) %>%
  dplyr::mutate(
    value = dplyr::if_else(
      as.character(.data$param) %in% params_attack_rate_prop,
      .data$value * 100,
      .data$value
    )
  )

# --- 4. Summary statistics by group and parameter ---
summary_tbl <- df_long %>%
  group_by(.data$group, .data$param, .data$param_label) %>%
  summarise(
    median = median(.data$value),
    p2.5   = quantile(.data$value, 0.025),
    p97.5  = quantile(.data$value, 0.975),
    mean   = mean(.data$value),
    .groups = "drop"
  )

print(summary_tbl, n = 100)

# --- 5. Histograms by group (display and save to PDF) ---
# (Run this script from the top, or source the whole file — mid-file runs skip data prep.)
preferred_group_order <- c(
  "Vaccine Risk",
  "Natural Disease Risk",
  "Vaccine Effectiveness",
  "Attack Rate / Exposure",
  "Symptomatic / Severity",
  "Life Expectancy",
  "DALY Components",
  "Regional Attack Rate"
)
groups <- preferred_group_order[preferred_group_order %in% unique(df_long$group)]

lhs_xlab_for_group <- function(grp) {
  if (identical(grp, "Regional Attack Rate")) {
    return("Attack rate (%)")
  }
  if (identical(grp, "Attack Rate / Exposure")) {
    return("Value (attack-rate scenarios: %; epidemic season: months; travel: days)")
  }
  "Simulated value"
}

subtitle_base <- paste0(
  "Latin hypercube sample, n = ", format(runs, big.mark = ","),
  "  |  Dashed: median; dotted: equal-tailed 2.5th and 97.5th percentiles"
)

for (grp in groups) {
  df_sub <- df_long %>% filter(.data$group == grp)
  if (grp == "Symptomatic / Severity") {
    df_sub <- df_sub %>% dplyr::filter(.data$param == "symp_overall")
  }

  subtitle_text <- if (grp == "Vaccine Risk") {
    paste0(
      subtitle_base,
      "  |  Vaccine death rate (18-64 y) fixed at 0 \n(no observed events; excluded from panels)"
    )
  } else if (grp == "Attack Rate / Exposure") {
    paste0(
      subtitle_base,
      "  |  Small/medium/large attack-rate scenarios use a percentage scale on the x-axis"
    )
  } else {
    subtitle_base
  }

  plot_title <- group_plot_titles[[grp]]
  if (is.null(plot_title) || length(plot_title) == 0 || is.na(plot_title) || !nzchar(plot_title)) {
    plot_title <- paste0("Probabilistic input distributions: ", grp)
  }

  vline_dat <- summary_tbl %>% dplyr::filter(.data$group == grp)
  if (grp == "Symptomatic / Severity") {
    vline_dat <- dplyr::filter(vline_dat, .data$param == "symp_overall")
  }

  p <- ggplot(df_sub, aes(x = .data$value)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 60,
      fill = group_colors[grp],
      color = "white",
      alpha = 0.85
    ) +
    geom_density(color = "grey30", linewidth = 0.5) +
    facet_wrap(
      ~param_label,
      scales = "free",
      ncol = 4,
      labeller = ggplot2::labeller(param_label = ggplot2::label_wrap_gen(width = 28))
    ) +
    labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = lhs_xlab_for_group(grp), y = "Density"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text       = element_text(face = "bold", size = 7.5, lineheight = 0.95),
      strip.text.x   = element_text(margin = margin(b = 3, t = 2)),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(color = "grey40", size = 8.5, lineheight = 1.1),
      plot.margin      = margin(10, 10, 14, 10),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 35, hjust = 1, size = 10)
    ) +
    geom_vline(
      data = vline_dat,
      aes(xintercept = .data$median),
      color = "grey20", linetype = "dashed", linewidth = 0.6
    ) +
    geom_vline(
      data = vline_dat,
      aes(xintercept = .data$p2.5),
      color = "grey20", linetype = "dotted", linewidth = 0.4
    ) +
    geom_vline(
      data = vline_dat,
      aes(xintercept = .data$p97.5),
      color = "grey20", linetype = "dotted", linewidth = 0.4
    )

  if (identical(grp, "Regional Attack Rate")) {
    p <- p + ggplot2::scale_x_continuous(labels = scales::label_number(suffix = "%"))
  }

  print(p)
}

output_dir <- "03_Manuscript/Submission/second_submission_lancet/suppl_figs/fig1_lhs_dist"
dir.create(output_dir, showWarnings = FALSE)

for (grp in groups) {
  df_sub      <- df_long %>% filter(.data$group == grp)
  if (grp == "Symptomatic / Severity") {
    df_sub <- df_sub %>% dplyr::filter(.data$param == "symp_overall")
  }
  summary_sub <- summary_tbl %>% filter(.data$group == grp)
  if (grp == "Symptomatic / Severity") {
    summary_sub <- summary_sub %>% dplyr::filter(.data$param == "symp_overall")
  }
  n_params    <- n_distinct(df_sub$param)

  n_cols <- min(4L, n_params)
  n_rows <- ceiling(n_params / n_cols)

  subtitle_text <- if (grp == "Vaccine Risk") {
    paste0(
      subtitle_base,
      "  |  Vaccine death rate (18-64 y) fixed at 0 \n(excluded from panels)"
    )
  } else if (grp == "Attack Rate / Exposure") {
    paste0(
      subtitle_base,
      "  |  Small/medium/large attack-rate scenarios use a percentage scale on the x-axis"
    )
  } else {
    subtitle_base
  }

  plot_title <- group_plot_titles[[grp]]
  if (is.null(plot_title) || length(plot_title) == 0 || is.na(plot_title) || !nzchar(plot_title)) {
    plot_title <- paste0("Probabilistic input distributions: ", grp)
  }

  p <- ggplot(df_sub, aes(x = .data$value)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins   = 60,
      fill   = group_colors[grp],
      color  = "white",
      alpha  = 0.85
    ) +
    geom_density(color = "grey30", linewidth = 0.5) +
    geom_vline(
      data = summary_sub,
      aes(xintercept = .data$median),
      color = "grey20", linetype = "dashed", linewidth = 0.6
    ) +
    geom_vline(
      data = summary_sub,
      aes(xintercept = .data$p2.5),
      color = "grey20", linetype = "dotted", linewidth = 0.4
    ) +
    geom_vline(
      data = summary_sub,
      aes(xintercept = .data$p97.5),
      color = "grey20", linetype = "dotted", linewidth = 0.4
    ) +
    facet_wrap(
      ~param_label,
      scales = "free",
      ncol = n_cols,
      labeller = ggplot2::labeller(param_label = ggplot2::label_wrap_gen(width = 28))
    ) +
    labs(
      title    = plot_title,
      subtitle = subtitle_text,
      x = lhs_xlab_for_group(grp), y = "Density"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text       = element_text(face = "bold", size = 8, lineheight = 0.95),
      strip.text.x   = element_text(margin = margin(b = 3, t = 2)),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(color = "grey40", size = 8.5, lineheight = 1.1),
      plot.margin      = margin(10, 10, 16, 10),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 35, hjust = 1, size = 10),
      axis.text.y      = element_text(size = 10)
    )

  if (identical(grp, "Regional Attack Rate")) {
    p <- p + ggplot2::scale_x_continuous(labels = scales::label_number(suffix = "%"))
  }

  fname <- file.path(
    output_dir,
    paste0(sprintf("%02d", which(groups == grp)), "_",
           gsub(" ", "_", gsub("/", "", grp)), ".pdf")
  )

  # Inches scale more reliably than cm for multi-facet PDFs; extra height avoids clipped strips/axes
  fig_w_in <- n_cols * 3.1
  fig_h_in <- n_rows * 2.65 + 1.35
  # Few panels (e.g. 03 Vaccine Effectiveness = single VE facet): bump size so title/strip/axes are not clipped
  if (n_params <= 2L) {
    fig_w_in <- max(fig_w_in, 6.2)
    fig_h_in <- max(fig_h_in, 5.2)
  }

  ggsave(
    fname, p,
    width = fig_w_in,
    height = fig_h_in,
    units = "in",
    dpi = 300,
    limitsize = FALSE,
    # capabilities() is in base, not grDevices::capabilities
    device = if (isTRUE(capabilities()[["cairo"]])) {
      grDevices::cairo_pdf
    } else {
      "pdf"
    }
  )

  cat("Saved:", fname, "\n")
}

cat("\nAll plots saved to:", output_dir, "\n")