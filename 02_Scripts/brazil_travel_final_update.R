setting_key <- c(
  "Ceará"             = "High",
  "Bahia"             = "Low",
  "Paraíba"           = "High",
  "Pernambuco"        = "Moderate",
  "Rio Grande do Norte" = "Low",
  "Piauí"             = "High",
  "Tocantins"         = "Moderate",
  "Alagoas"           = "High",
  "Minas Gerais"      = "Low",
  "Sergipe"           = "Low",
  "Goiás"             = "Low"
)

fn_br_space_benefit_by_setting <- function(psa_data,
                                           setting_key,
                                           na_rm = TRUE,
                                           days_levels = c("7d","14d","30d","90d")) {
  

  psa_data %>%
    mutate(
      days = factor(days, levels = days_levels),
      setting = unname(setting_key[state])
    ) %>%
    filter(!is.na(setting)) %>%
    pivot_longer(
      cols = c(
        excess_10k_sae, excess_10k_death,
        averted_10k_sae, averted_10k_death
      ),
      names_to = c(".value", "outcome"),
      names_pattern = "(excess_10k|averted_10k)_(sae|death)"
    ) %>%
    mutate(outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death")) %>%
    group_by(setting, age_group, days, outcome) %>%
    summarise(
      # x-axis: vaccine risk
      x_med = median(excess_10k, na.rm = na_rm),
      x_lo  = quantile(excess_10k, 0.025, na.rm = na_rm),
      x_hi  = quantile(excess_10k, 0.975, na.rm = na_rm),
      
      # y-axis: benefit averted
      y_med = median(averted_10k, na.rm = na_rm),
      y_lo  = quantile(averted_10k, 0.025, na.rm = na_rm),
      y_hi  = quantile(averted_10k, 0.975, na.rm = na_rm),
      
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(setting = factor(setting, levels = c("Low","Moderate","High")))
}
fn_br_space_daly_by_setting <- function(psa_data,
                                        setting_key,
                                        na_rm = TRUE,
                                        days_levels = c("7d", "14d", "30d", "90d")) {
  
  library(dplyr)
  
  psa_data %>%
    mutate(
      days = factor(days, levels = days_levels),
      outcome = "DALY",
      setting = unname(setting_key[state])
    ) %>%
    filter(!is.na(setting)) %>%
    group_by(setting, age_group, days, outcome) %>%
    summarise(
      # x-axis: vaccine harm (DALY from SAE)
      x_med = median(daly_sae, na.rm = na_rm),
      x_lo  = quantile(daly_sae, 0.025, na.rm = na_rm),
      x_hi  = quantile(daly_sae, 0.975, na.rm = na_rm),
      
      # y-axis: benefit (DALY averted from disease)
      y_med = median(daly_averted, na.rm = na_rm),
      y_lo  = quantile(daly_averted, 0.025, na.rm = na_rm),
      y_hi  = quantile(daly_averted, 0.975, na.rm = na_rm),
      
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(setting = factor(setting, levels = c("Low", "Moderate", "High")))
}

fn_panel_range <- function(br_representative_benefit,
                           group_var = "setting") {
  
  library(dplyr)
  
  g <- rlang::sym(group_var)
  
  panel_ranges <- br_representative_benefit %>%
    group_by(outcome, !!g, days) %>%
    summarise(
      x_min = min(c(0, x_lo), na.rm = TRUE),
      x_max = max(x_hi, na.rm = TRUE) * 1.05,
      y_min = min(c(0, y_lo), na.rm = TRUE),
      y_max = max(y_hi, na.rm = TRUE) * 1.05,
      .groups = "drop"
    )
  
  panel_ranges
}

fn_br_grid <- function(panel_ranges, group_var = "setting") {
  
  library(dplyr)
  library(rlang)
  
  g <- rlang::sym(group_var)
  
  bg_grid_optimized <- panel_ranges %>%
    rowwise() %>%
    do({
      panel_data <- .
      
      x_seq <- seq(panel_data$x_min, panel_data$x_max, length.out = 200)
      y_seq <- seq(panel_data$y_min, panel_data$y_max, length.out = 200)
      
      grid <- expand.grid(x = x_seq, y = y_seq)
      
      grid$brr <- with(grid, ifelse(x > 0, y / x, NA_real_))
      grid$log10_brr <- log10(grid$brr)
      
      grid$outcome <- panel_data$outcome
      grid$days    <- panel_data$days
      
      # group label(=setting or ar_category etc.)
      grid[[group_var]] <- panel_data[[group_var]]
      
      grid
    }) %>%
    ungroup()
  
  return(bg_grid_optimized)
}


fn_br_summ <- function(br_representative_benefit,
                       group_var = "setting",
                       age_var = NULL,
                       na_rm = TRUE){
  
  library(dplyr)
  library(rlang)
  
  g <- rlang::sym(group_var)
  
  if (is.null(age_var)) {
    if ("AgeCat" %in% names(br_representative_benefit)) {
      age_var <- "AgeCat"
    } else if ("age_group" %in% names(br_representative_benefit)) {
      age_var <- "age_group"
    } else {
      stop("Neither AgeCat nor age_group found in input data.")
    }
  }
  a <- rlang::sym(age_var)
  
  br_representative_benefit %>%
    group_by(outcome, !!g, days, !!a) %>%
    summarise(
      x_med = mean(x_med, na.rm = na_rm),
      y_med = mean(y_med, na.rm = na_rm),
      x_lo  = mean(x_lo,  na.rm = na_rm),
      x_hi  = mean(x_hi,  na.rm = na_rm),
      y_lo  = mean(y_lo,  na.rm = na_rm),
      y_hi  = mean(y_hi,  na.rm = na_rm),
      n = dplyr::n(),
      .groups = "drop"
    )
}
br_space_by_setting <- fn_br_space_benefit_by_setting(psa_df, setting_key)
br_space_daly_setting <- fn_br_space_daly_by_setting(psa_df, setting_key)

br_representative_benefit <- bind_rows(br_space_by_setting, br_space_daly_setting)

panel_ranges_benefit <- fn_panel_range(br_representative_benefit, group_var = "setting")

bg_grid_optimized <- fn_br_grid(panel_ranges_benefit, group_var = "setting")

br_summarized_setting <- fn_br_summ(br_representative_benefit, group_var = "setting")

log_min <- -2
log_max <- 2
log_range <- seq(log_min, log_max, by = 1)
brr_labels <- c("0.01", "0.1", "1", "10", "100")


plot_brr_outcome <- function(br_summarized, bg_grid_optimized,
                             target_outcome, title_text, color_val,
                             group_var = "setting",
                             shape_var = NULL,        
                             show_prop = TRUE,
                             eps_x = 1e-9,
                             keep_zero_axis = TRUE) {

  g <- rlang::sym(group_var)
  
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
  
  if (!(group_var %in% names(br_summarized))) {
    stop("Column `", group_var, "` is not found in br_summarized.")
  }
  if (!(group_var %in% names(bg_grid_optimized))) {
    stop("Column `", group_var, "` is not found in bg_grid_optimized.")
  }
  
  if (is.null(shape_var)) {
    if ("AgeCat" %in% names(br_summarized)) {
      shape_var <- "AgeCat"
    } else if ("age_group" %in% names(br_summarized)) {
      shape_var <- "age_group"
    } else {
      shape_var <- NULL
    }
  }
  s <- if (!is.null(shape_var)) rlang::sym(shape_var) else NULL
  
  plot_data <- br_summarized %>% filter(.data$outcome == target_outcome)
  plot_bg <- bg_grid_optimized %>%
    filter(.data$outcome == target_outcome, .data$x > eps_x)
  
  panel_prop <- plot_bg %>%
    mutate(is_fav = !is.na(log10_brr) & is.finite(log10_brr) & log10_brr > 0) %>%
    group_by(!!g, days) %>%
    summarise(prop_fav = mean(is_fav), .groups = "drop") %>%
    mutate(label = ifelse(prop_fav < 0.005, "BRR>1: <1%",
                          sprintf("BRR>1: %.0f%%", 100 * prop_fav)))
  
  panel_ranges <- plot_bg %>%
    group_by(!!g, days) %>%
    summarise(
      x_min = min(x, na.rm = TRUE),
      x_max = max(x, na.rm = TRUE),
      y_min = min(y, na.rm = TRUE),
      y_max = max(y, na.rm = TRUE),
      .groups = "drop"
    )
  
  panel_prop <- left_join(panel_prop, panel_ranges, by = c(group_var, "days"))
  
  facet_formula <- stats::as.formula(paste("~", group_var, "+ days"))
  
  p <- ggplot() +
    geom_raster(
      data = plot_bg,
      aes(x = x, y = y, fill = log10_brr),
      interpolate = FALSE, alpha = 0.85
    ) +
    scale_fill_gradient2(
      name = "Benefit–risk ratio",
      low = "#ca0020", mid = "#f7f7f7", high = "#0571b0",
      midpoint = 0, limits = c(log_min, log_max),
      breaks = log_range, labels = brr_labels,
      oob = scales::squish, na.value = "white"
    ) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", alpha = 0.4, linewidth = 0.9, colour = "grey35"
    ) +
    geom_errorbar(
      data = plot_data,
      aes(x = x_med, ymin = y_lo, ymax = y_hi),
      color = color_val, width = 0, linewidth = 0.5
    ) +
    geom_errorbarh(
      data = plot_data,
      aes(y = y_med, xmin = x_lo, xmax = x_hi),
      color = color_val, height = 0, linewidth = 0.5
    )
  
  if (!is.null(shape_var)) {
    p <- p +
      geom_point(
        data = plot_data,
        aes(x = x_med, y = y_med, shape = !!s),
        fill = "white", color = color_val, size = 1.5, stroke = 1
      )
  } else {
    p <- p +
      geom_point(
        data = plot_data,
        aes(x = x_med, y = y_med),
        fill = "white", color = color_val, size = 1.5, stroke = 1
      )
  }
  
  #if (show_prop) {
  #  p <- p +
  #    geom_label(
  #      data = panel_prop,
  #      aes(x = -Inf, y = Inf, label = label),
  #      hjust = 0, vjust = 1, size = 3,
  #      inherit.aes = FALSE,
  #      label.size = 0.25,
  #      fill = "white", alpha = 0.75, colour = "black"
  #    )
  #}
  
  p <- p +
    facet_wrap(facet_formula, scales = "free", ncol = 4) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = title_text,
      x = x_label,
      y = y_label
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "gray95"),
      legend.position = "right",
      panel.spacing = grid::unit(0.35, "lines")
    )
  
  if (!is.null(shape_var)) {
    p <- p + scale_shape_manual(
      values = c("1-11"=21, "12-17"=22, "18-64"=23, "65+"=24),
      name = "Age group"
    )
  }
  
  if (keep_zero_axis) {
    p <- p + coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), expand = FALSE)
  }
  
  return(p)
}

p_daly_mid  <- plot_brr_outcome(br_summarized_setting, bg_grid_optimized,
                                "DALY", "Benefit-Risk assessment: DALY", "#A23B72") + 
                theme(text = element_text(family = "Calibri"))+
                labs(tag = "B",
                     caption = "Note: Background colour indicates BRR = (DALYs averted by vaccination)/(DALYs attributable to vaccination) = y/x; dashed line indicates BRR = 1 (y = x).") +
                theme(
                plot.tag = element_text(face = "bold", size = 16),
                plot.tag.position = c(0, 1),
                plot.caption = element_text(hjust = 0, margin = margin(l = -8)),
                plot.caption.position = "plot",
                plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5)  
                 )


p_death_mid <- plot_brr_outcome(br_summarized_setting, bg_grid_optimized,
                                "Death", "Benefit-Risk assessment: Death", "#B8860B")+ 
                theme(text = element_text(family = "Calibri")) + 
                theme(text = element_text(family = "Calibri"))+
                labs(tag = "C",
                     caption = "Note: Background colour indicates BRR = (Deaths averted by vaccination)/(Deaths attributable to vaccination) = y/x; dashed line indicates BRR = 1 (y = x).") +
                theme(
                plot.tag = element_text(face = "bold", size = 16),
                plot.tag.position = c(0, 1),
                plot.caption = element_text(hjust = 0, margin = margin(l = -8)),
                plot.caption.position = "plot",
                plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 12)   
                )


p_sae_mid   <- plot_brr_outcome(br_summarized_setting, bg_grid_optimized,
                                "SAE",   "Benefit-Risk assessment: SAE",   "#1B7F1B")+ 
                theme(text = element_text(family = "Calibri")) +
                theme(text = element_text(family = "Calibri"))+
                labs(tag = "D",
                     caption = "Note: Background colour indicates BRR = (SAEs averted by vaccination)/(SAEs attributable to vaccination) = y/x; dashed line indicates BRR = 1 (y = x).") +
                theme(
                plot.tag = element_text(face = "bold", size = 16),
                plot.tag.position = c(0, 1),
                plot.caption = element_text(hjust = 0, margin = margin(l = -8)),
                plot.caption.position = "plot",
                plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5)  
                )


ggsave("06_Results/brr_travel_daly_mid.pdf", plot = p_daly_mid, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brr_travel_death_mid.pdf", plot = p_death_mid, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brr_travel_sae_mid.pdf", plot = p_sae_mid, width = 10, height = 8, device = cairo_pdf)


make_brr_long <- function(psa_df, setting_key) {
  psa_df %>%
    mutate(setting = unname(setting_key[state])) %>%
    filter(!is.na(setting)) %>%
    transmute(
      state, setting, days, age_group,
      brr_sae, brr_death, brr_daly
    ) %>%
    pivot_longer(
      cols = starts_with("brr_"),
      names_to = "outcome",
      values_to = "brr"
    ) %>%
    mutate(
      outcome = recode(outcome,
                       brr_sae = "SAE",
                       brr_death = "Death",
                       brr_daly = "DALY"),
      brr = as.numeric(brr)
    ) %>%
    filter(is.finite(brr), brr > 0)
}


make_brr_ceac <- function(brr_long,
                          thresholds = NULL,
                          group_vars = c("setting", "days", "age_group", "outcome"),
                          q_low = 0.001,          
                          q_high = 0.999,        
                          min_cap = 1e-3,        
                          max_cap = 1e3,          
                          step = 0.05             
) {
  library(dplyr)
  library(tidyr)
  
  if (is.null(thresholds)) {
    brr_vals <- brr_long$brr
    brr_vals <- brr_vals[is.finite(brr_vals) & brr_vals > 0]
    
    if (length(brr_vals) == 0) stop("No positive finite BRR values in brr_long.")
    
    lo <- as.numeric(quantile(brr_vals, q_low, na.rm = TRUE))
    hi <- as.numeric(quantile(brr_vals, q_high, na.rm = TRUE))
    
    lo <- max(lo, min_cap)
    hi <- min(hi, max_cap)
    
    lo_exp <- floor(log10(lo))
    hi_exp <- ceiling(log10(hi))
    
    thresholds <- 10^seq(lo_exp, hi_exp, by = step)
  }
  
  brr_long %>%
    tidyr::crossing(threshold = thresholds) %>%
    group_by(across(all_of(group_vars)), threshold) %>%
    summarise(
      p_accept = mean(brr > threshold),
      n = dplyr::n(),
      .groups = "drop"
    )
}

plot_brr_ceac <- function(ceac_df,
                          color_by = "age_group",
                          target_outcome = "DALY",
                          days_levels = c("7d","14d","30d","90d"),
                          setting_levels = c("Low","Moderate","High")) {
  

  df <- ceac_df %>%
    filter(outcome == target_outcome) %>%
    mutate(
      days = factor(days, levels = days_levels),
      setting = factor(setting, levels = setting_levels)
    )
  
  legend_title <- if (color_by == "age_group") "Age group" else color_by
  
  ggplot(df, aes(x = threshold, y = p_accept, colour = .data[[color_by]])) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.6) +
    scale_x_log10() +
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    facet_grid(setting ~ days) +   
    labs(
      x = "Benefit-risk ratio (BRR) threshold (t)",
      y = "Probability (BRR > t)",
      title = paste0("BRR acceptability curve: ", target_outcome),
      colour = legend_title
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
}
brr_long <- make_brr_long(psa_df, setting_key)

brr_long <- make_brr_long(psa_df, setting_key) %>%
  dplyr::mutate(days = factor(days, levels = c("7d","14d","30d","90d")))

ceac_df <- make_brr_ceac(
  brr_long
)

pr_ge1_tbl <- ceac_df %>%
  filter(threshold == 1) %>%
  transmute(
    setting, days, age_group, outcome,
    pr_brr_ge_1 = p_accept,
    n
  )

pr_gt1_wide <- ceac_df %>%
  mutate(
    setting   = factor(setting, levels = c("Low","Moderate","High")),
    days      = factor(days, levels = c("7d","14d","30d","90d")),
    age_group = ifelse(age_group == "65", "65+", age_group)
  ) %>%
  filter(threshold == 1) %>%
  transmute(outcome, setting, age_group, days,
            pr_brr_gt_1 = p_accept) %>%
  mutate(
    pr_fmt = sprintf("%.1f%%", 100 * pr_brr_gt_1)   
  ) %>%
  dplyr::select(-pr_brr_gt_1) %>%
  pivot_wider(
    names_from  = days,
    values_from = pr_fmt,
    names_glue  = "{days}_Pr(BRR>1)"
  )

p_daly_ceac <- plot_brr_ceac(ceac_df, target_outcome = "DALY") +
  theme(text = element_text(family = "Calibri"))+
  labs(tag = "E") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  )

p_death_ceac <- plot_brr_ceac(ceac_df, target_outcome = "Death")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "F") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  )

p_sae_ceac <- plot_brr_ceac(ceac_df, target_outcome = "SAE")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "G") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  )

ggsave("06_Results/brrac_daly_travel.pdf", plot = p_daly_ceac, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brrac_death_travel.pdf", plot = p_death_ceac, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brrac_sae_travel.pdf", plot = p_sae_ceac, width = 10, height = 8, device = cairo_pdf)


# summary of ceac_df
ceac_t1 <- ceac_df %>%
  filter(abs(log10(threshold)) < 1e-12) %>%   # == threshold=1
  mutate(
    days = factor(days, levels = c("7d","14d","30d","90d")),
    setting = factor(setting, levels = c("Low","Moderate","High"))
  ) %>%
  dplyr::select(outcome, setting, age_group, days, p_accept, n, threshold)

ceac_t1_wide <- ceac_t1 %>%
  mutate(p_fmt = sprintf("%.0f%%", 100 * p_accept)) %>%
  dplyr::select(outcome, setting, age_group, days, p_fmt, threshold) %>%
  pivot_wider(names_from = days, values_from = p_fmt) %>%
  arrange(outcome, setting, age_group)

write_xlsx(ceac_t1_wide, "06_Results/ceac_travel.xlsx")
write_xlsx(ceac_df, "06_Results/ceac_travel_all.xlsx")


## table

ar_summary_all <- psa_df %>%
  mutate(
    setting = unname(setting_key[state]),
    days = factor(days, levels = c("7d","14d","30d","90d")),
    age_group = ifelse(age_group == "65", "65+", age_group)
  ) %>%
  filter(!is.na(setting)) %>%
  pivot_longer(
    cols = c(brr_sae, brr_death, brr_daly),
    names_to = "outcome",
    values_to = "brr"
  ) %>%
  mutate(
    outcome = recode(outcome,
                     brr_sae   = "SAE",
                     brr_death = "Death",
                     brr_daly  = "DALY"),
    setting = factor(setting, levels = c("Low","Moderate","High"))
  ) %>%
  group_by(outcome, setting, age_group, days) %>%
  summarise(
    brr_med = median(brr, na.rm = TRUE),
    brr_lo  = quantile(brr, 0.025, na.rm = TRUE),
    brr_hi  = quantile(brr, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

ar_table_wide <- ar_summary_all %>%
  mutate(
    brr_formatted = sprintf("%.2f [%.2f–%.2f]", brr_med, brr_lo, brr_hi)
  ) %>%
  dplyr::select(outcome, setting, age_group, days, brr_formatted) %>%
  pivot_wider(names_from = days, values_from = brr_formatted) %>%
  arrange(outcome, setting, age_group)

idx_outcome <- table(ar_table_wide$outcome)

ar_table_wide2 <- ar_table_wide %>%
  left_join(pr_gt1_wide, by = c("outcome","setting","age_group"))

ar_table_wide2 <- ar_table_wide2 %>%
  relocate(`7d_Pr(BRR>1)`,  .after = `7d`)  %>%
  relocate(`14d_Pr(BRR>1)`, .after = `14d`) %>%
  relocate(`30d_Pr(BRR>1)`, .after = `30d`) %>%
  relocate(`90d_Pr(BRR>1)`, .after = `90d`)

kable(
  ar_table_wide,
  format = "html",
  col.names = c("Outcome", "Setting", "Age Group",
                "7 days", "14 days", "30 days", "90 days"),
  caption = "Benefit–Risk Ratio (BRR) by Outcome, Setting, Age Group, and Travel Duration",
  align = c("l", "c", "c", "r", "r", "r", "r")
) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    font_size = 12
  ) %>%
  pack_rows(index = idx_outcome) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2, bold = TRUE)