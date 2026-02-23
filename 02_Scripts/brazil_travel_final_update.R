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
    # 매핑 실패 방지
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
  
  # group_var가 "setting" 또는 "ar_category" 등일 때 모두 동작하게
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
  
  # AgeCat 있으면 그걸 쓰고, 없으면 age_group 사용
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
                             shape_var = NULL,        # <- "AgeCat" 또는 "age_group" 자동 선택
                             show_prop = TRUE,
                             eps_x = 1e-9,
                             keep_zero_axis = TRUE) {
  
  library(dplyr)
  library(ggplot2)
  library(rlang)
  
  g <- rlang::sym(group_var)
  
  if (!(group_var %in% names(br_summarized))) {
    stop("Column `", group_var, "` is not found in br_summarized.")
  }
  if (!(group_var %in% names(bg_grid_optimized))) {
    stop("Column `", group_var, "` is not found in bg_grid_optimized.")
  }
  
  # shape 변수 자동 선택
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
      color = color_val, width = 0, linewidth = 0.6
    ) +
    geom_errorbarh(
      data = plot_data,
      aes(y = y_med, xmin = x_lo, xmax = x_hi),
      color = color_val, height = 0, linewidth = 0.6
    )
  
  # ✅ point layer: shape_var가 있으면 shape mapped, 없으면 그냥 점
  if (!is.null(shape_var)) {
    p <- p +
      geom_point(
        data = plot_data,
        aes(x = x_med, y = y_med, shape = !!s),
        fill = "white", color = color_val, size = 3, stroke = 1
      )
  } else {
    p <- p +
      geom_point(
        data = plot_data,
        aes(x = x_med, y = y_med),
        fill = "white", color = color_val, size = 3, stroke = 1
      )
  }
  
  if (show_prop) {
    p <- p +
      geom_label(
        data = panel_prop,
        aes(x = -Inf, y = Inf, label = label),
        hjust = 0, vjust = 1, size = 3,
        inherit.aes = FALSE,
        label.size = 0.25,
        fill = "white", alpha = 0.75, colour = "black"
      )
  }
  
  p <- p +
    facet_wrap(facet_formula, scales = "free", ncol = 4) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = title_text,
      x = "Vaccine related excess outcome (per 10,000 vaccinated individuals)",
      y = "Outcomes averted by vaccination (per 10,000 vaccinated individuals)"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "gray95"),
      legend.position = "right",
      panel.spacing = grid::unit(0.35, "lines")
    )
  
  if (!is.null(shape_var)) {
    # age_group를 쓸 때도 기존 모양 매핑 유지
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
                                "DALY", "Benefit-Risk Assessment: DALY", "#A23B72")


p_death_mid <- plot_brr_outcome(br_summarized_setting, bg_grid_optimized,
                                "Death", "Benefit-Risk Assessment: Death", "#B8860B")

p_sae_mid   <- plot_brr_outcome(br_summarized_setting, bg_grid_optimized,
                                "SAE",   "Benefit-Risk Assessment: SAE",   "#1B7F1B")


library(dplyr)
library(tidyr)

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
                          thresholds = 10^seq(-2, 2, by = 0.05),  # 0.01 ~ 100
                          group_vars = c("setting", "days", "age_group", "outcome")) {
  
  library(dplyr)
  library(tidyr)
  
  # threshold grid를 붙여서 확률 계산
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
  
  library(dplyr)
  library(ggplot2)
  
  df <- ceac_df %>%
    filter(outcome == target_outcome) %>%
    mutate(
      days = factor(days, levels = days_levels),
      setting = factor(setting, levels = setting_levels)
    )
  
  ggplot(df, aes(x = threshold, y = p_accept, colour = .data[[color_by]])) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.6) +
    scale_x_log10() +
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    facet_grid(setting ~ days) +   # ✅ rows: Low→Moderate→High, cols: 7d→14d→30d→90d
    labs(
      x = "BRR threshold (t) [log scale]",
      y = "Probability(BRR > t)",
      title = paste0("BRR acceptability curve: ", target_outcome),
      colour = color_by
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
}
brr_long <- make_brr_long(psa_df, setting_key)

brr_long <- make_brr_long(psa_df, setting_key) %>%
  dplyr::mutate(days = factor(days, levels = c("7d","14d","30d","90d")))

ceac_df <- make_brr_ceac(
  brr_long,
  thresholds = 10^seq(-1, 1, by = 0.02)  
)

p_daly_ceac <- plot_brr_ceac(ceac_df, target_outcome = "DALY")
p_death_ceac <- plot_brr_ceac(ceac_df, target_outcome = "Death")
p_sae_ceac <- plot_brr_ceac(ceac_df, target_outcome = "SAE")



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