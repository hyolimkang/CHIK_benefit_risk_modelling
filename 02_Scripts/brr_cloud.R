library(dplyr)
library(tidyr)

make_cloud_long <- function(psa_df, setting_key,
                            days_levels = c("7d","14d","30d","90d")) {
  
  # Death cloud
  df_death <- psa_df %>%
    mutate(
      setting = unname(setting_key[state]),
      days = factor(days, levels = days_levels),
      outcome = "Death",
      caused_10k    = excess_10k_death,
      prevented_10k = averted_10k_death
    ) %>%
    dplyr::select(draw, state, setting, days, age_group, outcome, caused_10k, prevented_10k)
  
  # SAE cloud (optional but handy)
  df_sae <- psa_df %>%
    mutate(
      setting = unname(setting_key[state]),
      days = factor(days, levels = days_levels),
      outcome = "SAE",
      caused_10k    = excess_10k_sae,
      prevented_10k = averted_10k_sae
    ) %>%
    dplyr::select(draw, state, setting, days, age_group, outcome, caused_10k, prevented_10k)
  
  # DALY cloud (not per-10k; but we keep the same column names)
  df_daly <- psa_df %>%
    mutate(
      setting = unname(setting_key[state]),
      days = factor(days, levels = days_levels),
      outcome = "DALY",
      caused_10k    = daly_sae,
      prevented_10k = daly_averted
    ) %>%
    dplyr::select(draw, state, setting, days, age_group, outcome, caused_10k, prevented_10k)
  
  bind_rows(df_death, df_sae, df_daly) %>%
    filter(!is.na(setting)) %>%
    mutate(setting = factor(setting, levels = c("Low","Moderate","High")))
}


get_cloud_limits <- function(cloud_df, target_outcome, pad = 0.05) {
  df <- cloud_df %>% filter(outcome == target_outcome)
  
  x_max <- max(df$caused_10k, na.rm = TRUE)
  y_max <- max(df$prevented_10k, na.rm = TRUE)
  
  # 0부터 시작, 약간 여유(pad) 줌
  list(
    x_lim = c(0, x_max * (1 + pad)),
    y_lim = c(0, y_max * (1 + pad))
  )
}


library(ggplot2)

library(dplyr)
library(ggplot2)
library(scales)

plot_prob_cloud <- function(cloud_df, target_outcome,
                            days_levels = c("7d", "14d", "30d", "90d"),
                            add_summary = TRUE,
                            alpha_cloud = 0.03,
                            point_size = 0.9,
                            transform = c("pseudo_log", "none", "sqrt", "log10"),
                            eps = 1e-3,
                            sigma_pseudo = 0.1) {

  transform <- match.arg(transform)
  
  # Filter data and set factor levels
  df <- cloud_df %>%
    filter(outcome == target_outcome) %>%
    mutate(days = factor(days, levels = days_levels))
  
  # Calculate summary statistics
  if (add_summary) {
    summ <- df %>%
      group_by(setting, days, age_group) %>%
      summarise(
        x_med = median(caused_10k, na.rm = TRUE),
        x_lo  = quantile(caused_10k, 0.025, na.rm = TRUE),
        x_hi  = quantile(caused_10k, 0.975, na.rm = TRUE),
        y_med = median(prevented_10k, na.rm = TRUE),
        y_lo  = quantile(prevented_10k, 0.025, na.rm = TRUE),
        y_hi  = quantile(prevented_10k, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # Define axis labels
  label_suffix <- ifelse(target_outcome == "DALY", "", " (per 10,000 doses)")
  x_label <- paste0(target_outcome, "s caused by vaccination", label_suffix)
  y_label <- paste0(target_outcome, "s prevented by vaccination", label_suffix)
  
  # Handle log10 clipping for zeros
  if (transform == "log10") {
    df <- df %>% mutate(x_plot = pmax(caused_10k, eps), y_plot = pmax(prevented_10k, eps))
    if (add_summary) {
      summ <- summ %>% mutate(across(c(x_med, x_lo, x_hi), ~pmax(.x, eps)),
                              across(c(y_med, y_lo, y_hi), ~pmax(.x, eps)))
    }
  } else {
    df <- df %>% mutate(x_plot = caused_10k, y_plot = prevented_10k)
  }
  
  # Initialize plot
  p <- ggplot(df, aes(x = x_plot, y = y_plot, colour = age_group)) +
    geom_point(alpha = alpha_cloud, size = point_size) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    facet_grid(setting ~ days, scales = "free_y") +
    labs(title = paste0("Probabilistic cloud: ", target_outcome),
         x = x_label, y = y_label, colour = "Age group", shape = "Age group") +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  # Apply scales (pseudo_log is recommended for data with zeros)
  if (transform == "pseudo_log") {
    p <- p + scale_x_continuous(trans = pseudo_log_trans(sigma = sigma_pseudo)) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = sigma_pseudo))
  } else if (transform == "sqrt") {
    p <- p + scale_x_sqrt() + scale_y_sqrt()
  } else if (transform == "log10") {
    p <- p + scale_x_log10() + scale_y_log10()
  } else {
    lim <- get_cloud_limits(cloud_df, target_outcome)
    p <- p + coord_cartesian(xlim = lim$x_lim, expand = FALSE)
  }
  
  # Overlay summary with white-filled points and optimized error bars
  if (add_summary) {
    p <- p +
      geom_errorbar(data = summ, aes(x = x_med, ymin = y_lo, ymax = y_hi, colour = age_group),
                    width = 0, linewidth = 0.5, inherit.aes = FALSE) +
      geom_errorbarh(data = summ, aes(y = y_med, xmin = x_lo, xmax = x_hi, colour = age_group),
                     height = 0, linewidth = 0.5, inherit.aes = FALSE) +
      geom_point(data = summ, aes(x = x_med, y = y_med, colour = age_group, shape = age_group),
                 size = 2.8, fill = "white", stroke = 0.8, inherit.aes = FALSE) +
      scale_shape_manual(values = c(21, 24, 22, 23, 25))
  }
  
  return(p)
}


cloud_long <- make_cloud_long(psa_df, setting_key)

cloud_df_balanced <- cloud_long %>%
  group_by(setting, days, age_group, outcome) %>%
  slice_sample(n = 1000) %>% 
  ungroup()

p_cloud_death_log <- plot_prob_cloud(cloud_long, "Death", transform = "log10", eps = 1e-3) +
  theme(text = element_text(family = "Calibri"))

p_cloud_daly_log <- plot_prob_cloud(cloud_long, "DALY", transform = "log10", eps = 1e-3) +
  theme(text = element_text(family = "Calibri"))

p_cloud_sae_log <- plot_prob_cloud(cloud_long, "SAE", transform = "log10", eps = 1e-3) +
  theme(text = element_text(family = "Calibri"))

apc_long <- make_averted_per_caused_long(cloud_long)

apc_acc_dense <- make_apc_acceptability_dense(
  apc_long,
  t_min = 0.01, t_max = 100, n_grid = 300
)

p_apc_death <- plot_apc_acceptability(apc_acc_dense, target_outcome = "Death") +
  theme(text = element_text(family = "Calibri"))

######

library(dplyr)
library(ggplot2)

prep_cloud_pool_setting <- function(draw_level_xy_true,
                                    coverage_keep = "cov50",
                                    outcome_keep = c("SAE","Death","DALY")) {
  draw_level_xy_true %>%
    filter(Coverage == coverage_keep) %>%
    filter(outcome %in% outcome_keep) %>%
    mutate(
      # Scenario가 연령코드(1~4)라면 AgeCat으로 통일
      AgeCat = case_when(
        as.character(Scenario) == "1" ~ "1-11",
        as.character(Scenario) == "2" ~ "12-17",
        as.character(Scenario) == "3" ~ "18-64",
        as.character(Scenario) == "4" ~ "65+",
        TRUE ~ as.character(AgeCat)
      ),
      AgeCat = factor(AgeCat, levels = c("1-11","12-17","18-64","65+")),
      setting = factor(setting, levels = c("Low","Moderate","High")),
      VE_show = coalesce(VE_label, as.character(VE)),
      
      x = ifelse(outcome == "DALY", x_daly_10k, x_10k),
      y = y_10k
    ) %>%
    filter(is.finite(x), is.finite(y), x >= 0, y >= 0) %>%
    filter(!is.na(setting), !is.na(AgeCat), !is.na(VE_show))
}

plot_cloud_pool_auto_scale <- function(df,
                                       target_outcome = "SAE",
                                       # log10일 때만 eps 필요
                                       eps_log = 1e-2,
                                       # pseudo-log일 때 0 근처 완화 정도
                                       sigma_pseudo = 1e-3,
                                       alpha_cloud = 0.06,
                                       point_size = 0.7,
                                       add_summary = TRUE,
                                       drop_refline_for_SAE = FALSE) {
  
  library(dplyr)
  library(ggplot2)
  library(scales)
  
  dd <- df %>% filter(outcome == target_outcome)
  if (nrow(dd) == 0) stop("No rows left after filtering outcome/coverage.")
  
  # ====== outcome별 스케일 선택 ======
  mode <- if (target_outcome == "Death") "pseudo" else "log10"
  
  # plot 좌표 만들기
  if (mode == "log10") {
    dd <- dd %>% mutate(
      x_plot = pmax(x, eps_log),
      y_plot = pmax(y, eps_log)
    )
  } else { # pseudo
    dd <- dd %>% mutate(
      x_plot = x,
      y_plot = y
    )
  }
  
  # ====== 요약도 plot-space 기준으로 (중요) ======
  if (add_summary) {
    summ <- dd %>%
      group_by(setting, AgeCat, VE_show) %>%
      summarise(
        x_med = median(x_plot, na.rm = TRUE),
        x_lo  = quantile(x_plot, 0.025, na.rm = TRUE),
        x_hi  = quantile(x_plot, 0.975, na.rm = TRUE),
        y_med = median(y_plot, na.rm = TRUE),
        y_lo  = quantile(y_plot, 0.025, na.rm = TRUE),
        y_hi  = quantile(y_plot, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  x_label <- dplyr::case_when(
    target_outcome == "DALY"  ~ "DALYs attributable to vaccination (per 10,000 vaccinated)",
    target_outcome == "SAE"   ~ "SAEs attributable to vaccination (per 10,000 vaccinated)",
    target_outcome == "Death" ~ "Deaths attributable to vaccination (per 10,000 vaccinated)",
    TRUE ~ "Attributable to vaccination (per 10,000 vaccinated)"
  )
  y_label <- dplyr::case_when(
    target_outcome == "DALY"  ~ "DALYs averted by vaccination (per 10,000 vaccinated)",
    target_outcome == "SAE"   ~ "SAEs averted by vaccination (per 10,000 vaccinated)",
    target_outcome == "Death" ~ "Deaths averted by vaccination (per 10,000 vaccinated)",
    TRUE ~ "Averted by vaccination (per 10,000 vaccinated)"
  )
  
  p <- ggplot(dd, aes(x = x_plot, y = y_plot, colour = VE_show)) +
    geom_point(alpha = alpha_cloud, size = point_size) +
    facet_grid(setting ~ AgeCat) +
    labs(
      title = paste0("Probabilistic cloud: ", target_outcome,
                     " (", ifelse(mode == "log10","log10","pseudo-log"), " scale)"),
      x = x_label,
      y = y_label,
      colour = "Vaccine protection mechanism",
      shape  = "Vaccine protection mechanism"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  # 기준선(y=x): SAE에서 헷갈리면 옵션으로 제거
  if (!(drop_refline_for_SAE && target_outcome == "SAE")) {
    p <- p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40")
  }
  
  # 스케일 적용
  if (mode == "log10") {
    p <- p + scale_x_log10() + scale_y_log10()
  } else {
    p <- p +
      scale_x_continuous(trans = scales::pseudo_log_trans(sigma = sigma_pseudo)) +
      scale_y_continuous(trans = scales::pseudo_log_trans(sigma = sigma_pseudo))
  }
  
  # 요약 overlay: errorbar는 색만, point는 색+shape
  if (add_summary) {
    p <- p +
      geom_errorbar(
        data = summ,
        aes(x = x_med, ymin = y_lo, ymax = y_hi, colour = VE_show),
        width = 0, linewidth = 0.5, inherit.aes = FALSE
      ) +
      geom_errorbarh(
        data = summ,
        aes(y = y_med, xmin = x_lo, xmax = x_hi, colour = VE_show),
        height = 0, linewidth = 0.5, inherit.aes = FALSE
      ) +
      geom_point(
        data = summ,
        aes(x = x_med, y = y_med, colour = VE_show, shape = VE_show),
        size = 2.6, inherit.aes = FALSE, fill = "white", stroke = 0.9
      ) +
      scale_shape_manual(
        values = c(
          "Disease blocking only" = 21,
          "Disease and infection blocking" = 24
        )
      )
  }
  
  p
}
cloud_pool <- prep_cloud_pool_setting(draw_level_xy_true, coverage_keep = "cov50")

p_sae <- plot_cloud_pool_auto_scale(cloud_pool, "SAE", eps_log = 1e-2) +
  theme(text = element_text(family = "Calibri"))

p_daly <- plot_cloud_pool_auto_scale(cloud_pool, "DALY", eps_log = 1e-2) +
  theme(text = element_text(family = "Calibri"))

p_death <- plot_cloud_pool_auto_scale(cloud_pool, "Death", sigma_pseudo = 1e-2) +
  theme(text = element_text(family = "Calibri"))


p_death_sqrt <- plot_cloud_pool_death_linear(cloud_pool) +
  scale_x_continuous(trans = "sqrt") +
  scale_y_continuous(trans = "sqrt") +
  theme(text = element_text(family = "Calibri"))


