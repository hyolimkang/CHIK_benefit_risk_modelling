library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
risk_summary <- lhs_sample %>%
  summarise(
    sae_u65_med = median(p_sae_vacc_u65),
    sae_u65_lo  = quantile(p_sae_vacc_u65, 0.025),
    sae_u65_hi  = quantile(p_sae_vacc_u65, 0.975),
    sae_65_med  = median(p_sae_vacc_65),
    sae_65_lo   = quantile(p_sae_vacc_65, 0.025),
    sae_65_hi   = quantile(p_sae_vacc_65, 0.975),
    death_u65_med = median(p_death_vacc_u65),
    death_u65_lo  = quantile(p_death_vacc_u65, 0.025),
    death_u65_hi  = quantile(p_death_vacc_u65, 0.975),
    death_65_med  = median(p_death_vacc_65),
    death_65_lo   = quantile(p_death_vacc_65, 0.025),
    death_65_hi   = quantile(p_death_vacc_65, 0.975)
  )

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
df_both_ve <- combined_nnv_national_age_ixchiq %>%
  filter(VE %in% c("VE0", "VE98.9"), VC == "cov50") %>%
  mutate(
    AgeCat = ifelse(AgeGroup < 16, "18-64", "65+"),
    
    risk_sae_med    = ifelse(AgeGroup < 16, risk_summary$sae_u65_med,  risk_summary$sae_65_med),
    risk_sae_lo     = ifelse(AgeGroup < 16, risk_summary$sae_u65_lo,   risk_summary$sae_65_lo),
    risk_sae_hi     = ifelse(AgeGroup < 16, risk_summary$sae_u65_hi,   risk_summary$sae_65_hi),
    risk_death_med  = ifelse(AgeGroup < 16, risk_summary$death_u65_med, risk_summary$death_65_med),
    risk_death_lo   = ifelse(AgeGroup < 16, risk_summary$death_u65_lo,  risk_summary$death_65_lo),
    risk_death_hi   = ifelse(AgeGroup < 16, risk_summary$death_u65_hi,  risk_summary$death_65_hi),
    
    vax_sae_med   = tot_vacc * risk_sae_med,
    vax_sae_lo    = tot_vacc * risk_sae_lo,
    vax_sae_hi    = tot_vacc * risk_sae_hi,
    vax_death_med = tot_vacc * risk_death_med,
    vax_death_lo  = tot_vacc * risk_death_lo,
    vax_death_hi  = tot_vacc * risk_death_hi
  )

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
death_df <- df_both_ve %>%
  transmute(
    scenario, VE, VC, AgeCat,
    outcome = "Death",
    tot_vacc,
    
    x_med = vax_death_med,
    x_lo  = vax_death_lo,
    x_hi  = vax_death_hi,
    
    y_med = diff_fatal,
    y_lo  = diff_fatal_low,
    y_hi  = diff_fatal_hi
  )

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
hosp_df <- df_both_ve %>%
  transmute(
    scenario, VE, VC, AgeCat,
    outcome = "SAE",
    tot_vacc,
    
    # X축: 백신 SAE (excess)
    x_med = vax_sae_med,
    x_lo  = vax_sae_lo,
    x_hi  = vax_sae_hi,
    
    # Y축: 예방된 입원 (averted)
    y_med = diff_hosp,
    y_lo  = diff_hosp_lo,
    y_hi  = diff_hosp_hi
  )

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
summary_df <- bind_rows(death_df, hosp_df) %>%
  group_by(scenario, VE, VC, outcome, AgeCat) %>%
  summarise(
    x_med = sum(x_med, na.rm = TRUE),
    x_lo  = sum(x_lo,  na.rm = TRUE),
    x_hi  = sum(x_hi,  na.rm = TRUE),
    y_med = sum(y_med, na.rm = TRUE),
    y_lo  = sum(y_lo,  na.rm = TRUE),
    y_hi  = sum(y_hi,  na.rm = TRUE),
    tot_vacc_grp = sum(tot_vacc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_med = (x_med / tot_vacc_grp) * 1e4,
    x_lo  = (x_lo  / tot_vacc_grp) * 1e4,
    x_hi  = (x_hi  / tot_vacc_grp) * 1e4,
    y_med = (y_med / tot_vacc_grp) * 1e4,
    y_lo  = (y_lo  / tot_vacc_grp) * 1e4,
    y_hi  = (y_hi  / tot_vacc_grp) * 1e4,
    
    brr = y_med / x_med
  ) %>%
  filter(is.finite(x_med) & is.finite(y_med)) %>%
  dplyr::mutate(outcome = case_when(
    outcome == "SAE" ~ "Hospitalisation",
    TRUE ~ outcome  
  ))


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# 1. Data Preparation and Labeling
# ------------------------------------------------------------------------------

# Convert VE to factors with the requested descriptive labels
summary_df <- summary_df %>%
  mutate(VE_label = factor(VE, 
                           levels = c("VE0", "VE98.9"),
                           labels = c("Disease blocking only", "Disease and infection blocking")))

# Define scenario labels for facet headers
scenario_labels <- c(
  "Scenario_1" = "1–11 years",
  "Scenario_2" = "12–17 years",
  "Scenario_3" = "18–64 years",
  "Scenario_4" = "65+ years"
)

# ------------------------------------------------------------------------------
# 2. Global Limit Calculation with Buffer (To Prevent Clipping)
# ------------------------------------------------------------------------------

# Calculate global maximums and add a 15% buffer to prevent shapes from being cut off
x_max_buffered <- max(summary_df$x_hi, na.rm = TRUE) * 1.15
y_max_buffered <- max(summary_df$y_hi, na.rm = TRUE) * 1.15

# Set a small negative buffer for the origin to ensure points at 0 are fully visible
x_min_buffered <- -x_max_buffered * 0.01
y_min_buffered <- -y_max_buffered * 0.01


# ------------------------------------------------------------------------------
# 4. Final Visualization
# ------------------------------------------------------------------------------

create_br_plot <- function(data, target_outcome, log_min = -1, log_max = 3) {
  
  # 1. 해당 Outcome 데이터만 필터링
  plot_data <- data %>% filter(outcome == target_outcome)
  
  # 2. 패널별(AgeCat) 맞춤형 x, y 범위 계산 (Outcome specific)
  panel_limits <- plot_data %>%
    group_by(AgeCat) %>%
    summarise(
      x_max = max(x_hi, na.rm = TRUE) * 1.1,
      y_max = max(y_hi, na.rm = TRUE) * 1.1,
      .groups = "drop"
    )
  
  # 3. 패널별 배경 그리드 생성
  bg_grid_specific <- panel_limits %>%
    group_by(AgeCat) %>%
    reframe({
      x_seq <- seq(0, x_max, length.out = 100)
      y_seq <- seq(0, y_max, length.out = 100)
      grid <- expand.grid(x = x_seq, y = y_seq)
      grid$brr <- grid$y / (grid$x + 1e-9)
      grid$log10_brr <- pmax(pmin(log10(grid$brr), log_max), log_min)
      grid
    })
  
  # 4. ggplot 생성
  ggplot() +
    # 배경 (Outcome 전용 범위 적용)
    geom_raster(
      data = bg_grid_specific,
      aes(x = x, y = y, fill = log10_brr),
      alpha = 0.8, interpolate = TRUE
    ) +
    scale_fill_gradient2(
      name = "Benefit–risk ratio",
      low = "#ca0020", mid = "#f7f7f7", high = "#0571b0", midpoint = 0,
      limits = c(log_min, log_max),
      oob = scales::squish,
      labels = function(x) paste0("10^", x)
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.4) +
    
    # 에러바 및 포인트
    geom_errorbar(data = plot_data, aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome), width = 0) +
    geom_errorbarh(data = plot_data, aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome), height = 0) +
    geom_point(data = plot_data, aes(x = x_med, y = y_med, shape = VE_label, color = outcome), 
               size = 3, fill = "white", stroke = 1) +
    
    # 연령별 패널 분리 및 축 범위 개별화
    facet_wrap(~ AgeCat, scales = "free", ncol = 2) +
    
    # 색상 및 스타일 (Outcome별 고유 색상 유지)
    scale_color_manual(values = c("SAE" = "#1B7F1B", "Death" = "#B8860B", "DALY" = "#A23B72")) +
    scale_shape_manual(values = c("Disease blocking only" = 21, "Disease and infection blocking" = 24)) +
    
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    
    labs(
      title = paste("Benefit-risk Assessment:", target_outcome),
      x = "Excess adverse outcomes (per 10,000 vaccinated)",
      y = "Averted outcomes (per 10,000 vaccinated)"
    ) +
    theme_bw() +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
}

plot_sae   <- create_br_plot(summary_df, "SAE")
plot_death <- create_br_plot(summary_df, "Death")
plot_daly  <- create_br_plot(summary_df, "DALY")

plot_daly
plot_sae
plot_death

combined_plot <- plot_vaccine_mechanism_final / plot_brr_ori + 
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = '',
    tag_suffix = '.',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18, hjust = 0, vjust = 1)
    )
  )

ggsave("06_Results/brr_final_plot2.pdf", plot = combined_plot, width = 10, height = 11)
