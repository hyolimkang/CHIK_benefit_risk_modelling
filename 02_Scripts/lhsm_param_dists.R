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

# 기존 group_map에 합치기
group_map <- c(group_map, group_map_ar)

# group_colors에도 추가
group_colors <- c(
  group_colors,
  "Regional Attack Rate" = "#6DBF9E"
)
# ── 3. Long format 변환 ────────────────────────────────────────────────────
df_long <- as.data.frame(lhs_sample) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  mutate(group = group_map[param],
         param = factor(param, levels = cols))

# ── 4. 요약 통계 테이블 ─────────────────────────────────────────────────────
summary_tbl <- df_long %>%
  group_by(group, param) %>%
  summarise(
    median = median(value),
    p2.5   = quantile(value, 0.025),
    p97.5  = quantile(value, 0.975),
    mean   = mean(value),
    .groups = "drop"
  )

print(summary_tbl, n = 51)

# ── 5. 그룹별 히스토그램 저장 ─────────────────────────────────────────────
groups <- unique(df_long$group)
group_colors <- c(
  "Vaccine Risk"           = "#E63946",
  "Natural Disease Risk"   = "#457B9D",
  "Vaccine Effectiveness"  = "#2A9D8F",
  "Attack Rate / Exposure" = "#E9C46A",
  "Symptomatic / Severity" = "#F4A261",
  "Life Expectancy"        = "#A8DADC",
  "DALY Components"        = "#9B72CF"
)

for (grp in groups) {
  df_sub <- df_long %>% filter(group == grp)
  
  p <- ggplot(df_sub, aes(x = value)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 60,
      fill = group_colors[grp],
      color = "white",
      alpha = 0.85
    ) +
    geom_density(color = "grey30", linewidth = 0.5) +
    facet_wrap(~ param, scales = "free", ncol = 4) +
    labs(
      title    = paste0("LHS Parameter Distributions — ", grp),
      subtitle = paste0("n = ", format(runs, big.mark = ","), " samples  |  Median ± 95% CrI shown"),
      x = "Value", y = "Density"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text       = element_text(face = "bold", size = 8),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(color = "grey40", size = 9),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 30, hjust = 1, size = 7)
    ) +
    # median + 95% CrI 수직선
    geom_vline(
      data = summary_tbl %>% filter(group == grp),
      aes(xintercept = median),
      color = "grey20", linetype = "dashed", linewidth = 0.6
    ) +
    geom_vline(
      data = summary_tbl %>% filter(group == grp),
      aes(xintercept = p2.5),
      color = "grey20", linetype = "dotted", linewidth = 0.4
    ) +
    geom_vline(
      data = summary_tbl %>% filter(group == grp),
      aes(xintercept = p97.5),
      color = "grey20", linetype = "dotted", linewidth = 0.4
    )
  
  print(p) 
  
}