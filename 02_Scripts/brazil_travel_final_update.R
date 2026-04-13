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


make_brr_ceac_dt <- function(brr_long, 
                             thresholds = NULL, 
                             group_vars = c("setting", "days", "age_group", "outcome"),
                             q_low = 0.001, q_high = 0.999, 
                             min_cap = 1e-3, max_cap = 1e3, step = 0.05) {
  
  dt <- as.data.table(brr_long)
  
  if (is.null(thresholds)) {
    brr_vals <- dt[is.finite(brr) & brr > 0, brr]
    if (length(brr_vals) == 0) stop("No positive finite BRR values.")
    
    lo <- max(quantile(brr_vals, q_low), min_cap)
    hi <- min(quantile(brr_vals, q_high), max_cap)
    
    thresholds <- 10^seq(floor(log10(lo)), ceiling(log10(hi)), by = step)
  }
  
  thresh_dt <- data.table(threshold = thresholds)
  

  ceac_results <- dt[, {
    v <- .SD$brr
    lapply(thresholds, function(t) mean(v > t))
  }, by = group_vars]
  
  setnames(ceac_results, old = paste0("V", seq_along(thresholds)), new = as.character(thresholds))
  
  ceac_long <- melt(ceac_results, 
                    id.vars = group_vars, 
                    variable.name = "threshold", 
                    value.name = "p_accept")
  
  ceac_long[, threshold := as.numeric(as.character(threshold))]
  
  return(ceac_long)
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
      x = "Benefit-risk ratio (BRR)",
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

setDT(brr_long)

ceac_df <- make_brr_ceac_dt(brr_long)

pr_ge1_tbl <- ceac_df[threshold == 1, .(
  setting, 
  days, 
  age_group, 
  outcome, 
  pr_brr_ge_1 = p_accept,
  n = .N  
)]

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
  # Define pure averted and pure caused first
  mutate(
    # Benefit
    av_sae = averted_10k_sae, av_death = averted_10k_death, av_daly = daly_averted,
    # Risk (Pure vaccine only)
    ca_sae = excess_10k_sae, ca_death = excess_10k_death, ca_daly = daly_sae
  ) %>%
  pivot_longer(
    cols = c(contains("av_"), contains("ca_")),
    names_to = c(".value", "outcome"),
    names_sep = "_"
  ) %>%
  mutate(
    outcome = toupper(outcome),
    ca  = ifelse(ca == 0, NA_real_, ca),
    brr = av / ca
  ) %>%
  group_by(outcome, setting, age_group, days) %>%
  summarise(
    brr_med = median(brr, na.rm = TRUE),
    brr_lo  = quantile(brr, 0.025, na.rm = TRUE),
    brr_hi  = quantile(brr, 0.975, na.rm = TRUE),
    av_med = median(av, na.rm = TRUE),
    av_lo  = quantile(av, 0.025, na.rm = TRUE),
    av_hi  = quantile(av, 0.975, na.rm = TRUE),
    ca_med = median(ca, na.rm = TRUE),
    ca_lo  = quantile(ca, 0.025, na.rm = TRUE),
    ca_hi  = quantile(ca, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# ── 1. Format strings ──────────────────────────────────────────────────────────
ar_long_fmt <- ar_summary_all %>%
  mutate(
    ca_med = ifelse(ca_med < 1e-10, NA_real_, ca_med),
    ca_lo  = ifelse(ca_lo  < 1e-10, NA_real_, ca_lo),
    ca_hi  = ifelse(ca_hi  < 1e-10, NA_real_, ca_hi),
    brr_med = ifelse(is.infinite(brr_med) | brr_med > 1e10, NA_real_, brr_med),
    brr_lo  = ifelse(is.infinite(brr_lo)  | brr_lo  > 1e10, NA_real_, brr_lo),
    brr_hi  = ifelse(is.infinite(brr_hi)  | brr_hi  > 1e10, NA_real_, brr_hi),
    Benefit = sprintf("%.2f\n(%.2f–%.2f)", av_med, av_lo, av_hi),
    Risk    = sprintf("%.2f\n(%.2f–%.2f)", ca_med, ca_lo, ca_hi),
    BRR     = sprintf("%.2f\n(%.2f–%.2f)", brr_med, brr_lo, brr_hi)
  ) %>%
  dplyr::select(outcome, setting, age_group, days, Benefit, Risk, BRR)

# ── 2. Pivot wide by days ──────────────────────────────────────────────────────
ar_wide_ft <- ar_long_fmt %>%
  pivot_wider(
    names_from  = days,
    values_from = c(Benefit, Risk, BRR),
    names_glue  = "{days}_{.value}"
  ) %>%
  rename(
    Outcome     = outcome,
    Setting     = setting,
    `Age group` = age_group
  ) %>%
  mutate(
    Setting = factor(Setting, levels = c("High", "Moderate", "Low"))
  ) %>%
  arrange(Outcome, Setting, `Age group`)

# ── 3. Reorder columns ────────────────────────────────────────────────────────
day_cols <- paste0(rep(c("7d","14d","30d","90d"), each = 3),
                   c("_Benefit","_Risk","_BRR"))

ar_wide_ft <- ar_wide_ft %>%
  dplyr::select(Outcome, Setting, `Age group`, all_of(day_cols))

# ── 4. Replace NA strings → "beneficial" ─────────────────────────────────────
na_patterns <- c(
  "NA\n(NA–NA)", "NA\n(NA-NA)", "NA (NA–NA)", "NA (NA-NA)", "NA\n(NA NA)"
)

ar_wide_ft <- ar_wide_ft %>%
  dplyr::mutate(
    dplyr::across(
      all_of(day_cols),
      ~ dplyr::case_when(.x %in% na_patterns ~ "beneficial", TRUE ~ .x)
    )
  )

# ── 5. col_labels  ─────────────────────────────────────
col_labels <- c(
  Outcome = "Outcome", Setting = "Setting", `Age group` = "Age group"
)
for (d in c("7d","14d","30d","90d")) {
  col_labels[paste0(d, "_Benefit")] <- "Benefit:\nAverted\n(per 10,000)"
  col_labels[paste0(d, "_Risk")]    <- "Risk:\nAttributable\n(per 10,000)"
  col_labels[paste0(d, "_BRR")]     <- "BRR\n(Prevented\nper 1 caused)"
}


# ── 6. p_accept_wide 
p_accept_wide <- ceac_df %>%
  filter(threshold == 1) %>%   
  mutate(
    days      = factor(days, levels = c("7d","14d","30d","90d")),
    age_group = ifelse(age_group == "65", "65+", age_group),
    outcome   = toupper(outcome),
    setting   = stringr::str_to_title(setting)
  ) %>%
  mutate(prob_fmt = sprintf("%.1f%%", p_accept * 100)) %>%
  dplyr::select(outcome, setting, age_group, days, prob_fmt)

p_accept_wide_spread <- ceac_df %>%
  filter(threshold == 1) %>%
  mutate(
    days      = as.character(days),  # factor → character
    age_group = ifelse(age_group == "65", "65+", age_group),
    outcome   = toupper(outcome),
    setting   = stringr::str_to_title(setting)
  ) %>%
  mutate(prob_fmt = sprintf("%.1f%%", p_accept * 100)) %>%
  dplyr::select(outcome, setting, age_group, days, prob_fmt) %>%
  # 중복 있으면 첫번째만
  distinct(outcome, setting, age_group, days, .keep_all = TRUE) %>%
  rename(Outcome = outcome, Setting = setting, `Age group` = age_group) %>%
  pivot_wider(
    names_from  = days,
    values_from = prob_fmt,
    names_glue  = "{days}_Prob"
  )

# 확인 — 15행, list-col 없어야 함
glimpse(p_accept_wide_spread)

# join
ar_wide_ft <- ar_wide_ft %>%
  left_join(p_accept_wide_spread, by = c("Outcome", "Setting", "Age group"))

# ── 8. 컬럼 순서 재정렬 — days별로 Benefit/Risk/BRR/Prob 순서
day_cols_full <- paste0(
  rep(c("7d","14d","30d","90d"), each = 4),
  c("_Benefit","_Risk","_BRR","_Prob")
)

ar_wide_ft <- ar_wide_ft %>%
  dplyr::select(Outcome, Setting, `Age group`, all_of(day_cols_full))

# ar_wide_ft → long 변환 (p_accept 이미 포함돼 있음)
ar_long_table <- ar_wide_ft %>%
  pivot_longer(
    cols      = all_of(day_cols_full),
    names_to  = c("days", ".value"),
    names_sep = "_"
  ) %>%
  # NA 처리 — Prob 컬럼 NA → "N/A"
  mutate(
    Prob    = tidyr::replace_na(Prob, "N/A"),
    days    = factor(days, levels = c("7d","14d","30d","90d")),
    Setting = factor(Setting, levels = c("High","Moderate","Low"))
  ) %>%
  arrange(Outcome, Setting, `Age group`, days) %>%
  rename(`Travel duration` = days)

ar_long_table <- ar_long_table %>%
  mutate(Prob = ifelse(is.na(Prob) | Prob == "N/A", "beneficial", Prob))

ft_ar <- flextable::flextable(ar_long_table) %>%
  flextable::set_header_labels(
    Outcome           = "Outcome",
    Setting           = "Setting",
    `Age group`       = "Age group",
    `Travel duration` = "Travel\nduration\n(days)",
    Benefit           = "Benefit:\nOutcomes averted\n(per 10,000)",
    Risk              = "Risk:\nOutcomes attributable\n(per 10,000)",
    BRR               = "BRR:\nPrevented per 1\noutcome caused",
    Prob              = "Probability\n(BRR > 1)\n(%)"
  ) %>%
  flextable::theme_booktabs() %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:4, align = "left", part = "all") %>%
  flextable::merge_v(j = c("Outcome", "Setting", "Age group")) %>%
  flextable::valign(j = c("Outcome", "Setting", "Age group"), valign = "top") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit()

ft_ar 

doc <- officer::read_docx() %>%
  officer::body_add_par(
    "Benefit-Risk Summary by Outcome, Setting, Age Group, and Travel Duration",
    style = "heading 2"
  ) %>%
  flextable::body_add_flextable(ft_ar)

print(doc, target = "06_Results/BRR_table_travel.docx")