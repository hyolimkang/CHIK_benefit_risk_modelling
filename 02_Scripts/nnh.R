library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# 1) LHS에서 백신 AE 확률 -> NNH로 변환 (NNH = 1/p)
nnh_df <- lhs_sample %>%
  transmute(
    `18-64_SAE`   = 1 / p_sae_vacc_u65,
    `65+_SAE`     = 1 / p_sae_vacc_65,
    `18-64_Death` = 1 / p_death_vacc_u65,
    `65+_Death`   = 1 / p_death_vacc_65
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "key",
    values_to = "NNH"
  ) %>%
  separate(key, into = c("age_group", "outcome"), sep = "_") %>%
  mutate(
    age_group = recode(age_group, `65+` = "65+"),
    outcome   = recode(outcome, SAE = "SAE", Death = "Death")
  ) %>%
  filter(is.finite(NNH), NNH > 0)

# 2) (선택) 점선으로 표시할 기준값: 각 분포의 중앙값(또는 너가 가진 외부 문헌 범위로 대체 가능)
vline_df <- nnh_df %>%
  group_by(outcome, age_group) %>%
  summarise(nnh_med = median(NNH, na.rm = TRUE), .groups = "drop")

# 3) Plot
p_nnh_age <- ggplot(nnh_df, aes(x = NNH, colour = age_group)) +
  stat_ecdf(geom = "step", linewidth = 1.1) +
  facet_wrap(~ outcome, scales = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(labels = comma) +
  labs(
    x = "Number vaccinated per 1 excess case (NNH)",
    y = "Cumulative probability",
    colour = "Age group"
    #title = "Uncertainty in vaccine AE risk expressed as NNH (CDF)",
    #subtitle = "Solid = CDF of NNH from LHS draws; Dashed = median NNH (per age group)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

p_nnh_age

pow10_breaks <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(NULL)
  k <- floor(log10(min(x))) : ceiling(log10(max(x)))
  10^k
}

p_nnh_age +
  scale_x_log10(
    breaks = pow10_breaks,
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(x = "NNH (log10 scale)")

###############################################################################

nnh_long <- psa_df %>%
  select(
    draw, age_group, AR_total_pct, days,
    risk_nv_10k_sae, risk_v_10k_sae,
    risk_nv_10k_death, risk_v_10k_death
  ) %>%
  # wide -> long: outcome(Death/SAE)별로 묶기
  pivot_longer(
    cols = c(risk_nv_10k_sae, risk_v_10k_sae,
             risk_nv_10k_death, risk_v_10k_death),
    names_to = c("risk_type", "outcome"),
    names_pattern = "risk_(nv|v)_10k_(sae|death)",
    values_to = "risk_10k"
  ) %>%
  mutate(
    outcome = recode(outcome,
                     "sae"   = "SAE",
                     "death" = "Death"),
    scenario = if_else(risk_type == "nv", "no_vacc", "vacc")
  ) %>%
  select(-risk_type) %>%
  # 다시 wide: no_vacc vs vacc
  pivot_wider(
    names_from = scenario,
    values_from = risk_10k
  ) %>%
  mutate(
    # risk difference (10,000명당)
    delta_10k = vacc - no_vacc,
    # contextual NNH (10,000 / |Δrisk|)
    nnh = ifelse(delta_10k != 0, 1e4 / abs(delta_10k), NA_real_),
    direction = case_when(
      delta_10k < 0 ~ "benefit",  
      delta_10k > 0 ~ "harm",     
      TRUE          ~ NA_character_
    )
  )

nnh_summary_10pct <- nnh_long %>%
  mutate(
    AR_pct_10 = round(AR_total_pct / 10) * 10,
    days = factor(days, levels = c("7d", "14d", "30d", "90d"))
  ) %>%
  group_by(age_group, outcome, AR_pct_10, days, direction) %>%
  summarise(
    nnh_med = median(nnh, na.rm = TRUE),
    nnh_lo  = quantile(nnh, 0.025, na.rm = TRUE),
    nnh_hi  = quantile(nnh, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ar_label   = paste0(AR_pct_10, "%"),
    days_label = recode(days,
                        "7d"  = "7 days",
                        "14d" = "14 days",
                        "30d" = "30 days",
                        "90d" = "90 days"),
    nnh_fmt = sprintf("%.0f (%.0f, %.0f)", nnh_med, nnh_lo, nnh_hi)
  )

nnh_net_long <- psa_df %>%
  select(
    draw, age_group, AR_total_pct, days,
    risk_nv_10k_death, risk_v_10k_death,
    risk_nv_10k_sae,   risk_v_10k_sae
  ) %>%
  pivot_longer(
    cols = starts_with("risk_"),
    names_to = c("arm", "outcome"),
    names_pattern = "^risk_(nv|v)_10k_(death|sae)$",
    values_to = "risk_10k"
  ) %>%
  mutate(
    arm     = if_else(arm == "nv", "no_vacc", "vacc"),
    outcome = recode(outcome, death = "Death", sae = "SAE")
  ) %>%
  pivot_wider(
    names_from = arm,
    values_from = risk_10k
  ) %>%
  mutate(
    delta_10k = vacc - no_vacc,            # net risk difference per 10k
    net_harm  = delta_10k > 0,             # TRUE면 net harm
    nnh_net   = ifelse(delta_10k > 0, 1e4 / delta_10k, NA_real_)
  )

nnh_net_summary_10pct <- nnh_net_long %>%
  mutate(
    AR_pct_10 = round(AR_total_pct / 10) * 10,
    days = factor(days, levels = c("7d","14d","30d","90d"))
  ) %>%
  group_by(age_group, outcome, AR_pct_10, days) %>%
  summarise(
    pr_net_harm = mean(net_harm, na.rm = TRUE),                 # P(Δrisk>0)
    nnh_med = median(nnh_net, na.rm = TRUE),                    # harm일 때만 존재
    nnh_lo  = quantile(nnh_net, 0.025, na.rm = TRUE),
    nnh_hi  = quantile(nnh_net, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ar_label = paste0(AR_pct_10, "%"),
    days_label = recode(as.character(days),
                        "7d"="7 days","14d"="14 days","30d"="30 days","90d"="90 days"),
    nnh_fmt = ifelse(is.finite(nnh_med),
                     sprintf("%.0f (%.0f, %.0f)", nnh_med, nnh_lo, nnh_hi),
                     "—")
  )

plot_df <- nnh_net_summary_10pct %>%
  mutate(
    days = factor(days, levels = c("7d","14d","30d","90d"))
  )

ggplot(plot_df, aes(x = days, y = AR_pct_10, fill = nnh_med)) +
  # 1) 값 있는 칸만 색으로 그림 (NA칸은 여기서 제외!)
  geom_tile(data = plot_df %>% filter(is.finite(nnh_med))) +
  
  # 2) NA 칸만 빗살무늬로 덮어그리기
  geom_tile_pattern(
    data = plot_df %>% filter(!is.finite(nnh_med)),
    aes(x = days, y = AR_pct_10),
    inherit.aes = FALSE,
    pattern = "stripe",
    pattern_angle = 45,
    pattern_density = 0.25,
    pattern_spacing = 0.05,
    pattern_colour = "grey50",
    fill = "grey90",
    colour = NA
  ) +
  facet_grid(outcome ~ age_group) +
  scale_x_discrete(
    limits = c("7d","14d","30d","90d"),
    labels = c("7d"="7 days","14d"="14 days","30d"="30 days","90d"="90 days")
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 10),
    labels = function(x) paste0(x, "%")
  ) +
  scale_fill_viridis_c(trans = "log10") +
  labs(
    x = "Travel duration",
    y = "Attack rate over epidemic period",
    fill = "NNH (median)",
  ) +
  theme_minimal()


plot_pr <- nnh_long %>%
  mutate(
    AR_pct_10 = floor(AR_total_pct/10)*10,
    days = factor(days, levels=c("7d","14d","30d","90d"))
  ) %>%
  filter(AR_pct_10 >= 10) %>%
  group_by(age_group, outcome, AR_pct_10, days) %>%
  summarise(pr_net_harm = mean(delta_10k > 0, na.rm=TRUE), .groups="drop") %>%
  mutate(
    ar_label = factor(paste0(AR_pct_10, "%"), levels=paste0(seq(10,100,10),"%")),
    days_label = factor(recode(as.character(days),
                               "7d"="7 days","14d"="14 days","30d"="30 days","90d"="90 days"),
                        levels=c("7 days","14 days","30 days","90 days"))
  )

ggplot(plot_pr, aes(days_label, ar_label, fill = pr_net_harm)) +
  geom_tile() +
  facet_grid(outcome ~ age_group) +
  scale_fill_viridis_c(
    direction = -1,   # ← 핵심
    limits = c(0, 1)
  ) +
  labs(
    x = "Travel duration",
    y = "Attack rate",
    fill = "Pr(net harm)",
  ) +
  theme_minimal()
