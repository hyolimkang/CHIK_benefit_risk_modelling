
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)

# log scale
plot_combined_risk_acceptability <- function(psa_df) {
  
  # 1. 데이터 통합 및 Long Format 변환
  psa_long <- psa_df %>%
    dplyr::select(draw, age_group, excess_10k_death, excess_10k_sae, daly_sae) %>%
    pivot_longer(
      cols = c(excess_10k_death, excess_10k_sae, daly_sae),
      names_to = "outcome",
      values_to = "val_to_plot"
    ) %>%
    mutate(outcome = case_when(
      outcome == "excess_10k_death" ~ "Death",
      outcome == "excess_10k_sae"   ~ "SAE",
      outcome == "daly_sae"         ~ "DALY",
      TRUE ~ outcome
    )) %>%
    mutate(val_to_plot = pmax(0, val_to_plot, na.rm = TRUE))
  
  # 2. 아웃컴별 Acceptability 계산 (안전한 중첩 방식)
  accept_results <- psa_long %>%
    group_by(outcome, age_group) %>%
    nest() %>% # 아웃컴/연령별로 데이터를 묶음
    mutate(results = map(data, function(df_sub) {
      
      max_val <- max(df_sub$val_to_plot, na.rm = TRUE)
      # X축 범위 설정 (1 in 10 ~ 1 in 1,000,000+)
      upper_limit <- if(max_val > 0) max(1000000, 10000 / (max_val / 50)) else 1000000
      denom_seq <- exp(seq(log(10), log(upper_limit), length.out = 100))
      
      # 각 분모(denom)에 대해 확률 계산
      map_dfr(denom_seq, function(d) {
        threshold_val <- 10000 / d
        data.frame(
          denom = d,
          prob_acceptable = mean(df_sub$val_to_plot <= threshold_val) * 100
        )
      })
    })) %>%
    dplyr::select(-data) %>%
    unnest(results)
  
  # 3. 시각화
  p <- ggplot(accept_results, aes(x = denom, y = prob_acceptable, color = age_group)) +
    geom_line(linewidth = 1.1, na.rm = TRUE) +
    facet_wrap(~outcome, scales = "free_x", ncol = 3) + 
    scale_x_log10(labels = label_comma(), breaks = trans_breaks("log10", function(x) 10^x)) + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 2)) +
    labs(
      title = "Vaccine Safety Acceptability Curves",
      subtitle = paste("Based on", length(unique(psa_df$draw)), "Probabilistic Draws"),
      x = "Acceptable Risk Threshold (1 in X doses)",
      y = "% of Runs with Risk <= Threshold",
      color = "Age Group"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90"),
      text = element_text(family = "Calibri")
    )
  
  return(p)
}

p_combined <- plot_combined_risk_acceptability(psa_df)
print(p_combined)





plot_combined_risk_acceptability_linear <- function(psa_df) {
  
  # 1. 데이터 통합 및 Long Format 변환
  psa_long <- psa_df %>%
    dplyr::select(draw, age_group, excess_10k_death, excess_10k_sae, daly_sae) %>%
    pivot_longer(
      cols = c(excess_10k_death, excess_10k_sae, daly_sae),
      names_to = "outcome_raw",
      values_to = "val_to_plot"
    ) %>%
    mutate(outcome = case_when(
      outcome_raw == "excess_10k_death" ~ "Death",
      outcome_raw == "excess_10k_sae"   ~ "SAE",
      outcome_raw == "daly_sae"         ~ "DALY",
      TRUE ~ outcome_raw
    )) %>%
    mutate(val_to_plot = pmax(0, val_to_plot, na.rm = TRUE))
  
  # 2. 아웃컴별 Acceptability 계산 (자연 스케일 최적화)
  accept_results <- psa_long %>%
    group_by(outcome, age_group) %>%
    nest() %>%
    # map2를 사용하여 데이터와 그룹명(outcome)을 동시에 전달 (에러 방지)
    mutate(results = map2(data, outcome, function(df_sub, out_name) {
      
      # 데이터의 상위 리스크 확인 (곡선이 떨어지는 지점을 찾기 위함)
      max_val <- max(df_sub$val_to_plot, na.rm = TRUE)
      q99_val <- quantile(df_sub$val_to_plot, 0.99, na.rm = TRUE)
      
      # 자연 스케일에서 곡선이 바닥으로 내려가는 것을 보기 위한 분모(X) 계산
      # 리스크가 0에 가깝다면 매우 큰 분모가 필요함
      target_upper <- if(q99_val > 0) (10000 / (q99_val / 5)) else 100000
      
      # 아웃컴별 최소 보장 범위 설정 (자연 스케일 가독성용)
      upper_limit <- case_when(
        out_name == "Death" ~ max(200000, target_upper), # 사망은 최소 20만까지
        out_name == "SAE"   ~ max(20000, target_upper),  # SAE는 최소 2만까지
        TRUE                ~ max(10000, target_upper)   # DALY 등
      )
      
      # 선형 간격으로 분모 생성
      denom_seq <- seq(100, upper_limit, length.out = 120)
      
      map_dfr(denom_seq, function(d) {
        threshold_val <- 10000 / d
        data.frame(
          denom = d,
          prob_acceptable = mean(df_sub$val_to_plot <= threshold_val) * 100
        )
      })
    })) %>%
    dplyr::select(-data) %>%
    unnest(results)
  
  # 3. 시각화 (Linear Scale & 스타일 적용)
  p <- ggplot(accept_results, aes(x = denom, y = prob_acceptable, color = age_group)) +
    # 선 굵기를 0.7~0.8 정도로 얇게 조절
    geom_line(linewidth = 0.75, na.rm = TRUE) +
    # 자연 스케일 적용 (scales = "free_x"로 아웃컴별 최적화)
    facet_wrap(~outcome, scales = "free_x", ncol = 3) + 
    scale_x_continuous(labels = label_comma(), expand = c(0.02, 0)) + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 2)) +
    labs(
      title = "Vaccine Safety Acceptability Curves",
      subtitle = paste("Linear Scale | Based on", length(unique(psa_df$draw)), "Probabilistic Draws"),
      x = "Acceptable Risk Threshold (1 in X doses)",
      y = "% of Simulations Satisfying Threshold",
      color = "Age Group"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      # 레전드 위치 오른쪽
      legend.position = "right",
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold"),
      text = element_text(family = "Calibri")
    )
  
  return(p)
}

# 함수 실행
p_combined <- plot_combined_risk_acceptability_linear(psa_df)
print(p_combined)

