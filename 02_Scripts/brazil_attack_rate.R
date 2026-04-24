get_total_S0 <- function(sim_results, t0 = 1) {
  out <- data.frame()
  
  for (region in names(sim_results)) {
    for (VE in names(sim_results[[region]])) {
      for (VC in names(sim_results[[region]][[VE]])) {
        
        scenario_list <- sim_results[[region]][[VE]][[VC]]
        
        for (sc in seq_along(scenario_list)) {
          S <- scenario_list[[sc]]$sim_out$S   # matrix: age x time
          total_S0 <- sum(S[, t0], na.rm = TRUE)
          
          out <- rbind(out, data.frame(
            region = region,
            VE = VE,
            VC = VC,
            scenario = paste0("Scenario_", sc),
            t0 = t0,
            total_S0 = total_S0
          ))
        }
      }
    }
  }
  
  out
}

baseline_cov <- "cov50"
baseline_ve <- "VE0"
baseline_scenario <- 1L
use_linear_ar_primary <- TRUE
ar_linear_cap <- 0.99

S0_df <- get_total_S0(sim_results_vc_ixchiq_model, t0 = 1)

S0_df <- S0_df %>% distinct(region, total_S0)

build_drawkey_preburden <- function(all_draws_df, baseline_cov, baseline_ve, baseline_scenario) {
  pre_check <- all_draws_df %>%
    dplyr::group_by(.data$Region, .data$draw_id) %>%
    dplyr::summarise(n_pre = dplyr::n_distinct(.data$total_pre), .groups = "drop")

  if (any(pre_check$n_pre != 1L)) {
    stop("Non-unique total_pre detected for some Region/draw_id; cannot build deterministic pre-burden.")
  }

  row_key <- all_draws_df %>%
    dplyr::mutate(
      scenario_chr = as.character(.data$Scenario),
      is_baseline = .data$Coverage == baseline_cov &
        .data$VE == baseline_ve &
        .data$scenario_chr == as.character(baseline_scenario)
    ) %>%
    dplyr::arrange(.data$Region, .data$draw_id, dplyr::desc(.data$is_baseline), .data$Coverage, .data$VE, .data$scenario_chr) %>%
    dplyr::group_by(.data$Region, .data$draw_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  pre_key <- row_key %>%
    dplyr::transmute(Region, draw_id, total_pre_key = .data$total_pre)

  if ("rho" %in% names(row_key)) {
    rho_key <- row_key %>%
      dplyr::transmute(Region, draw_id, rho_key = as.numeric(.data$rho))
  } else if ("total_pre_true" %in% names(row_key)) {
    rho_key <- row_key %>%
      dplyr::transmute(
        Region,
        draw_id,
        rho_key = .data$total_pre / pmax(as.numeric(.data$total_pre_true), 1e-12)
      )
  } else {
    stop("Neither rho nor total_pre_true available to reconstruct draw-key pre-burden.")
  }

  out <- pre_key %>%
    dplyr::left_join(rho_key, by = c("Region", "draw_id")) %>%
    dplyr::mutate(total_pre_true_key = .data$total_pre_key / .data$rho_key)

  if (any(!is.finite(out$total_pre_true_key) | out$total_pre_true_key <= 0, na.rm = TRUE)) {
    stop("Invalid total_pre_true_key generated for some Region/draw_id.")
  }

  out
}

drawkey_pre <- build_drawkey_preburden(
  all_draws_df = all_draws_ix_true,
  baseline_cov = baseline_cov,
  baseline_ve = baseline_ve,
  baseline_scenario = baseline_scenario
)

## total symptomatic draws (rho-adjusted) rebuilt by deterministic draw key
tot_symp_by_region_draw <- drawkey_pre %>%
  dplyr::transmute(Region, draw_id, tot_symp = .data$total_pre_true_key)

# draw-level symptomatic fraction (convert symptomatic burden to inferred true infections)
symp_by_draw <- lhs_sample_young %>%  
  transmute(
    draw_id   = row_number(),
    symp_prop = as.numeric(symp_overall)
  ) %>%
  mutate(
    symp_prop = pmin(pmax(symp_prop, 1e-6), 1 - 1e-6)
  )

true_inf_by_region_draw <- tot_symp_by_region_draw %>%
  left_join(symp_by_draw, by = "draw_id") %>%
  mutate(
    tot_inf = tot_symp / symp_prop
  )

# Sanity check: pre-burden should generally be invariant across VE/Coverage/Scenario.
pre_invariance_check <- all_draws_ix_true %>%
  dplyr::group_by(Region, draw_id) %>%
  dplyr::summarise(n_pre_values = dplyr::n_distinct(total_pre), .groups = "drop")

if (any(pre_invariance_check$n_pre_values > 1, na.rm = TRUE)) {
  warning("total_pre varies across VE/Coverage/Scenario for some Region/draw_id; check source construction before AR estimation.")
}

# attack rates by draw
# H_S0 = I/S0: cumulative infection count / initial susceptibles (linear fraction; also used as
#   hazard input to Poisson survival map). AR_S0 = 1 - exp(-H_S0); gap_H_minus_AR_S0 = H_S0 - AR_S0 >= 0.
AR_draw_by_region <- true_inf_by_region_draw %>%
  left_join(S0_df, by = c("Region" = "region")) %>%
  left_join(pop_by_state, by = c("Region" = "region")) %>%
  mutate(
    # Linear cumulative incidence using initial susceptibles as denominator.
    # This can exceed 1 when inferred infections are larger than S0.
    AR_lin_S0 = .data$tot_inf / .data$total_S0,
    AR_lin_S0_pct = 100 * .data$AR_lin_S0,
    flag_inf_gt_S0 = .data$tot_inf > .data$total_S0,

    # Hazard-style transform (bounded in [0,1]); retained for backward compatibility.
    H_S0 = .data$AR_lin_S0,
    AR_S0 = 1 - exp(-.data$H_S0),
    AR_S0_pct = 100 * .data$AR_S0,
    gap_H_minus_AR_S0 = .data$H_S0 - .data$AR_S0,

    # Population attack rate with total population denominator.
    AR_lin_pop = .data$tot_inf / .data$tot_pop,
    AR_lin_pop_pct = 100 * .data$AR_lin_pop,
    AR_pop = 1 - exp(-.data$AR_lin_pop),
    AR_pop_pct = 100 * .data$AR_pop
  )

# Empirical AR distributions for LHS (setup.R)
# - draws_list_lin_S0: linear I/S0
# - draws_list_prob_S0: bounded transform 1-exp(-I/S0)
draws_list_lin_S0 <- split(AR_draw_by_region$AR_lin_S0, AR_draw_by_region$Region)
draws_list_prob_S0 <- split(AR_draw_by_region$AR_S0, AR_draw_by_region$Region)


eps <- 1e-12
draws_list_lin_S0 <- lapply(draws_list_lin_S0, function(x){
  x <- x[is.finite(x)]
  x <- x[x > 0]
  pmin(x, ar_linear_cap)
})

draws_list_prob_S0 <- lapply(draws_list_prob_S0, function(x){
  x <- x[is.finite(x)]
  x <- x[x > 0]                 
  pmin(pmax(x, eps), 1 - eps)
})

# Backward-compatible object name used downstream
draws_list <- if (use_linear_ar_primary) draws_list_lin_S0 else draws_list_prob_S0


# Regional summaries of symptomatic burden and implied infections (same draw-level symp_prop as above)
true_inf_region <- true_inf_by_region_draw %>%
  dplyr::group_by(.data$Region) %>%
  dplyr::summarise(
    tot_symp_mid = stats::median(.data$tot_symp, na.rm = TRUE),
    tot_symp_lo  = stats::quantile(.data$tot_symp, 0.025, na.rm = TRUE, names = FALSE),
    tot_symp_hi  = stats::quantile(.data$tot_symp, 0.975, na.rm = TRUE, names = FALSE),
    true_inf_mid = stats::median(.data$tot_inf, na.rm = TRUE),
    true_inf_lo  = stats::quantile(.data$tot_inf, 0.025, na.rm = TRUE, names = FALSE),
    true_inf_hi  = stats::quantile(.data$tot_inf, 0.975, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  dplyr::rename(region = Region)


ar_target_region <- true_inf_region %>%
  dplyr::left_join(pop_by_state, by = "region") %>%
  dplyr::left_join(S0_df, by = "region") %>%
  dplyr::mutate(
    AR_mid_totpop = 1 - exp(-.data$true_inf_mid / .data$tot_pop),
    AR_lo_totpop  = 1 - exp(-.data$true_inf_lo  / .data$tot_pop),
    AR_hi_totpop  = 1 - exp(-.data$true_inf_hi  / .data$tot_pop),
    # Same functional form as draw-level AR_S0 = 1 - exp(-tot_inf / total_S0)
    H_mid_S0 = .data$true_inf_mid / .data$total_S0,
    H_lo_S0  = .data$true_inf_lo  / .data$total_S0,
    H_hi_S0  = .data$true_inf_hi  / .data$total_S0,
    AR_mid_S0 = 1 - exp(-.data$true_inf_mid / .data$total_S0),
    AR_lo_S0  = 1 - exp(-.data$true_inf_lo  / .data$total_S0),
    AR_hi_S0  = 1 - exp(-.data$true_inf_hi  / .data$total_S0),
    S_prop    = .data$total_S0 / .data$tot_pop
  )

AR_summary_by_region <- AR_draw_by_region %>%
  dplyr::filter(
    !is.na(.data$AR_S0),
    is.finite(.data$AR_S0),
    .data$AR_S0 > 0,
    .data$AR_S0 < 1
  ) %>%
  dplyr::group_by(.data$Region) %>%
  dplyr::summarise(
    AR_lin_mid_S0 = stats::median(.data$AR_lin_S0, na.rm = TRUE),
    AR_lin_lo_S0  = stats::quantile(.data$AR_lin_S0, 0.025, na.rm = TRUE, names = FALSE),
    AR_lin_hi_S0  = stats::quantile(.data$AR_lin_S0, 0.975, na.rm = TRUE, names = FALSE),
    H_mid_S0 = stats::median(.data$H_S0, na.rm = TRUE),
    H_lo_S0  = stats::quantile(.data$H_S0, 0.025, na.rm = TRUE, names = FALSE),
    H_hi_S0  = stats::quantile(.data$H_S0, 0.975, na.rm = TRUE, names = FALSE),
    AR_mid_S0 = stats::median(.data$AR_S0, na.rm = TRUE),
    AR_lo_S0  = stats::quantile(.data$AR_S0, 0.025, na.rm = TRUE, names = FALSE),
    AR_hi_S0  = stats::quantile(.data$AR_S0, 0.975, na.rm = TRUE, names = FALSE),
    gap_mid_H_minus_AR = stats::median(.data$gap_H_minus_AR_S0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(region = Region)

AR_by_state <- list(
  mid = setNames(pmin(AR_summary_by_region$AR_lin_mid_S0, ar_linear_cap), AR_summary_by_region$region),
  lo  = setNames(pmin(AR_summary_by_region$AR_lin_lo_S0, ar_linear_cap), AR_summary_by_region$region),
  hi  = setNames(pmin(AR_summary_by_region$AR_lin_hi_S0, ar_linear_cap), AR_summary_by_region$region)
)

AR_by_state_prob <- lapply(AR_by_state, function(v) 1 - exp(-v))


## Supplementary figure: empirical AR_S0 = 1 - exp(-tot_inf / total_S0)
## (Avoid theme(text = element_text(family = ...)) with showtext in RStudio — tiny device + large pt.)
p <- ggplot(AR_draw_by_region, aes(x = .data$AR_S0)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30) +
  geom_density(linewidth = 0.8) +
  facet_wrap(~Region, scales = "free_y") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Attack rate (AR_S0)",
    y = "Density",
    caption = "H_S0 = tot_inf/total_S0 (linear I/S0, hazard proxy); AR_S0 = 1-exp(-H_S0); gap = H_S0-AR_S0. cov50, VE0, scenario 1."
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text       = element_text(size = 8.5, face = "bold", margin = margin(b = 2, t = 2)),
    axis.title       = element_text(size = 10, face = "bold"),
    axis.text        = element_text(size = 8),
    axis.text.x      = element_text(angle = 30, hjust = 1, size = 7.5),
    plot.caption     = element_text(size = 7, colour = "grey35", hjust = 0, lineheight = 1.15),
    panel.grid.minor = element_blank(),
    panel.spacing    = grid::unit(0.35, "lines")
  )

ggsave("06_Results/suppl_attack_rates.pdf", plot = p, width = 8, height = 6, device = cairo_pdf)


scale_foi_to_target_ar <- function(foi_daily, ar_target) {
  H0 <- sum(foi_daily, na.rm = TRUE)
  Ht <- -log(1 - ar_target)
  m  <- Ht / H0
  foi_daily * m
}

foi_daily_adj_by_state <- lapply(names(foi_daily_by_state), function(st) {
  foi0 <- foi_daily_by_state[[st]]
  list(
    mid = scale_foi_to_target_ar(foi0, AR_by_state_prob$mid[[st]]),
    lo  = scale_foi_to_target_ar(foi0, AR_by_state_prob$lo[[st]]),
    hi  = scale_foi_to_target_ar(foi0, AR_by_state_prob$hi[[st]])
  )
})
names(foi_daily_adj_by_state) <- names(foi_daily_by_state)

foi_daily_by_state_mid <- lapply(foi_daily_adj_by_state, `[[`, "mid")

states_to_run <- names(foi_daily_by_state_mid)

## end ------------------------------------------------------------------------