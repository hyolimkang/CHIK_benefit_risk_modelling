# ----------------------------------------------------------
# Load packages and source functions/other scripts
# ----------------------------------------------------------

# setwd
setwd("C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/CHIK_benefit_risk")
setwd("/Users/hyolimkang/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/CHIK_benefit_risk")

# load packages
pacman::p_load(
  dplyr, tidyr, tidyverse, ggplot2, patchwork, purrr, flextable, sf, raster, officer, viridis,
  cowplot, scales, ggpubr, lhs, reshape, truncnorm, knitr, kableExtra, glue, ggpattern
)

# Clean up from previous code / runs
rm(list = ls(all = TRUE))
options(scipen = 999)


# Data
load("01_Data/lhs_sample_young.RData")
load("01_Data/le_sample.RData")
load("01_Data/chikv_fatal_hosp_rate.RData")
load("01_Data/sim_results_vc_ixchiq_model.RData")
load("01_Data/combined_nnv_national_age_ixchiq.RData")
load("01_Data/combined_nnv_df_region_coverage_model.RData")
load("01_Data/mortality_chikv.RData")
load("01_Data/rho_df.RData")
load("01_Data/pop_by_state.RData")
load("01_Data/posterior_list.RData")
load("01_Data/preui_all.RData")
load("01_Data/postsim_vc_ixchiq_model.RData")

all_risk <- read.csv("01_Data/all_risk_four_age.csv")

# ----------------------------------------------------------
# Functions 
# ----------------------------------------------------------

# compute attack rate
compute_ar <- function(lambda, s0 = 1, days){
  s0 * (1 - exp(- lambda * days/365))
}

# compute outcome metrics
compute_outcome <- function(AR, p_hosp, p_death, p_sae_vacc, p_death_vacc, VE_hosp, VE_death) {
  
  # --- Infection outcomes per 10,000 ---
  risk_nv_hosp  <- 1e4 * AR * p_hosp
  risk_nv_death <- 1e4 * AR * p_death
  risk_nv_sae   <- risk_nv_hosp + risk_nv_death
  
  risk_v_hosp   <- 1e4 * AR * p_hosp  * (1 - VE_hosp)
  risk_v_death  <- 1e4 * AR * p_death * (1 - VE_death)
  risk_v_sae    <- risk_v_hosp + risk_v_death
  
  # --- Vaccine-related outcomes per 10,000 ---
  vacc_sae_10k   <- 1e4 * p_sae_vacc
  vacc_death_10k <- 1e4 * p_death_vacc
  
  # --- Averted (infection SAE/death only) ---
  averted_10k_sae   <- risk_nv_sae   - risk_v_sae
  averted_10k_death <- risk_nv_death - risk_v_death
  
  # --- Excess (vaccine-related) ---
  excess_10k_sae   <- vacc_sae_10k + vacc_death_10k   # SAE = hosp + death
  excess_10k_death <- vacc_death_10k                  # death-only
  
  # --- BRR ---
  brr_sae   <- ifelse(excess_10k_sae   == 0, NA, averted_10k_sae   / excess_10k_sae)
  brr_death <- ifelse(excess_10k_death == 0, NA, averted_10k_death / excess_10k_death)
  
  list(
    # needed by DALY block
    risk_nv_hosp  = risk_nv_hosp,
    risk_nv_death = risk_nv_death,
    risk_v_hosp   = risk_v_hosp,
    risk_v_death  = risk_v_death,
    
    # used for summaries
    averted_10k_sae   = averted_10k_sae,
    excess_10k_sae    = excess_10k_sae,
    brr_sae           = brr_sae,
    
    averted_10k_death = averted_10k_death,
    excess_10k_death  = excess_10k_death,
    brr_death         = brr_death
  )
}

make_weights_foi <- function(L) {
  t <- seq(0 , 1, length.out = L) # L is the total lengths of an outbreak and the sum should be 1 for FOI 
  w <- sin(pi * t) # weight (shape similar to bell-shape curve outbreak)
  w_norm <- w / sum(w) # normalise
  return(w_norm)
}

weekly_to_daily_foi <- function(phi_weekly) {
  phi_weekly <- as.numeric(phi_weekly)
  rep(phi_weekly / 7, each = 7) # 길이 364
}

compute_daily_foi <- function(AR_total, L) {
  foi_total <- -log(1 - AR_total) # cumulative FOI over the outbreak 
  w <- make_weights_foi(L)
  foi_daily <- foi_total * w
  return(foi_daily)
}

compute_ar_travel <- function(foi_daily, entry_day, D) {
  L <- length(foi_daily)
  D <- max(1L, round(D))
  
  if (D >= L) {
    foi_travel <- sum(foi_daily)
  } else {
    entry_day <- max(1L, min(entry_day, L - D + 1L))
    idx <- entry_day:(entry_day + D - 1L)
    foi_travel <- sum(foi_daily[idx])
  }
  
  1 - exp(-foi_travel)
}

# functions
extract_phi_long <- function(sim_results) {
  out <- list()
  k <- 1L
  
  states <- names(sim_results)
  for (st in states) {
    ve_names <- names(sim_results[[st]])
    for (ve in ve_names) {
      cov_names <- names(sim_results[[st]][[ve]])
      for (cv in cov_names) {
        reps <- sim_results[[st]][[ve]][[cv]]
        for (r in seq_along(reps)) {
          phi <- reps[[r]]$sim_out$phi
          
          out[[k]] <- data.frame(
            state = st,
            VE    = ve,
            cov   = cv,
            rep   = r,
            week  = seq_along(phi),
            phi   = as.numeric(phi)
          )
          k <- k + 1L
        }
      }
    }
  }
  do.call(rbind, out)
}

fn_br_space_benefit <- function(psa_data, 
                                na_rm = TRUE,
                                days_levels = c("7d", "14d", "30d", "90d")){
  
  br_space_data_benefit <- psa_data %>%
    mutate(days = factor(days, levels = c("7d", "14d", "30d", "90d"))) %>%
    pivot_longer(
      cols = c(
        excess_10k_sae, excess_10k_death,
        averted_10k_sae, averted_10k_death
      ),
      names_to = c(".value", "outcome"),
      names_pattern = "(excess_10k|averted_10k)_(sae|death)"
    ) %>%
    mutate(outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death")) %>%
    group_by(AR_total_pct, age_group, days, outcome) %>%
    summarise(
      # x-axis: vaccine risk
      x_med = median(excess_10k),
      x_lo  = quantile(excess_10k, 0.025),
      x_hi  = quantile(excess_10k, 0.975),
      
      # y-axis: PURE benefit (infection-related outcomes averted)
      y_med = median(averted_10k),
      y_lo  = quantile(averted_10k, 0.025),
      y_hi  = quantile(averted_10k, 0.975),
      
      .groups = "drop"
    )
  
  return(br_space_data_benefit = br_space_data_benefit)
  
}


fn_br_space_event <- function(psa_data, 
                              na_rm = TRUE,
                              days_levels = c("7d", "14d", "30d", "90d")) {
  
  br_space_sd <- psa_data %>%
    mutate(days = factor(days, levels = c("7d", "14d", "30d", "90d"))) %>%
    pivot_longer(
      cols = c(
        excess_10k_sae, excess_10k_death,
        averted_10k_sae, averted_10k_death
      ),
      names_to = c(".value", "outcome"),
      names_pattern = "(excess_10k|averted_10k)_(sae|death)"
    ) %>%
    mutate(outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death")) %>%
    group_by(AR_total_pct, age_group, days, outcome) %>%
    summarise(
      x_med = median(excess_10k, na.rm = TRUE),
      x_lo = quantile(excess_10k, 0.025, na.rm = TRUE),
      x_hi = quantile(excess_10k, 0.975, na.rm = TRUE),
      
      
      y_med = median(averted_10k, na.rm = TRUE),
      y_lo = quantile(averted_10k, 0.025, na.rm = TRUE),
      y_hi = quantile(averted_10k, 0.975, na.rm = TRUE),
      
      
      .groups = "drop"
    )
  
  return(br_space_sd = br_space_sd)
  
} 


fn_br_space_daly <- function(psa_data, 
                             na_rm = TRUE,
                             days_levels = c("7d", "14d", "30d", "90d")) {
  
  br_space_daly <- psa_data %>%
    mutate(days = factor(days, levels = c("7d", "14d", "30d", "90d")),
           outcome = "DALY") %>%
    group_by(AR_total_pct, age_group, days, outcome) %>%
    summarise(
      x_med = median(daly_sae, na.rm = TRUE),
      x_lo = quantile(daly_sae, 0.025, na.rm = TRUE),
      x_hi = quantile(daly_sae, 0.975, na.rm = TRUE),
      
      
      y_med = median(daly_averted, na.rm = TRUE),
      y_lo = quantile(daly_averted, 0.025, na.rm = TRUE),
      y_hi = quantile(daly_averted, 0.975, na.rm = TRUE),
      
      
      .groups = "drop"
    )
  
  return(br_space_daly = br_space_daly)
  
}


fn_panel_range <- function(br_representative_benefit){
  
  panel_ranges <- br_representative_benefit %>%
    group_by(outcome, ar_category, days) %>%
    summarise(
      x_min = min(c(0, x_lo), na.rm = TRUE),
      x_max = max(x_hi, na.rm = TRUE) * 1.05,
      y_min = min(c(0, y_lo), na.rm = TRUE),
      y_max = max(y_hi, na.rm = TRUE) * 1.05,
      .groups = "drop"
    )
  
  return(panel_ranges = panel_ranges)
  
}


fn_br_grid  <- function (panel_ranges) {
  
  bg_grid_optimized <- panel_ranges %>%
    rowwise() %>%
    do({
      panel_data <- .
      
      x_seq <- seq(panel_data$x_min, panel_data$x_max, length.out = 200)
      y_seq <- seq(panel_data$y_min, panel_data$y_max, length.out = 200)
      
      grid <- expand.grid(x = x_seq, y = y_seq)
      
      grid$brr <- with(grid, ifelse(x > 0, y / x, NA_real_))
      grid$log10_brr <- log10(grid$brr)
      
      grid$ar_level    <- panel_data$ar_level
      grid$outcome     <- panel_data$outcome
      grid$ar_category <- panel_data$ar_category
      grid$days        <- panel_data$days
      
      grid
    }) %>%
    ungroup()
  
  return(bg_grid_optimized = bg_grid_optimized)
  
}


fn_br_summ <- function(br_representative_benefit){
  
  br_summarized <- br_representative_benefit %>%
    group_by(outcome, ar_category, days, AgeCat) %>%
    summarise(
      x_med = mean(x_med, na.rm = TRUE),
      y_med = mean(y_med, na.rm = TRUE),
      
      x_lo = mean(x_lo, na.rm = TRUE),
      x_hi = mean(x_hi, na.rm = TRUE),
      y_lo = mean(y_lo, na.rm = TRUE),
      y_hi = mean(y_hi, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  return(br_summarized = br_summarized)
}

##  function for DALY estimation
# compute brr
compute_daly_one <- function(age_group,
                             deaths_10k = 0,
                             hosp_10k = 0,
                             nonhosp_symp_10k = 0,
                             symp_10k = NULL,      
                             sae_10k = 0,
                             deaths_sae_10k = 0,
                             draw_pars) {
  
  life_expectancy <- dplyr::case_when(
    age_group == "1-11"  ~ draw_pars$le_lost_1_11,
    age_group == "12-17" ~ draw_pars$le_lost_12_17,
    age_group == "18-64" ~ draw_pars$le_lost_18_64, 
    age_group == "65+"   ~ draw_pars$le_lost_65,
    TRUE ~ NA_real_  
  )
  
  # ---- split weights (acute/subacute/chronic durations) ----
  p_acute  <- draw_pars$acute
  p_subac  <- draw_pars$subac
  p_chr6m  <- draw_pars$chr6m
  p_chr12m <- draw_pars$chr12m
  p_chr30m <- draw_pars$chr30m
  
  total_p <- p_acute + p_subac + p_chr6m + p_chr12m + p_chr30m
  p_acute  <- p_acute  / total_p
  p_subac  <- p_subac  / total_p
  p_chr6m  <- p_chr6m  / total_p
  p_chr12m <- p_chr12m / total_p
  p_chr30m <- p_chr30m / total_p
  
  # ---- DW & durations ----
  dw_hosp      <- draw_pars$dw_hosp
  dw_nonhosp   <- draw_pars$dw_nonhosp
  dw_subac     <- draw_pars$dw_subac
  dw_chronic   <- draw_pars$dw_chronic
  
  dur_acute    <- draw_pars$dur_acute
  dur_nonhosp  <- draw_pars$dur_nonhosp
  dur_subac    <- draw_pars$dur_subac
  dur_6m       <- draw_pars$dur_6m
  dur_12m      <- draw_pars$dur_12m
  dur_30m      <- draw_pars$dur_30m
  
  # =========================
  # Disease (natural infection)
  # =========================
  yll_dz <- deaths_10k * life_expectancy
  
  # acute YLD (separate by hosp vs nonhosp)
  yld_dz_hosp_acute    <- hosp_10k * dw_hosp    * dur_acute
  yld_dz_nonhosp_acute <- nonhosp_symp_10k * dw_nonhosp * dur_nonhosp
  
  # chronic/subacute YLD (apply to symptomatic pool; if symp_10k not given, infer it)
  if (is.null(symp_10k)) symp_10k <- hosp_10k + nonhosp_symp_10k
  
  dz_subac_10k <- symp_10k * p_subac
  dz_6m_10k    <- symp_10k * p_chr6m
  dz_12m_10k   <- symp_10k * p_chr12m
  dz_30m_10k   <- symp_10k * p_chr30m
  
  yld_dz_subac <- dz_subac_10k * dw_subac   * dur_subac
  yld_dz_6m    <- dz_6m_10k    * dw_chronic * dur_6m
  yld_dz_12m   <- dz_12m_10k   * dw_chronic * dur_12m
  yld_dz_30m   <- dz_30m_10k   * dw_chronic * dur_30m
  
  yld_dz <- yld_dz_hosp_acute + yld_dz_nonhosp_acute +
    yld_dz_subac + yld_dz_6m + yld_dz_12m + yld_dz_30m
  
  daly_dz <- yll_dz + yld_dz
  
  # =========================
  # Vaccine SAE
  # =========================
  yll_sae <- deaths_sae_10k * life_expectancy
  
  sae_acute_10k  <- sae_10k * p_acute
  sae_subac_10k  <- sae_10k * p_subac
  sae_6m_10k     <- sae_10k * p_chr6m
  sae_12m_10k    <- sae_10k * p_chr12m
  sae_30m_10k    <- sae_10k * p_chr30m
  
  yld_sae_acute <- sae_acute_10k * dw_hosp    * dur_acute
  yld_sae_subac <- sae_subac_10k * dw_subac   * dur_subac
  yld_sae_6m    <- sae_6m_10k    * dw_chronic * dur_6m
  yld_sae_12m   <- sae_12m_10k   * dw_chronic * dur_12m
  yld_sae_30m   <- sae_30m_10k   * dw_chronic * dur_30m
  
  yld_sae <- yld_sae_acute + yld_sae_subac + yld_sae_6m + yld_sae_12m + yld_sae_30m
  daly_sae <- yll_sae + yld_sae
  
  list(
    daly_dz = daly_dz,
    yll_dz  = yll_dz,
    yld_dz  = yld_dz,
    daly_sae = daly_sae,
    yll_sae  = yll_sae,
    yld_sae  = yld_sae,
    components = list(
      # disease
      yld_dz_hosp_acute = yld_dz_hosp_acute,
      yld_dz_nonhosp_acute = yld_dz_nonhosp_acute,
      yld_dz_subac = yld_dz_subac,
      yld_dz_6m = yld_dz_6m,
      yld_dz_12m = yld_dz_12m,
      yld_dz_30m = yld_dz_30m,
      # sae
      yld_sae_acute = yld_sae_acute,
      yld_sae_subac = yld_sae_subac,
      yld_sae_6m = yld_sae_6m,
      yld_sae_12m = yld_sae_12m,
      yld_sae_30m = yld_sae_30m
    )
  )
}

## br space traveller function
plot_brr_outcome <- function(br_summarized, bg_grid_optimized, 
                             target_outcome, title_text, color_val,
                             show_prop = TRUE,
                             eps_x = 1e-9,          # Avoid x==0 in background grid (prevents undefined BRR/log)
                             keep_zero_axis = TRUE  # Force axes to start at 0 (prevents "looks negative")
) {
  
  plot_data <- br_summarized %>% 
    dplyr::filter(.data$outcome == target_outcome)
  
  plot_bg <- bg_grid_optimized %>% 
    dplyr::filter(.data$outcome == target_outcome) %>%
    # Remove x==0 line only for the background raster to avoid undefined BRR/log at the boundary
    dplyr::filter(.data$x > eps_x)
  
  panel_prop <- plot_bg %>%
    dplyr::mutate(is_fav = !is.na(log10_brr) & is.finite(log10_brr) & log10_brr > 0) %>%
    dplyr::group_by(ar_category, days) %>%
    dplyr::summarise(prop_fav = mean(is_fav), .groups = "drop") %>%
    dplyr::mutate(label = ifelse(prop_fav < 0.005, "BRR>1: <1%",
                                 sprintf("BRR>1: %.0f%%", 100 * prop_fav)))
  
  panel_ranges <- plot_bg %>%
    dplyr::group_by(ar_category, days) %>%
    dplyr::summarise(
      x_min = min(x, na.rm = TRUE),
      x_max = max(x, na.rm = TRUE),
      y_min = min(y, na.rm = TRUE),
      y_max = max(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      x_lab = x_min + 0.02 * (x_max - x_min),
      y_lab = y_max - 0.04 * (y_max - y_min)
    )
  
  panel_prop <- panel_prop %>%
    dplyr::left_join(panel_ranges, by = c("ar_category", "days"))
  
  p <- ggplot() +
    geom_raster(
      data = plot_bg,
      aes(x = x, y = y, fill = log10_brr),
      interpolate = FALSE,  # Turn off interpolation to reduce boundary artefacts near x≈0
      alpha = 0.85
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
    # Use identity positioning (no dodging) to avoid shifting points below/left of 0
    geom_errorbar(
      data = plot_data,
      aes(x = x_med, ymin = y_lo, ymax = y_hi),
      color = color_val, width = 0, linewidth = 0.6,
      position = "identity"
    ) +
    geom_errorbarh(
      data = plot_data,
      aes(y = y_med, xmin = x_lo, xmax = x_hi),
      color = color_val, height = 0, linewidth = 0.6,
      position = "identity"
    ) +
    geom_point(
      data = plot_data,
      aes(x = x_med, y = y_med, shape = AgeCat),
      fill = "white", color = color_val, size = 3, stroke = 1,
      position = "identity"
    ) +
    {if (show_prop)
      geom_label(
        data = panel_prop,
        aes(x = -Inf, y = Inf, label = label),
        hjust = 0, vjust = 1,
        size = 3,
        inherit.aes = FALSE,
        label.size = 0.25,
        fill = "white", alpha = 0.75,
        colour = "black"
      )
    } +
    facet_wrap(~ ar_category + days, scales = "free", ncol = 4) +
    scale_shape_manual(values = c("1-11"=21, "12-17"=22, "18-64"=23, "65+"=24),
                       name = "Age group") +
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
  
  if (keep_zero_axis) {
    p <- p + coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), expand = FALSE)
  }
  
  return(p)
}


############################################ end of function ###################

##------------------------------------------------------------------------------
## LHS Samples 
##------------------------------------------------------------------------------
set.seed(123)
runs = 1000
region_key <- c(
  "Ceará"="ce","Bahia"="bh","Paraíba"="pa","Pernambuco"="pn",
  "Rio Grande do Norte"="rg","Piauí"="pi","Tocantins"="tc",
  "Alagoas"="ag","Minas Gerais"="mg","Sergipe"="se","Goiás"="go"
)

A <- randomLHS(n = runs, k = 51 + length(region_key))
lhs_sample <- matrix(NA_real_, nrow = nrow(A), ncol = ncol(A))


# vacc sae and death
lhs_sample [,1]   <- qbeta(A[,1], shape1 = 6+0.5, shape2 = 32949-6+0.5) # conservative values
lhs_sample [,2]   <- qbeta(A[,2], shape1 = 19+0.5, shape2 = 18445-19+0.5) # conservative values
lhs_sample [,3]   <- qbeta(A[,3], shape1 = 0+0.5, shape2 = 32949-0+0.5) # conservative values
lhs_sample [,4]   <- qbeta(A[,4], shape1 = 1+0.5, shape2 = 18445-1+0.5) # conservative values

# natural hospitalisation (case + 1 / n - case + 1)
lhs_sample [,5]   <- qbeta(A[,5], shape1 = 3703+1, shape2 = 67683-3703+1)
lhs_sample [,6]   <- qbeta(A[,6], shape1 = 1568+1, shape2 = 61390-1568+1)
lhs_sample [,7]   <- qbeta(A[,7], shape1 = 11344+1, shape2 = 696504-11344+1)
lhs_sample [,8]   <- qbeta(A[,8], shape1 = 3690+1, shape2 = 113736-3690+1)

# natural death
lhs_sample [,9]   <- qbeta(A[,9], shape1 = 32+1, shape2 = 67683-32+1)
lhs_sample [,10]   <- qbeta(A[,10], shape1 = 22+1, shape2 = 61390-22+1)
lhs_sample [,11]   <- qbeta(A[,11], shape1 = 312+1, shape2 = 696504-312+1)
lhs_sample [,12]   <- qbeta(A[,12], shape1 = 534+1, shape2 = 113736-534+1)

# ve, ar, travel duration 
lhs_sample [,13]   <- qtruncnorm(A[,13], a = 0.967, b = 0.998, mean = 0.989, sd   = (0.998 - 0.967) / (2 * qnorm(0.975)))
lhs_sample [,14]  <- qunif(A[,14], min = 0.1 * 0.9, max = 0.1 * 1.1)
lhs_sample [,15]  <- qunif(A[,15], min = 0.2 * 0.9, max = 0.2 * 1.1)
lhs_sample [,16]  <- qunif(A[,16], min = 0.3 * 0.9, max = 0.3 * 1.1)
lhs_sample [,17]  <- qunif(A[,17], min = 9 * 0.9, max = 9 * 1.1)
lhs_sample [,18]  <- qunif(A[,18], min = 7 * 0.9, max = 7 * 1.1)
lhs_sample [,19]  <- qunif(A[,19], min = 14 * 0.9, max = 14 * 1.1)
lhs_sample [,20]  <- qunif(A[,20], min = 30 * 0.9, max = 30 * 1.1)
lhs_sample [,21]  <- qunif(A[,21], min = 90 * 0.9, max = 90 * 1.1)

# symp and hosps 
lhs_sample [,22]   <- qbeta (p = A[,22], shape1 = 49.14034, shape2 = 34.14837, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,23]   <- qbeta (p = A[,23], shape1 = 34.21298, shape2 = 31.963, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,24]   <- qbeta (p = A[,24], shape1 = 36.77819, shape2 = 32.7458, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,25]   <- qbeta (p = A[,25], shape1 = 35.84287, shape2 = 32.55955, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,26]   <- qbeta (p = A[,26], shape1 = 539.2823, shape2 = 14152.25, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,27]   <- qbeta (p = A[,27], shape1 = 58.96698, shape2 = 1415.207, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,28]   <- qbeta (p = A[,28], shape1 = 115.4225, shape2 = 111.383, ncp=0, lower.tail = TRUE, log.p = FALSE)

# life expectancy by age group
lhs_sample [,29]   <- qlnorm (p = A[,29], meanlog = 4.26127, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,30]   <- qlnorm (p = A[,30], meanlog = 4.119037, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,31]   <- qlnorm (p = A[,31], meanlog = 3.583519, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,32]   <- qlnorm (p = A[,32], meanlog = 1.808289, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)

# dws and durations 
lhs_sample [,33]   <- qbeta (p = A[,33], shape1 = 4.581639, shape2 = 9.87148, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,34]  <- qlnorm (p = A[,34], meanlog = -0.6301724, sdlog = 0.0852, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,35]  <- qbeta (p = A[,35], shape1 = 22.51835, shape2 = 146.7926, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,36]  <- qlnorm (p = A[,36], meanlog =  -3.734278, sdlog =  0.3877361, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,37]  <- qbeta (p = A[,37], shape1 = 21.45106, shape2 = 399.158, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,38]  <- qlnorm (p = A[,38], meanlog = -4.108138, sdlog =  0.5998406, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,39]   <- qbeta (p = A[,39], shape1 = 393.9252, shape2 = 547.0268, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,40]   <- qbeta (p = A[,40], shape1 = 17875.92, shape2 = 42754.63, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,41]   <- qbeta (p = A[,41], shape1 = 799.2911, shape2 = 3306.976, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,42]   <- qbeta (p = A[,42], shape1 = 77.56885, shape2 = 836.687, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,43]   <- qbeta (p = A[,43], shape1 = 7.605944, shape2 = 1074.944, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,44]  <- qlnorm (p = A[,44], meanlog = -2.145581, sdlog =  0.1815621, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,45]  <- qlnorm (p = A[,45], meanlog = -0.5430045, sdlog =  0.154684, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,46]  <- qlnorm (p = A[,46], meanlog = -2.262311, sdlog =  0.4746817, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,47]  <- qlnorm (p = A[,47], meanlog = -1.148854, sdlog = 0.1815042, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,48]  <- qlnorm (p = A[,48], meanlog =  -0.6931472, sdlog =  0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,49]  <- qlnorm (p = A[,49], meanlog = 0, sdlog =  0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,50]  <- qlnorm (p = A[,50], meanlog = 0.6931472, sdlog = 0.08380206, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,51]   <- qbeta (p = A[,51], shape1 = 164.6044, shape2 = 618335.4, ncp=0, lower.tail = TRUE, log.p = FALSE)


cols <- c(
  "p_sae_vacc_u65",
  "p_sae_vacc_65",
  "p_death_vacc_u65",
  "p_death_vacc_65",
  
  "p_sae_nat_11",
  "p_sae_nat_17",
  "p_sae_nat_64",
  "p_sae_nat_65",
  
  "p_death_nat_11",
  "p_death_nat_17",
  "p_death_nat_64",
  "p_death_nat_65",
  
  "ve",
  "ar_small",
  "ar_med",
  "ar_large",
  "epi_months",
  "trav_7d",
  "trav_14d",
  "trav_30d",
  "trav_90d",
  
  "symp_asia", 
  "symp_africa", 
  "symp_america", 
  "symp_overall", 
  "fatal_hosp", 
  "hosp", 
  "lt", 
  
  "le_lost_1_11",
  "le_lost_12_17",
  "le_lost_18_64",
  "le_lost_65",
  
  "dw_chronic", 
  "dur_chronic", 
  "dw_hosp", 
  "dur_acute", 
  "dw_nonhosp", 
  "dur_nonhosp",
  "acute", 
  "subac", 
  "chr6m", 
  "chr12m", 
  "chr30m", 
  "dw_chronic_mild", 
  "dw_chronic_severe", 
  "dur_subac", 
  "dw_subac", 
  "dur_6m", 
  "dur_12m", 
  "dur_30m", 
  "fatal_nonhosp"
)

draws_list <- split(AR_draw_by_region$AR_S0, AR_draw_by_region$Region)

eps <- 1e-12
draws_list <- lapply(draws_list, function(x){
  x <- x[is.finite(x)]
  pmin(pmax(x, eps), 1 - eps)
})

regs <- names(region_key)
start_col <- 52

for (i in seq_along(regs)) {
  reg <- regs[i]
  u   <- A[, 51 + i]   
  
  lhs_sample[, start_col + i - 1] <- as.numeric(
    quantile(draws_list[[reg]], probs = u, type = 8, na.rm = TRUE)
  )
}

cols_ar <- paste0("ar_", unname(region_key))  # ar_ce ... ar_go
cols2 <- c(cols, cols_ar)

colnames(lhs_sample) <- cols2
lhs_sample <- as.data.frame(lhs_sample)


# tables
digits_prob_default <- 1
digits_prob_small   <- 4
digits_other <- 2
summ_one_fmt <- function(x, name,
                         digits_prob_default = 1,
                         digits_prob_small   = 4,
                         digits_other        = 4) {
  x <- x[is.finite(x)]
  qs <- quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE, type = 8)
  
  is_pct <- grepl("^(p_|ve$|ar_)", name)
  
  small_prob <- grepl("^p_.*_vacc_", name) || grepl("^p_.*vacc", name) || grepl("^p_death_", name)
  
  if (is_pct) {
    qs <- qs * 100
    d <- if (small_prob) digits_prob_small else digits_prob_default
    fmt <- function(v) formatC(v, format = "f", digits = d)
    paste0(fmt(qs[2]), "% (95% UI: ", fmt(qs[1]), "% - ", fmt(qs[3]), "%)")
  } else {
    fmt <- function(v) formatC(v, format = "f", digits = digits_other)
    paste0(fmt(qs[2]), " (95% UI: ", fmt(qs[1]), " - ", fmt(qs[3]), ")")
  }
}

summary_tbl <- data.frame(
  parameter = names(lhs_sample),
  value = vapply(
    names(lhs_sample),
    FUN = function(nm) summ_one_fmt(lhs_sample[[nm]], nm,
                                    digits_prob_default, digits_prob_small, digits_other),
    FUN.VALUE = character(1)
  ),
  row.names = NULL,
  stringsAsFactors = FALSE
)

subset(summary_tbl, grepl("^p_death_", parameter))

write.csv(summary_tbl, "06_Results/lhs_sample_summary_95ui.csv", row.names = FALSE)