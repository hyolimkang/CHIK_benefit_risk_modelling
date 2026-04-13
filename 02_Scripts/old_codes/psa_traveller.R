
# uncertain params and distributions -------------------------------------------
# p_sae_vacc            | beta (uninformative)
# p_death_vacc          | beta (uninformative)
# hosp_rate <65         | beta (uninformative)
# hosp_rate >65         | beta (uninformative)
# death_rate <65        | beta (uninformative)
# death_rate >65        | beta (uninformative)
# vaccine efficacy      | beta (uninformative)
# duration of travel    | uniform (uninformative)
# duration of outbreak  | uniform (uninformative)
# entry time            | uniform (uninformative)
# exit time             | uniform (uninformative)
# attack rate small     | uniform (uninformative)
# attack rate medium    | uniform (uninformative)
# attack rate large     | uniform (uninformative)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# lhs sampling 

set.seed(123)
runs = 1000
A <- randomLHS (n = runs, 
                k = 17) 

lhs_sample <- matrix (nrow = nrow(A), ncol = ncol(A))

lhs_sample [,1]   <- qbeta(A[,1], shape1 = 7+1, shape2 = 32949-7+1)
lhs_sample [,2]   <- qbeta(A[,2], shape1 = 25+1, shape2 = 18445-25+1)
lhs_sample [,3]   <- qbeta(A[,3], shape1 = 0+1, shape2 = 32949-0+1)
lhs_sample [,4]   <- qbeta(A[,4], shape1 = 3+1, shape2 = 18445-3+1)
lhs_sample [,5]   <- qbeta(A[,5], shape1 = 11974+1, shape2 = 550749-11974+1)
lhs_sample [,6]   <- qbeta(A[,6], shape1 = 3064+1, shape2 = 110578-3064+1)
lhs_sample [,7]   <- qbeta(A[,7], shape1 = 296+1, shape2 = 550749-296+1)
lhs_sample [,8]   <- qbeta(A[,8], shape1 = 428+1, shape2 = 110578-428+1)
lhs_sample [,9]   <- qtruncnorm(A[,9], a = 0.972, b = 0.983, mean = 0.978, sd   = (0.983 - 0.972) / (2 * qnorm(0.975)))
lhs_sample [,10]  <- qunif(A[,10], min = 0.1 * 0.9, max = 0.1 * 1.1)
lhs_sample [,11]  <- qunif(A[,11], min = 0.2 * 0.9, max = 0.2 * 1.1)
lhs_sample [,12]  <- qunif(A[,12], min = 0.3 * 0.9, max = 0.3 * 1.1)
lhs_sample [,13]  <- qunif(A[,13], min = 4 * 0.9, max = 4 * 1.1)
lhs_sample [,14]  <- qunif(A[,14], min = 7 * 0.9, max = 7 * 1.1)
lhs_sample [,15]  <- qunif(A[,15], min = 14 * 0.9, max = 14 * 1.1)
lhs_sample [,16]  <- qunif(A[,16], min = 30 * 0.9, max = 30 * 1.1)
lhs_sample [,17]  <- qunif(A[,17], min = 90 * 0.9, max = 90 * 1.1)

cols <- c(
  "p_sae_vacc_u65",
  "p_sae_vacc_65",
  "p_death_vacc_u65",
  "p_death_vacc_65",
  "p_sae_nat_u65",
  "p_sae_nat_65",
  "p_death_nat_u65",
  "p_death_nat_65",
  "ve",
  "ar_small",
  "ar_med",
  "ar_large",
  "epi_months",
  "trav_7d",
  "trav_14d",
  "trav_30d",
  "trav_90d"
)
colnames (lhs_sample) <- cols
lhs_sample   <- as.data.frame(lhs_sample)

#check 
par(mfrow = c(3, 4))  
for (nm in names(lhs_sample)) {
  hist(lhs_sample[[nm]], main = nm, xlab = "", col = "skyblue")
}
par(mfrow = c(1, 1))

## functions to make daily FOI --------------------------------------------------
# make time varying FOi (daily) by the outbreak curve shape 
make_weights_foi <- function(L) {
  t <- seq(0 , 1, length.out = L) # L is the total lengths of an outbreak and the sum should be 1 for FOI 
  w <- sin(pi * t) # weight (shape similar to bell-shape curve outbreak)
  w_norm <- w / sum(w) # normalise
  return(w_norm)
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
  if(D >= L) {
    cum_foi_travel <- sum(foi_daily) # if D is large enough to cover total duration, then sum all daily FOIs
  }else{
    entry_day <- max(1L, min(entry_day, L - D + 1)) 
    idx <- entry_day: (entry_day + D - 1L)  # if entry day is too late, take the maximum entry 
    foi_travel <- sum(foi_daily[idx])
  }
  1 - exp(-foi_travel)
}
  
### Start PSA ------------------------------------------------------------------
psa_out_list <- list()

for (d in seq_len(nrow(lhs_sample))) { 
  
  # sample from each random draw
  ve_d <- lhs_sample$ve[d]
  
  p_sae_vacc_u65    <- lhs_sample$p_sae_vacc_u65[d]
  p_sae_vacc_65     <- lhs_sample$p_sae_vacc_65[d]
  p_death_vacc_u65  <- lhs_sample$p_death_vacc_u65[d]
  p_death_vacc_65   <- lhs_sample$p_death_vacc_65[d]
  
  p_sae_nat_u65     <- lhs_sample$p_sae_nat_u65[d]
  p_sae_nat_65      <- lhs_sample$p_sae_nat_65[d]
  p_death_nat_u65   <- lhs_sample$p_death_nat_u65[d]
  p_death_nat_65    <- lhs_sample$p_death_nat_65[d]
  
  ar_small          <- lhs_sample$ar_small[d]
  ar_med            <- lhs_sample$ar_med[d]
  ar_large          <- lhs_sample$ar_large[d]
  
  epi_months        <- lhs_sample$epi_months[d]
  trav_7d           <- lhs_sample$trav_7d[d]
  trav_14d          <- lhs_sample$trav_14d[d]
  trav_30d          <- lhs_sample$trav_30d[d]
  trav_90d          <- lhs_sample$trav_90d[d]
  
  # compute FOIs for each outbreak scenario (where total duration of an outbreak is also treated as an uncertain variable)
  # ar_small: 10% of total attack rate during an outbreak
  # ar_med  : 20% of total attack rate during an outbreak
  # ar_large: 30% of total attack rate during an outbreak
  #foi_levels <- c(
  #  small = -log(1 - ar_small) / (epi_months / 12),
  #  med   = -log(1 - ar_med)   / (epi_months / 12),
  #  large = -log(1 - ar_large) / (epi_months / 12)
  #)
  
  travel_days <- list(
    "7d" = trav_7d,
    "14d" = trav_14d,
    "30d" = trav_30d,
    "90d" = trav_90d
  )
  
  L_days <- max(1L, round(epi_months) * 30.5) # total days of outbreak
  
  foi_daily_list <- list(
    small = compute_daily_foi(AR_total = ar_small, L = L_days),
    med   = compute_daily_foi(AR_total = ar_med,   L = L_days),
    large = compute_daily_foi(AR_total = ar_large, L = L_days)
  )
  
  # loop through age group
  for (i in seq_len(nrow(all_risk))) {
    group <- all_risk[i, ]
    age   <- group$age_group
    
    if (grepl("65", age) && !grepl("under", age)) {
      p_sae_vacc   <- p_sae_vacc_65
      p_death_vacc <- p_death_vacc_65
      p_sae_nat    <- p_sae_nat_65
      p_death_nat  <- p_death_nat_65
    } else {
      p_sae_vacc   <- p_sae_vacc_u65
      p_death_vacc <- p_death_vacc_u65
      p_sae_nat    <- p_sae_nat_u65
      p_death_nat  <- p_death_nat_u65
    }
    
    # loop through outbreak scenario
    for (foi_label in names(foi_daily_list)) {
      
      #foi <- foi_levels[[foi_label]]
      foi_daily <- foi_daily_list[[foi_label]]
      
      AR_total <- switch(
        foi_label,
        "small" = ar_small,
        "med"   = ar_med,
        "large" = ar_large,
        NA_real_
      )
      
      # loop through different travel durations (days of exposure)
      for (days_label in names(travel_days)) {
        
        #days_value <- travel_days[[days_label]] 
        D <- round(travel_days[[days_label]]) # rounding travel days 
        
        entry_day <- sample.int(L_days - D + 1L, 1L) # minimum it can start from day1, but it can go up to L-d+1
        
        AR_travel <- compute_ar_travel(foi_daily, entry_day, D) # total attack rate during the travel window
        
        #AR_travel  <- compute_ar(lambda = foi, days = days_value)
        
        # compute SAE outcomes 
        sae <- compute_outcome(AR_travel, p_nat = p_sae_nat,   p_vacc = p_sae_vacc,   VE = ve_d)
        death <- compute_outcome(AR_travel, p_nat = p_death_nat, p_vacc = p_death_vacc, VE = ve_d)
        
        # store results
        psa_out_list[[length(psa_out_list) + 1]] <- data.frame(
          draw              = d,
          age_group         = age,
          ar_scenario       = foi_label,
          days              = days_label,
          entry_day         = entry_day,
          AR_total          = AR_total,
          AR                = AR_travel,
          risk_nv_10k_sae   = sae$risk_nv_10k,
          risk_v_10k_sae    = sae$risk_v_10k,
          averted_10k_sae   = sae$averted_10k,
          excess_10k_sae    = sae$excess_10k,
          brr_sae           = sae$brr,
          net_10k_sae       = sae$net_10k,
          risk_nv_10k_death = death$risk_nv_10k,
          risk_v_10k_death  = death$risk_v_10k,
          averted_10k_death = death$averted_10k,
          excess_10k_death  = death$excess_10k,
          brr_death         = death$brr,
          net_10k_death     = death$net_10k
        )
      }
    }
  }
}

# combined results 
psa_df <- do.call(rbind, psa_out_list)

# ---------------------------------------------------end of PSA-----------------
#-------------------------------------------------------------------------------
# plots 

ve_by_draw <- data.frame(draw = seq_len(nrow(lhs_sample)),
                         ve   = lhs_sample$ve)
psa_df <- psa_df %>% left_join(ve_by_draw, by = "draw")

# Death
death_df <- psa_df %>%
  transmute(
    draw, age_group,
    foi = ar_scenario,
    days = days,                    # "14d", "30d", "50d"
    outcome = "Death",
    x = excess_10k_death,
    y = risk_nv_10k_death,
    ve
  )

# SAE
sae_df <- psa_df %>%
  transmute(
    draw, age_group,
    foi = ar_scenario,
    days = days,
    outcome = "SAE",
    x = excess_10k_sae,
    y = risk_nv_10k_sae,
    ve
  )

plot_df <- bind_rows(death_df, sae_df)

# factors
plot_df$age_group <- factor(plot_df$age_group, levels = c("18-64", "65"))
plot_df$outcome   <- factor(plot_df$outcome,   levels = c("Death", "SAE"))

# labels for facets
foi_map <- c(small = "Small outbreak", med = "Medium outbreak", large = "Large outbreak")
plot_df <- plot_df %>%
  mutate(
    foi_lab  = factor(dplyr::recode(foi, !!!foi_map), levels = foi_map),
    days_lab = factor(days, levels = c("7d", "14d","30d","90d"),
                      labels = c("Days: 7","Days: 14","Days: 30", "Days: 90"))
  )

summary_df <- plot_df %>%
  group_by(age_group, outcome, foi_lab, days_lab) %>%
  summarise(
    x_med = median(x, na.rm = TRUE),
    x_lo  = quantile(x, 0.025, na.rm = TRUE),
    x_hi  = quantile(x, 0.975, na.rm = TRUE),
    y_med = median(y, na.rm = TRUE),
    y_lo  = quantile(y, 0.025, na.rm = TRUE),
    y_hi  = quantile(y, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


x_max <- max(plot_df$x, na.rm = TRUE)
y_max <- max(plot_df$y, na.rm = TRUE)

x_seq <- seq(0, x_max, length.out = 150)
y_seq <- seq(0, y_max, length.out = 150)

ve_q    <- quantile(lhs_sample$ve, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
ve_low  <- unname(ve_q[1]); ve_med <- unname(ve_q[2]); ve_high <- unname(ve_q[3])

# facet grid
bg_grid <- expand.grid(
  x = x_seq,
  y = y_seq,
  foi_lab  = levels(plot_df$foi_lab),
  days_lab = levels(plot_df$days_lab)
) %>%
  mutate(
    net_low  = y * ve_low  - x,
    net_med  = y * ve_med  - x,
    net_high = y * ve_high - x
  )

net_rng <- range(bg_grid$net_med, na.rm = TRUE)
max_abs <- max(abs(net_rng))

# threshold line
thr_df <- data.frame(
  x = seq(0, x_max * 1.05, length.out = 400),
  y = seq(0, x_max * 1.05, length.out = 400) / ve_med
)

ggplot() +
  # Background: median net benefit (VE uncertainty reflected)
  geom_raster(data = bg_grid, aes(x = x, y = y, fill = net_med), alpha = 0.8) +
  # Uncertainty boundaries where net=0 at VE low/high
  geom_contour(data = bg_grid, aes(x = x, y = y, z = net_low),
               breaks = 0, linetype = "dotted", linewidth = 0.5, color = "black") +
  geom_contour(data = bg_grid, aes(x = x, y = y, z = net_high),
               breaks = 0, linetype = "dotted", linewidth = 0.5, color = "black") +
  # Threshold line (median VE): y = x / VE_med
  geom_line(data = thr_df, aes(x = x, y = y),
            linetype = "dashed", linewidth = 0.7, color = "black") +
  
  # Raw PSA cloud
  #geom_point(data = plot_df,
  #         aes(x = x, y = y, shape = age_group, color = outcome),
  #           alpha = 0.08, size = 1.4, stroke = 0.2) +
  
  # 95% UI error bars around the medians
  geom_errorbar(data = summary_df,
                aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome),
                width = 0.0005, alpha = 0.8) +
  geom_errorbarh(data = summary_df,
                 aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome),
                 height = 0.0005, alpha = 0.8) +
  # Median point
  geom_point(data = summary_df,
             aes(x = x_med, y = y_med, shape = age_group, color = outcome),
             size = 3, stroke = 0.7) +
  
  facet_grid(rows = vars(foi_lab), cols = vars(days_lab), scales = "free") +
  scale_shape_manual(values = c("18-64" = 22, "65" = 24)) +
  scale_fill_gradientn(
    name   = "Outcomes averted\nper 10,000 (median)",
    colours = c("red", "white", "skyblue"),
    values  = scales::rescale(c(-max_abs, 0, max_abs)),
    limits  = c(-max_abs, max_abs)
  ) +
  labs(
    x = "Vaccine-attributable severe outcome per 10,000",
    y = "Infection-attributable severe outcome per 10,000",
    shape = "Age group",
    color = "Outcome",
    caption = "Lines show break-even (net = 0) at median VE (dashed)"
  ) +
  theme_bw(base_size = 12)

# brr
brr_long <- psa_df %>%
  transmute(
    draw,
    age_group,
    foi  = ar_scenario,     # "small" / "med" / "large"
    days = days,            # "14d" / "30d" / "50d"  
    Death = brr_death,
    SAE   = brr_sae
  ) %>%
  pivot_longer(cols = c(Death, SAE),
               names_to = "outcome", values_to = "BRR") %>%
  mutate(
    outcome   = factor(outcome, levels = c("Death","SAE")),
    age_group = factor(age_group, levels = c("18-64","65")),
    # facet label
    foi_lab = factor(dplyr::recode(foi, small = "Small outbreak",
                            med   = "Medium outbreak",
                            large = "Large outbreak"),
                     levels = c("Small outbreak","Medium outbreak","Large outbreak")),
    # x axis fix
    days = factor(days, levels = c("7d","14d","30d", "90d"))
  )

# 2) summary and 95%UIs
brr_sum <- brr_long %>%
  group_by(age_group, outcome, foi_lab, days) %>%
  summarise(
    BRR_med = median(BRR, na.rm = TRUE),
    BRR_lo  = quantile(BRR, 0.025, na.rm = TRUE),
    BRR_hi  = quantile(BRR, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(brr_sum, aes(x = days, y = BRR_med, color = age_group, group = age_group)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbar(aes(ymin = BRR_lo, ymax = BRR_hi),
                width = 0.15, alpha = 0.8, position = position_dodge(width = 0.0)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(outcome ~ foi_lab, scales = "free_y") +
  labs(
    x = "Days of exposure",
    y = "Benefit–Risk Ratio",
    color = "Age group",
    caption = "Points = median BRR, bars = 95% uncertainty intervals across PSA draws."
  ) +
  my_theme

################################################################################
### old functions 
################################################################################
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






