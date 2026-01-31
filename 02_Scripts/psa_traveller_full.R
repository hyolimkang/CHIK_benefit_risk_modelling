
# Define continuous attack rate grid
ar_grid <- seq(0.005, 0.80, by = 0.005)  # 0.5% to 80%, 0.5% increments
# Total: 160 different AR values

# Or more focused range
ar_grid <- seq(0.01, 0.80, by = 0.005)  # 1% to 80%, 0.5% increments
# Total: 79 different AR values

cat("Total AR scenarios:", length(ar_grid), "\n")
cat("AR range:", min(ar_grid) * 100, "% to", max(ar_grid) * 100, "%\n")

###  PSA with continuous AR ----------------------------------------

psa_out_list <- list()

# Note: Remove ar_small, ar_med, ar_large from LHS
# They are now deterministic scenarios, not uncertain parameters

# Number of entry day samples per combination (to average over entry timing uncertainty)
# Higher values = more stable estimates but slower computation
n_entry_samples <- 50  # Sample 50 random entry days per draw/age/AR/duration combination

for (d in seq_len(nrow(lhs_sample))) { 
  
  # Sample uncertain parameters from LHS
  ve_d <- lhs_sample$ve[d]
  
  p_sae_vacc_u65    <- lhs_sample$p_sae_vacc_u65[d]
  p_sae_vacc_65     <- lhs_sample$p_sae_vacc_65[d]
  p_death_vacc_u65  <- lhs_sample$p_death_vacc_u65[d]
  p_death_vacc_65   <- lhs_sample$p_death_vacc_65[d]
  
  p_sae_nat_u65     <- lhs_sample$p_sae_nat_u65[d]
  p_sae_nat_65      <- lhs_sample$p_sae_nat_65[d]
  p_death_nat_u65   <- lhs_sample$p_death_nat_u65[d]
  p_death_nat_65    <- lhs_sample$p_death_nat_65[d]
  
  epi_months        <- lhs_sample$epi_months[d]
  trav_7d           <- lhs_sample$trav_7d[d]
  trav_14d          <- lhs_sample$trav_14d[d]
  trav_30d          <- lhs_sample$trav_30d[d]
  trav_90d          <- lhs_sample$trav_90d[d]
  
  travel_days <- list(
    "7d" = trav_7d,
    "14d" = trav_14d,
    "30d" = trav_30d,
    "90d" = trav_90d
  )
  
  L_days <- max(1L, round(epi_months) * 30.5)
  
  # Loop through age groups
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
    
    # NEW: Loop through continuous AR grid
    for (AR_total in ar_grid) {
      
      # Compute daily FOI for this specific AR
      foi_daily <- compute_daily_foi(AR_total = AR_total, L = L_days)
      
      # Loop through travel durations
      for (days_label in names(travel_days)) {
        
        D <- round(travel_days[[days_label]])
        max_entry <- max(1L, L_days - D + 1L)
        
        # Sample multiple entry days to average over entry timing uncertainty
        # This reduces variance and better represents expected benefit-risk
        entry_days <- sample.int(max_entry, size = n_entry_samples, replace = TRUE)
        
        # Initialize storage for this combination
        temp_results <- list()
        
        # Compute outcomes for each entry day
        for (entry_idx in seq_along(entry_days)) {
          entry_day <- entry_days[entry_idx]
          
          # Travel-specific AR for this entry day
          AR_travel <- compute_ar_travel(foi_daily, entry_day, D)
          
          # Compute outcomes
          sae <- compute_outcome(AR_travel, p_nat = p_sae_nat, 
                                 p_vacc = p_sae_vacc, VE = ve_d)
          death <- compute_outcome(AR_travel, p_nat = p_death_nat, 
                                   p_vacc = p_death_vacc, VE = ve_d)
          
          # Store intermediate results
          temp_results[[entry_idx]] <- list(
            entry_day = entry_day,
            AR_travel = AR_travel,
            sae = sae,
            death = death
          )
        }
        
        # Average across entry day samples (using median for robustness)
        # Extract vectors for averaging
        AR_travel_vec <- sapply(temp_results, function(x) x$AR_travel)
        risk_nv_10k_sae_vec <- sapply(temp_results, function(x) x$sae$risk_nv_10k)
        risk_v_10k_sae_vec <- sapply(temp_results, function(x) x$sae$risk_v_10k)
        averted_10k_sae_vec <- sapply(temp_results, function(x) x$sae$averted_10k)
        excess_10k_sae_vec <- sapply(temp_results, function(x) x$sae$excess_10k)
        brr_sae_vec <- sapply(temp_results, function(x) x$sae$brr)
        net_10k_sae_vec <- sapply(temp_results, function(x) x$sae$net_10k)
        
        risk_nv_10k_death_vec <- sapply(temp_results, function(x) x$death$risk_nv_10k)
        risk_v_10k_death_vec <- sapply(temp_results, function(x) x$death$risk_v_10k)
        averted_10k_death_vec <- sapply(temp_results, function(x) x$death$averted_10k)
        excess_10k_death_vec <- sapply(temp_results, function(x) x$death$excess_10k)
        brr_death_vec <- sapply(temp_results, function(x) x$death$brr)
        net_10k_death_vec <- sapply(temp_results, function(x) x$death$net_10k)
        
        # Compute median across entry day samples (more robust than mean)
        # Store single averaged result per combination
        psa_out_list[[length(psa_out_list) + 1]] <- data.frame(
          draw              = d,
          age_group         = age,
          AR_total          = AR_total,          # Continuous AR
          AR_total_pct      = AR_total * 100,    # For easy reading
          days              = days_label,
          entry_day         = median(entry_days),  # Representative entry day
          AR                = median(AR_travel_vec, na.rm = TRUE),
          risk_nv_10k_sae   = median(risk_nv_10k_sae_vec, na.rm = TRUE),
          risk_v_10k_sae    = median(risk_v_10k_sae_vec, na.rm = TRUE),
          averted_10k_sae   = median(averted_10k_sae_vec, na.rm = TRUE),
          excess_10k_sae    = median(excess_10k_sae_vec, na.rm = TRUE),
          brr_sae           = median(brr_sae_vec, na.rm = TRUE),
          net_10k_sae       = median(net_10k_sae_vec, na.rm = TRUE),
          risk_nv_10k_death = median(risk_nv_10k_death_vec, na.rm = TRUE),
          risk_v_10k_death  = median(risk_v_10k_death_vec, na.rm = TRUE),
          averted_10k_death = median(averted_10k_death_vec, na.rm = TRUE),
          excess_10k_death  = median(excess_10k_death_vec, na.rm = TRUE),
          brr_death         = median(brr_death_vec, na.rm = TRUE),
          net_10k_death     = median(net_10k_death_vec, na.rm = TRUE)
        )
      }
    }
  }
  
  # Progress indicator (optional)
  if (d %% 10 == 0) {
    cat("Completed", d, "of", nrow(lhs_sample), "PSA draws\n")
  }
}

# Combine results
psa_df <- bind_rows(psa_out_list)

cat("\nTotal simulations:", nrow(psa_df), "\n")
cat("Structure:", nrow(lhs_sample), "draws ×", 
    length(ar_grid), "AR values ×",
    nrow(all_risk), "age groups ×",
    length(travel_days), "durations\n")
cat("Note: Each combination averaged over", n_entry_samples, "random entry day samples\n")
cat("Total entry day samples computed:", nrow(psa_df) * n_entry_samples, "\n")


ar_summary_all <- psa_df %>%
  pivot_longer(
    cols = c(brr_sae, brr_death, net_10k_sae, net_10k_death),
    names_to = c(".value", "outcome"),
    names_pattern = "(brr|net_10k)_(sae|death)"
  ) %>%
  mutate(
    outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death"),
    days = factor(days, levels = c("7d", "14d", "30d", "90d"))
  ) %>%
  group_by(AR_total_pct, age_group, days, outcome) %>%
  summarise(
    brr_med = median(brr, na.rm = TRUE),
    brr_lo = quantile(brr, 0.025, na.rm = TRUE),
    brr_hi = quantile(brr, 0.975, na.rm = TRUE),
    net_med = median(net_10k, na.rm = TRUE),
    net_lo = quantile(net_10k, 0.025, na.rm = TRUE),
    net_hi = quantile(net_10k, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Plot with unified Y-axis within each outcome
plot_brr_travel <- 
ggplot(ar_summary_all, 
       aes(x = AR_total_pct, y = brr_med, 
           color = age_group, fill = age_group)) +
  
  geom_ribbon(aes(ymin = brr_lo, ymax = brr_hi), 
              alpha = 0.2, color = NA) +
  
  geom_line(linewidth = 1.2) +
  
  geom_hline(yintercept = 1, linetype = "dashed", 
             color = "gray50", linewidth = 0.7) +
  
  # Fixed Y-axis within rows (outcome), but allow different scales between outcomes
  facet_grid(rows = vars(outcome), 
             cols = vars(days),
             scales = "free_y") +  # This keeps Y consistent within each row
  
  scale_color_manual(
    values = c("18-64" = "#E69F00", "65" = "#56B4E9"),
    name = "Age group"
  ) +
  scale_fill_manual(
    values = c("18-64" = "#E69F00", "65" = "#56B4E9"),
    name = "Age group"
  ) +
  
  labs(
    x = "Total outbreak attack rate (%)",
    y = "Benefit-Risk Ratio (BRR)"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.caption = element_text(hjust = 0, size = 9, color = "gray40"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right"
  ) + scale_y_continuous(
    trans = pseudo_log_trans(base = 10),
    breaks = c(0, 1, 10, 100, 1000),
    labels = scales::comma_format()
  )



### Benefit Harm trade off -----------------------------------------------------
# Prepare data
br_space_data <- psa_df %>%
  mutate(
    days = factor(days, levels = c("7d", "14d", "30d", "90d"))
  ) %>%
  pivot_longer(
    cols = c(excess_10k_sae, excess_10k_death,
             averted_10k_sae, averted_10k_death,
             net_10k_sae, net_10k_death),
    names_to = c(".value", "outcome"),
    names_pattern = "(excess_10k|averted_10k|net_10k)_(sae|death)"
  ) %>%
  mutate(
    outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death")
  ) %>%
  group_by(AR_total_pct, age_group, days, outcome) %>%
  summarise(
    x_med = median(excess_10k),
    y_med = median(averted_10k),
    x_lo = quantile(excess_10k, 0.025),
    x_hi = quantile(excess_10k, 0.975),
    y_lo = quantile(averted_10k, 0.025),
    y_hi = quantile(averted_10k, 0.975),
    net_med = median(net_10k),
    .groups = "drop"
  )

# Key AR points (with cross marks)
key_ar <- c(10, 20, 30, 40)
br_key_points <- br_space_data %>%
  filter(AR_total_pct %in% key_ar)

# REDUCED labels (only show some, not all!)
label_ar <- c(10, 40)  # Only label first and last
br_label_points <- br_key_points %>%
  filter(AR_total_pct %in% label_ar)

# Background grid
create_bg_grid <- function(data) {
  x_max <- max(data$x_hi, na.rm = TRUE) * 1.1
  y_max <- max(data$y_hi, na.rm = TRUE) * 1.1
  y_min <- min(data$y_lo, na.rm = TRUE)
  y_min <- min(y_min, -5)
  
  expand.grid(
    x = seq(0, x_max, length.out = 100),
    y = seq(y_min, y_max, length.out = 100)
  ) %>%
    mutate(net = y - x)
}

# Plot function with legend
plot_br_final <- function(data, outcome_type) {
  
  plot_data_all <- data %>% filter(outcome == outcome_type)
  plot_data_key <- br_key_points %>% filter(outcome == outcome_type)
  plot_data_label <- br_label_points %>% filter(outcome == outcome_type)
  
  plot_key_18_64 <- plot_data_key %>% filter(age_group == "18-64")
  plot_key_65 <- plot_data_key %>% filter(age_group == "65")
  
  bg_grid <- create_bg_grid(plot_data_all)
  max_abs <- max(abs(bg_grid$net), na.rm = TRUE)
  
  ggplot() +
    # Background
    geom_raster(data = bg_grid,
                aes(x = x, y = y, fill = net),
                alpha = 0.7) +
    
    scale_fill_gradientn(
      name = "Net benefit\nper 10,000",
      colours = c("#d73027", "#fee090", "white", "#e0f3f8", "#4575b4"),
      values = scales::rescale(c(-max_abs, -max_abs/2, 0, max_abs/2, max_abs)),
      limits = c(-max_abs, max_abs),
      na.value = "gray90"
    ) +
    
    # Break-even line
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", linewidth = 0.8, color = "black") +
    
    # Full trajectory (faint)
    geom_path(data = plot_data_all,
              aes(x = x_med, y = y_med, group = age_group),
              linewidth = 0.3, alpha = 0.2, color = "gray50") +
    
    # Key trajectory
    geom_path(data = plot_key_18_64,
              aes(x = x_med, y = y_med),
              color = "#E69F00", linewidth = 1.5, alpha = 0.9) +
    geom_path(data = plot_key_65,
              aes(x = x_med, y = y_med),
              color = "#56B4E9", linewidth = 1.5, alpha = 0.9) +
    
    # THINNER cross error bars (0.5 instead of 0.8)
    geom_errorbar(data = plot_key_18_64,
                  aes(x = x_med, ymin = y_lo, ymax = y_hi),
                  width = 0, linewidth = 0.5, alpha = 0.7, color = "#E69F00") +
    geom_errorbarh(data = plot_key_18_64,
                   aes(y = y_med, xmin = x_lo, xmax = x_hi),
                   height = 0, linewidth = 0.5, alpha = 0.7, color = "#E69F00") +
    
    geom_errorbar(data = plot_key_65,
                  aes(x = x_med, ymin = y_lo, ymax = y_hi),
                  width = 0, linewidth = 0.5, alpha = 0.7, color = "#56B4E9") +
    geom_errorbarh(data = plot_key_65,
                   aes(y = y_med, xmin = x_lo, xmax = x_hi),
                   height = 0, linewidth = 0.5, alpha = 0.7, color = "#56B4E9") +
    
    # SMALLER points (size 3.5 instead of 4)
    geom_point(data = plot_key_18_64,
               aes(x = x_med, y = y_med),
               shape = 21, size = 3.5, 
               fill = "#E69F00", color = "black", stroke = 0.7) +
    geom_point(data = plot_key_65,
               aes(x = x_med, y = y_med),
               shape = 24, size = 3.5, 
               fill = "#56B4E9", color = "black", stroke = 0.7) +
    
    # FEWER labels (only 10% and 40%)
    geom_text(data = plot_data_label,
              aes(x = x_med, y = y_med, label = paste0(AR_total_pct, "%")),
              nudge_y = max(bg_grid$y) * 0.04,
              size = 3.5, fontface = "bold", color = "black") +
    
    # Facets
    facet_wrap(~days, nrow = 1, scales = "free") +
    
    # Labels
    labs(
      title = paste0("Benefit-Risk trajectory: ", outcome_type),
      subtitle = "Trajectory shows AR 1-40% with key points at 10%, 20%, 30%, 40%. Crosses show 95% uncertainty.",
      x = "Vaccine-attributable severe outcome per 10,000",
      y = "Infection-attributable severe outcome per 10,000",
      caption = "Background: red = net harm, white = break-even, blue = net benefit. Diagonal line = break-even (net = 0).\nCircles = 18-64 years, Triangles = 65+ years. Labeled points show AR at 10% and 40%."
    ) +
    
    # Theme
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 11),
      plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray95"),
      legend.position = "right",
      panel.spacing = unit(1.2, "lines")
    )
}

# Create plots
plot_sae_final <- plot_br_final(br_space_data, "SAE")
plot_death_final <- plot_br_final(br_space_data, "Death")

print(plot_sae_final)
print(plot_death_final)

# Save
ggsave("br_trajectory_sae_final.png", plot_sae_final,
       width = 16, height = 5, dpi = 300, bg = "white")
ggsave("br_trajectory_death_final.png", plot_death_final,
       width = 16, height = 5, dpi = 300, bg = "white")



## v4. representative
br_space_data <- psa_df %>%
  mutate(
    days = factor(days, levels = c("7d", "14d", "30d", "90d"))
  ) %>%
  pivot_longer(
    cols = c(excess_10k_sae, excess_10k_death,
             averted_10k_sae, averted_10k_death,
             net_10k_sae, net_10k_death),
    names_to = c(".value", "outcome"),
    names_pattern = "(excess_10k|averted_10k|net_10k)_(sae|death)"
  ) %>%
  mutate(
    outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death")
  ) %>%
  group_by(AR_total_pct, age_group, days, outcome) %>%
  summarise(
    x_med = median(excess_10k),
    y_med = median(averted_10k),
    x_lo = quantile(excess_10k, 0.025),
    x_hi = quantile(excess_10k, 0.975),
    y_lo = quantile(averted_10k, 0.025),
    y_hi = quantile(averted_10k, 0.975),
    net_med = median(net_10k),
    .groups = "drop"
  )

# SELECT representative AR values (like original)
representative_ar <- c(10, 30, 50, 80)  # 10%, 30%, 50%, 80%
# Or use: c(10, 50, 80) if you have data up to 80%

br_representative <- br_space_data %>%
  filter(AR_total_pct %in% representative_ar) %>%
  mutate(
    ar_category = factor(
      paste0("AR = ", AR_total_pct, "%"),
      levels = paste0("AR = ", representative_ar, "%")
    )
  )

# Create background grid for each AR
create_bg_grid_by_ar <- function(data, ar_value) {
  data_subset <- data %>% filter(AR_total_pct == ar_value)
  
  x_max <- max(data_subset$x_hi, na.rm = TRUE) * 1.1
  y_max <- max(data_subset$y_hi, na.rm = TRUE) * 1.1
  y_min <- min(data_subset$y_lo, na.rm = TRUE)
  y_min <- min(y_min, -5)
  
  expand.grid(
    x = seq(0, x_max, length.out = 100),
    y = seq(y_min, y_max, length.out = 100),
    ar_category = paste0("AR = ", ar_value, "%")
  ) %>%
    mutate(
      net = y - x,
      ar_category = factor(ar_category, 
                           levels = paste0("AR = ", representative_ar, "%"))
    )
}

# Create background for all AR values
bg_grid_list <- lapply(representative_ar, function(ar) {
  create_bg_grid_by_ar(br_space_data, ar)
})
bg_grid_all <- bind_rows(bg_grid_list)

# Calculate global max for consistent color scale
max_abs <- max(abs(bg_grid_all$net), na.rm = TRUE)

# Plot
ggplot() +
  # Background gradient
  geom_raster(data = bg_grid_all,
              aes(x = x, y = y, fill = net),
              alpha = 0.9) +
  
  scale_fill_gradientn(
    name = "Outcomes averted\nper 10,000 (median)",
    colours = c("#d73027", "#fee090", "white", "#e0f3f8", "#4575b4"),
    values = scales::rescale(c(-max_abs, -max_abs/2, 0, max_abs/2, max_abs)),
    limits = c(-max_abs, max_abs),
    na.value = "gray90"
  ) +
  
  # Break-even line
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 0.7, color = "black") +
  
  # Error bars (crosses)
  geom_errorbar(data = br_representative,
                aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome),
                width = 0, linewidth = 0.6, alpha = 0.7) +
  
  geom_errorbarh(data = br_representative,
                 aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome),
                 height = 0, linewidth = 0.6, alpha = 0.7) +
  
  # Points
  geom_point(data = br_representative,
             aes(x = x_med, y = y_med, 
                 shape = age_group, color = outcome),
             size = 3, stroke = 0.7, fill = "white") +
  
  # Facets: AR (rows) × Days (columns)
  facet_grid(rows = vars(ar_category),
             cols = vars(days),
             scales = "free") +
  
  # Scales
  scale_color_manual(
    values = c("SAE" = "#00BFC4", "Death" = "#F8766D"),
    name = "Outcome"
  ) +
  
  scale_shape_manual(
    values = c("18-64" = 22, "65" = 24),  # Square and triangle
    name = "Age group"
  ) +
  
  # Labels
  labs(
    x = "Vaccine-attributable severe outcome per 10,000",
    y = "Infection-attributable severe outcome per 10,000",
    caption = "Lines show break-even (net = 0) at median VE (dashed).\nBackground gradient shows net benefit."
  ) +
  
  # Theme
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "right",
    legend.box = "vertical",
    panel.spacing = unit(0.8, "lines")
  )

# Risk space -------------------------------------------------------------------
# Data preparation
br_space_data_correct <- psa_df %>%
  mutate(
    days = factor(days, levels = c("7d", "14d", "30d", "90d"))
  ) %>%
  pivot_longer(
    cols = c(risk_v_10k_sae, risk_v_10k_death,
             risk_nv_10k_sae, risk_nv_10k_death,
             net_10k_sae, net_10k_death),
    names_to = c(".value", "outcome"),
    names_pattern = "(risk_v_10k|risk_nv_10k|net_10k)_(sae|death)"
  ) %>%
  mutate(
    outcome = dplyr::recode(outcome, "sae" = "SAE", "death" = "Death")
  ) %>%
  group_by(AR_total_pct, age_group, days, outcome) %>%
  summarise(
    x_med = median(risk_v_10k),
    x_lo = quantile(risk_v_10k, 0.025),
    x_hi = quantile(risk_v_10k, 0.975),
    y_med = median(risk_nv_10k),
    y_lo = quantile(risk_nv_10k, 0.025),
    y_hi = quantile(risk_nv_10k, 0.975),
    net_med = median(net_10k),
    .groups = "drop"
  )

representative_ar <- c(10, 30, 50, 80)

br_representative_correct <- br_space_data_correct %>%
  filter(AR_total_pct %in% representative_ar) %>%
  mutate(
    ar_category = factor(
      paste0("AR = ", AR_total_pct, "%"),
      levels = paste0("AR = ", representative_ar, "%")
    )
  ) %>%
  mutate(age_group = ifelse(age_group == "65", "65+", age_group))

br_representative_correct <- br_representative_correct %>%
  dplyr::rename(AgeCat = age_group)

# SIMPLE: Create ONE large background for all facets
# Get global ranges
x_max_global <- max(br_representative_correct$x_hi, na.rm = TRUE) * 1.1
y_max_global <- max(br_representative_correct$y_hi, na.rm = TRUE) * 1.1

cat("Global ranges:\n")
cat("  X: [0,", round(x_max_global, 2), "]\n")
cat("  Y: [0,", round(y_max_global, 2), "]\n")

# Create single background grid
bg_grid_unified <- expand.grid(
  x = seq(0, x_max_global, length.out = 300),
  y = seq(0, y_max_global, length.out = 300)
) %>%
  mutate(outcomes_averted = y - x)

max_averted <- max(abs(bg_grid_unified$outcomes_averted), na.rm = TRUE)

cat("\nBackground grid created:\n")
cat("  Points:", nrow(bg_grid_unified), "\n")
cat("  Outcomes averted range: [", 
    round(min(bg_grid_unified$outcomes_averted), 2), ",",
    round(max(bg_grid_unified$outcomes_averted), 2), "]\n")

# Plot
plot_risk_travel <- ggplot() +
  # Unified background (same for all facets)
  geom_raster(
    data = bg_grid_unified,
    aes(x = x, y = y, fill = outcomes_averted),
    alpha = 0.85,
    interpolate = TRUE
  ) +
  
  scale_fill_gradientn(
    name = "Outcomes averted\nper 10,000 (median)",
    colours = c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0"),
    values = scales::rescale(c(-max_averted, -max_averted/4, 0, 
                               max_averted/4, max_averted)),
    limits = c(-max_averted, max_averted),
    na.value = "gray90"
  ) +
  
  # Break-even line
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 0.7, color = "black") +
  
  # Error bars
  geom_errorbar(
    data = br_representative_correct,
    aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome),
    width = 0, linewidth = 0.6, alpha = 0.7
  ) +
  geom_errorbarh(
    data = br_representative_correct,
    aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome),
    height = 0, linewidth = 0.6, alpha = 0.7
  ) +
  
  # Points
  geom_point(
    data = br_representative_correct,
    aes(x = x_med, y = y_med, shape = AgeCat, color = outcome),
    size = 2, stroke = 0.7, fill = "white"
  ) +
  
  # CRITICAL: Free X only, Y is fixed!
  facet_grid(
    rows = vars(ar_category),
    cols = vars(days),
    scales = "free_x"  # Only X is free, Y is shared
  ) +
  
  # Scales
  scale_color_manual(
    values = c("Death" = "#F8766D", "SAE" = "#00BFC4"),
    name = "Chikungunya\nadverse events"
  ) +
  
  scale_shape_manual(
    values = c("18-64" = 22, "65+" = 24),
    name = "Age group"
  ) +
  
  guides(
    fill = guide_colorbar(order = 1, title.position = "top"),
    color = guide_legend(order = 2, title.position = "top"),
    shape = guide_legend(order = 3, title.position = "top")
  ) + 
  
  # Set Y limits explicitly (optional, for cleaner look)
  coord_cartesian(ylim = c(0, y_max_global)) +
  
  # Labels
  labs(
    x = "Risk of vaccine associated adverse events (per 10,000 individuals)",
    y = "Risk of severe outcomes following infection (per 10,000 individuals)",
  ) +
  
  # Theme
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "right",
    legend.box = "vertical",
    panel.spacing = unit(0.8, "lines")
  )

plot_risk_travel <- plot_risk_travel + 
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = '',
    tag_suffix = '.',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18, hjust = 0, vjust = 1)
    )
  )



# BRR tables -------------------------------------------------------------------
# Prepare summary with representative ARs only
representative_ar <- c(10, 30, 50, 80)

ar_summary_representative <- ar_summary_all %>%
  filter(AR_total_pct %in% representative_ar) %>%
  mutate(
    # Format BRR with glue
    brr_formatted = glue("{format(round(brr_med, 2), nsmall = 2)} ({format(round(brr_lo, 2), nsmall = 2)}, {format(round(brr_hi, 2), nsmall = 2)})"),
    # Order outcome: Death first, then SAE
    outcome = factor(outcome, levels = c("Death", "SAE")),
    # Format AR
    ar_label = paste0(AR_total_pct, "%")
  ) %>%
  arrange(outcome, AR_total_pct, age_group, days)

# Create wide format table
summary_table_wide <- ar_summary_representative %>%
  select(outcome, ar_label, age_group, days, brr_formatted) %>%
  pivot_wider(
    names_from = days,
    values_from = brr_formatted
  ) %>%
  arrange(outcome, ar_label, age_group)

# Display table
kable(summary_table_wide, 
      format = "html",
      col.names = c("Outcome", "Attack Rate", "Age Group", 
                    "7 days", "14 days", "30 days", "90 days"),
      caption = "Benefit-Risk Ratio by Outbreak Intensity and Travel Duration",
      align = c("l", "c", "c", "r", "r", "r", "r")) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    font_size = 12
  ) %>%
  pack_rows(index = table(summary_table_wide$outcome)) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2, bold = TRUE)

# Save as file
write.csv(summary_table_wide, "06_Results/brr_summary_table.csv", row.names = FALSE)


