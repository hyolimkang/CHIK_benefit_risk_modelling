# setting ----------------------------------------------------------------------
setwd("/Users/hyolimkang/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/CHIK_benefit_risk")
setwd("C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/CHIK_benefit_risk")
options(scipen = 999)
load("01_Data/mortality_chikv.RData")
all_risk <- read.csv("01_Data/all_risk.csv")

my_theme = theme_bw() +
  theme(
    legend.text = element_text(size = 12),
    legend.spacing.y = unit(0, "cm"),
    legend.spacing.x = unit(0, "cm"),
    legend.title = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.caption = element_text(hjust = 0, size = 9) # left sided
  ) 

# Basic parameters and assumptions 
# La reunion outbreak duration: Jan 2024 - May 2025 (4 months)

all_risk <- read.csv("01_Data/all_risk.csv")
ve         <- 0.989

foi_epidemic_small <- -log(1-0.1)/(4/12) # 4 month attack rate = 0.1
foi_epidemic_med   <- -log(1-0.2)/(4/12) # 4 month attack rate = 0.2
foi_epidemic_large <- -log(1-0.3)/(4/12) # 4 month attack rate = 0.3

foi_level  <- c(foi_epidemic_small, foi_epidemic_med, foi_epidemic_large)
T_days_vec     <- c(14, 30, 50)

# endemic 
T_days_vec_endemic <- 365
foi_level_endemic  <- c(0.01, 0.05, 0.1)

# Ixchiq benefit-risk analysis -------------------------------------------------

# Data -------------------------------------------------------------------------
# 1. Vaccine induced deaths (x_vax_death)
# 2. Vaccine induced severe adverse events (x_vax_sae)
# 3. Total number of vaccinated (n_dose)

# Senario ----------------------------------------------------------------------
# 1. Traveller vaccine
## assumes everyone is susceptible
## travel duration: 14/30/50 d
## destination outbreak size: small/medium/large (10%, 20%, 30% attack rate during 4 months)

# 2. Outbreak response vaccine
## assumes underlying immunity: low/medium/high
## assumes outbreak duration for 4-6 months (full exposure)
## vaccination duration for 10 weeks 

# 3. Endemic scenario 

# Unceratiny -------------------------------------------------------------------
# 1. parameter uncertainty
# 2. traveller's visiting point during an outbreak (visit-exit)
# 3. outbreak immunisation start-end

# Analysis ---------------------------------------------------------------------
# 1. long-term average annual FOI to Attack rate 
## ar = S0 * (1 - exp(- sum(lambda_t)))
## for traveller vaccine, assume that S0 = 1 (fully susceptible to CHIKV)

# 2. deaths averted by vaccination (per 10,000)
# ve
# ar
# IFR (infection fatality)
# hospitalisation (infection - hospitalisation)

# 3. excess deaths (per 10,000)
# deaths attributable to vaccine (= x_vax_death)

# 4. benefit risk ratio (BRR)
# BRR = deaths averted per 10k / excess deaths per 10k
# BRR > 1 = benefit greater than risk 
# if brr=1 

# 5. NNV to avert a death vs. NNH (vacc induced death)

################################################################################
# functions 

# compute attack rate
compute_ar <- function(lambda, s0 = 1, days){
  s0 * (1 - exp(- lambda * days/365))
}

# compute outcome metrics
compute_outcome <- function(AR, p_nat, p_vacc, VE) {
  
  # 1. Risk in unvaccinated individuals (per 10,000)
  risk_nv_10k <- 1e4 * AR * p_nat
  
  # 2. Risk in vaccinated individuals (per 10,000)
  risk_v_10k  <- 1e4 * (AR * p_nat * (1 - VE)) + 1e4 * p_vacc
  
  # 3. Averted outcomes (infection prevented thanks to vaccine)
  averted_10k <- risk_nv_10k - risk_v_10k   
  
  # 4. Excess risk caused by vaccine
  excess_10k  <- 1e4 * p_vacc
  
  # 5. Benefit-Risk Ratio
  brr         <- ifelse(excess_10k == 0, NA, averted_10k / excess_10k)
  
  # 6. Net benefit (difference in cases prevented vs caused)
  net_10k     <- averted_10k - excess_10k
  
  list(
    risk_nv_10k = risk_nv_10k,
    risk_v_10k  = risk_v_10k,
    averted_10k = averted_10k,
    excess_10k  = excess_10k,
    brr         = brr,
    net_10k     = net_10k
  )
}

# result
out <- list()

# Create result list
out_list <- list()

for (i in seq_len(nrow(all_risk))) {
  group <- all_risk[i, ]
  age <- group$age_group
  
  for (foi in foi_level) {
    for (days in T_days_vec) {
      
      AR <- compute_ar(foi, days = days)
      
      # SAE outcome
      sae <- compute_outcome(AR, group$p_sae_nat, group$p_sae_vacc, VE = ve)
      
      # Death outcome
      death <- compute_outcome(AR, group$p_death_nat, group$p_death_vacc, VE = ve)
      
      # Save (add 'days' column)
      out_list[[length(out_list) + 1]] <- data.frame(
        age_group         = age,
        foi               = foi,
        days              = days,
        AR                = AR,
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
# Combine results
result_df <- do.call(rbind, out_list)


## graphs
# ---- 1) Prepare long data for plotting (use your existing result_df) ----
# x = vaccine-attributable risk per 10k (excess), y = infection risk per 10k (unvaccinated)
death_df <- result_df %>%
  transmute(
    age_group, foi, days,
    outcome = "Death",
    x = excess_10k_death,
    y = risk_nv_10k_death
  )

sae_df <- result_df %>%
  transmute(
    age_group, foi, days,
    outcome = "SAE",
    x = excess_10k_sae,
    y = risk_nv_10k_sae
  )

plot_df <- bind_rows(death_df, sae_df)

# Ensure factor ordering if needed
plot_df$age_group <- factor(plot_df$age_group, levels = c("18-64", "65"))
plot_df$outcome   <- factor(plot_df$outcome,   levels = c("Death", "SAE"))

foi_levels <- sort(unique(plot_df$foi))
foi_labels <- setNames(c("Small outbreak","Medium outbreak","Large outbreak"),
                       as.character(foi_levels))
plot_df <- plot_df %>%
  mutate(foi_lab = factor(foi, 
                          levels = foi_levels,
                          labels = foi_labels))%>%
  mutate(days_lab = factor(paste0("Days: ", days)))

# ---- 2) Build threshold line y = x / VE across full x-range ----
# Use your VE defined earlier
x_max <- max(plot_df$x, na.rm = TRUE)
y_max <- max(plot_df$y, na.rm = TRUE)

bg_grid <- expand.grid(
  x = seq(0, x_max, length.out = 150),
  y = seq(0, y_max, length.out = 150),
  foi_lab = levels(plot_df$foi_lab),
  days    = sort(unique(plot_df$days))
) %>%
  mutate(net_per10k = ve * y - x)%>%
  mutate(days_lab = factor(paste0("Days: ", days)))

# threshold: benefit = risk
# benefit from vaccine = VE * y 
# risk from vaccine = x 
# VE * y = x
# y = x / VE 
thr_df <- data.frame(
  x = seq(0, max(plot_df$x, na.rm = TRUE) * 1.05, length.out = 400)
) %>%
  mutate(y = x / ve)

# ---- 3) Plot: FOI facets, threshold line, shapes by age, colors by outcome ----
net_rng  <- range(bg_grid$net_per10k, na.rm = TRUE)
max_abs  <- max(abs(net_rng))

ggplot()+
  geom_raster(data = bg_grid, aes(x = x, y = y, fill = net_per10k), alpha = 0.75)+
  geom_line(data = thr_df, aes(x = x, y = y),
            linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_point(data = plot_df,
             aes(x = x, y = y, shape = age_group, color = outcome),
             size = 3, stroke = 0.7)+
  facet_grid(rows = vars(foi_lab), cols = vars(days_lab), scales = "free", labeller = label_value)+
  scale_shape_manual(values = c("18-64" = 22, "65" = 24)) + 
  scale_fill_gradientn(
    name   = "Outcomes averted\nper 10,000",
    colours = c("red", "white", "skyblue"),    
    values  = scales::rescale(c(-max_abs, 0, max_abs)),
    limits  = c(-max_abs, max_abs)
  ) +
  labs(
    x = "Vaccine-attributable severe outcome per 10,000",
    y = "Infection-attributable severe outcome per 10,000",
    shape = "Age group",
    color = "Outcome",
    caption = "Epidemic situation assumed 6 months of outbreak duration. Small outbreak assumed 6 months attack rate = 10%, Medium outbreak for 20%, Large outbreak for 30%."
  )+
  my_theme
  
# for brr
brr_long <- result_df %>%
  transmute(
    age_group, foi, days,
    `Death` = brr_death,
    `SAE`   = brr_sae
  ) %>%
  pivot_longer(cols = c(`Death`, `SAE`), names_to = "outcome", values_to = "BRR") %>%
  mutate(
    outcome = factor(outcome, levels = c("Death", "SAE")),
    foi_lab = factor(foi,
                     levels = sort(unique(foi)),
                     labels = c("Small outbreak","Medium outbreak","Large outbreak"))
  )

ggplot(brr_long, aes(x = days, y = BRR, color = age_group)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(outcome ~ foi_lab, scales = "free_y") +
  labs(x = "Days of exposure", y = "Benefit–Risk Ratio (log scale)",
       color = "Age group") +
  my_theme
### endemic analysis------------------------------------------------------------
out_list_endemic <- list()

for (i in seq_len(nrow(all_risk))) {
  group <- all_risk[i, ]
  age <- group$age_group
  
  for (foi in foi_level_endemic) {
    for (days in T_days_vec_endemic) {
      
      AR <- compute_ar(foi, days = days)
      
      # SAE outcome
      sae <- compute_outcome(AR, group$p_sae_nat, group$p_sae_vacc, VE = ve)
      
      # Death outcome
      death <- compute_outcome(AR, group$p_death_nat, group$p_death_vacc, VE = ve)
      
      # Save (add 'days' column)
      out_list_endemic[[length(out_list_endemic) + 1]] <- data.frame(
        age_group         = age,
        foi               = foi,
        days              = days,
        AR                = AR,
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
# Combine results
result_df_endemic <- do.call(rbind, out_list_endemic)


## graphs
# ---- 1) Prepare long data for plotting (use your existing result_df) ----
# x = vaccine-attributable risk per 10k (excess), y = infection risk per 10k (unvaccinated)
death_df_endemic <- result_df_endemic %>%
  transmute(
    age_group, foi, days,
    outcome = "Death",
    x = excess_10k_death,
    y = risk_nv_10k_death
  )

sae_df_endemic <- result_df_endemic %>%
  transmute(
    age_group, foi, days,
    outcome = "SAE",
    x = excess_10k_sae,
    y = risk_nv_10k_sae
  )

plot_df_endemic <- bind_rows(death_df_endemic, sae_df_endemic)

# Ensure factor ordering if needed
plot_df_endemic$age_group <- factor(plot_df_endemic$age_group, levels = c("18-64", "65"))
plot_df_endemic$outcome   <- factor(plot_df_endemic$outcome,   levels = c("Death", "SAE"))

foi_levels <- sort(unique(plot_df_endemic$foi))
foi_labels <- setNames(c("Small outbreak","Medium outbreak","Large outbreak"),
                       as.character(foi_levels))
plot_df_endemic <- plot_df_endemic %>%
  mutate(foi_lab = factor(foi, 
                          levels = foi_levels,
                          labels = foi_labels))%>%
  mutate(days_lab = factor(paste0("Days: ", days)))

# ---- 2) Build threshold line y = x / VE across full x-range ----
# Use your VE defined earlier
x_max <- max(plot_df_endemic$x, na.rm = TRUE)
y_max <- max(plot_df_endemic$y, na.rm = TRUE)

bg_grid_endemic <- expand.grid(
  x = seq(0, x_max, length.out = 150),
  y = seq(0, y_max, length.out = 150),
  foi_lab = levels(plot_df_endemic$foi_lab),
  days    = sort(unique(plot_df_endemic$days))
) %>%
  mutate(net_per10k = ve * y - x)%>%
  mutate(days_lab = factor(paste0("Days: ", days)))

# threshold: benefit = risk
# benefit from vaccine = VE * y 
# risk from vaccine = x 
# VE * y = x
# y = x / VE 
thr_df_endemic <- data.frame(
  x = seq(0, max(plot_df_endemic$x, na.rm = TRUE) * 1.05, length.out = 400)
) %>%
  mutate(y = x / ve)

# ---- 3) Plot: FOI facets, threshold line, shapes by age, colors by outcome ----
net_rng  <- range(bg_grid_endemic$net_per10k, na.rm = TRUE)
max_abs  <- max(abs(net_rng))

ggplot()+
  geom_raster(data = bg_grid_endemic, aes(x = x, y = y, fill = net_per10k), alpha = 0.75)+
  geom_line(data = thr_df_endemic, aes(x = x, y = y),
            linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_point(data = plot_df_endemic,
             aes(x = x, y = y, shape = age_group, color = outcome),
             size = 3, stroke = 0.7)+
  facet_wrap(~foi_lab)+
  scale_shape_manual(values = c("18-64" = 22, "65" = 24)) + 
  scale_fill_gradientn(
    name   = "Outcomes averted\nper 10,000",
    colours = c("red", "white", "skyblue"),    
    values  = scales::rescale(c(-max_abs, 0, max_abs)),
    limits  = c(-max_abs, max_abs)
  ) +
  labs(
    x = "Vaccine-attributable severe outcome per 10,000",
    y = "Infection-attributable severe outcome per 10,000",
    shape = "Age group",
    color = "Outcome",
    caption = "Endemic situation assumed 365 days of exposure. Low endemic FOI = 0.01, Medium endemic FOI = 0.05, High endemic assumed FOI = 0.1."
  )+
  my_theme


