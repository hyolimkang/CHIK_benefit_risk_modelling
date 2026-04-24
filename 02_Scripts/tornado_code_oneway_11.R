## ============================================================
## Deterministic One-way Sensitivity Tornado (BRR, DALY)
## - Center point: BRR using median of all parameters
## - Side points : BRR when one parameter is changed by +/-10%
## ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(showtext)
library(sysfonts)

# Optional font setup (safe if files exist in working directory)
font_add("Calibri", regular = "calibri.ttf", bold = "calibrib.ttf")
showtext_auto()
showtext_opts(dpi = 300)

# ------------------------------------------------------------
# 0. Inputs expected in environment
# ------------------------------------------------------------
required_objects <- c(
  "lhs_sample", "foi_daily_by_state", "states_to_run",
  "compute_outcome", "compute_daly_one", "compute_ar_travel"
)
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", "),
    ". Run setup script first (e.g., source('02_Scripts/setup.R'))."
  )
}

compute_outcome_fn <- get("compute_outcome", mode = "function")
compute_daly_one_fn <- get("compute_daly_one", mode = "function")
compute_ar_travel_fn <- get("compute_ar_travel", mode = "function")
utils::globalVariables(c(
  "scenario", "param_label", "low", "high",
  "brr_minus10", "brr_base", "brr_plus10"
))

params_of_interest <- c(
  "p_sae_vacc_u65", "p_sae_vacc_65",
  "p_death_vacc_u65", "p_death_vacc_65",
  "p_sae_nat_64", "p_sae_nat_65",
  "p_death_nat_64", "p_death_nat_65",
  "ve",
  "ar_small", "ar_med", "ar_large",
  "hosp", "fatal_hosp", "fatal_nonhosp",
  "dw_hosp", "dw_nonhosp", "dw_chronic",
  "dur_acute", "dur_subac", "dur_6m", "dur_12m", "dur_30m",
  "acute", "chr6m", "chr12m", "chr30m",
  "le_lost_18_64", "le_lost_65",
  "symp_overall"
)

param_labels <- c(
  p_sae_vacc_u65   = "Vaccine SAE rate (18-64)",
  p_sae_vacc_65    = "Vaccine SAE rate (65+)",
  p_death_vacc_u65 = "Vaccine death rate (18-64)",
  p_death_vacc_65  = "Vaccine death rate (65+)",
  p_sae_nat_64     = "Natural infection SAE rate (18-64)",
  p_sae_nat_65     = "Natural infection SAE rate (65+)",
  p_death_nat_64   = "Natural infection death rate (18-64)",
  p_death_nat_65   = "Natural infection death rate (65+)",
  ve               = "Vaccine efficacy",
  ar_small         = "Attack rate (low transmission)",
  ar_med           = "Attack rate (moderate transmission)",
  ar_large         = "Attack rate (high transmission)",
  hosp             = "Hospitalisation rate given infection",
  fatal_hosp       = "CFR (hospitalised)",
  fatal_nonhosp    = "CFR (non-hospitalised)",
  dw_hosp          = "Disability weight (hospitalised acute)",
  dw_nonhosp       = "Disability weight (non-hospitalised acute)",
  dw_chronic       = "Disability weight (chronic)",
  dur_acute        = "Duration (acute)",
  dur_subac        = "Duration (subacute)",
  dur_6m           = "Duration (chronic 6 months)",
  dur_12m          = "Duration (chronic 12 months)",
  dur_30m          = "Duration (chronic 30 months)",
  acute            = "Proportion acute only",
  chr6m            = "Proportion chronic 6m",
  chr12m           = "Proportion chronic 12m",
  chr30m           = "Proportion chronic 30m",
  le_lost_18_64    = "Life expectancy lost (18-64)",
  le_lost_65       = "Life expectancy lost (65+)",
  symp_overall     = "Symptomatic fraction"
)

# ------------------------------------------------------------
# 1. Build median baseline parameter vector
# ------------------------------------------------------------
lhs_df <- as.data.frame(lhs_sample)
median_pars <- vapply(
  names(lhs_df),
  function(p) median(lhs_df[[p]], na.rm = TRUE),
  numeric(1)
)
median_pars <- as.list(median_pars)

bound_prob <- function(x) pmin(pmax(x, 0), 1)

vary_param_10pct <- function(pars, param, multiplier) {
  out <- pars
  val <- as.numeric(out[[param]])
  if (!is.finite(val)) return(out)

  out[[param]] <- val * multiplier

  # Keep probabilities within [0, 1]
  prob_like <- grepl("^(p_|ve$|symp_|hosp$|fatal_|dw_|acute$|subac$|chr)", param)
  if (prob_like) out[[param]] <- bound_prob(out[[param]])

  out
}

scenarios <- list(
  list(label = "OWSA BRR (DALY) | 18-64, low transmission", age = "18-64", ar = "low"),
  list(label = "OWSA BRR (DALY) | 65+, low transmission", age = "65+", ar = "low"),
  list(label = "OWSA BRR (DALY) | 18-64, moderate transmission", age = "18-64", ar = "moderate"),
  list(label = "OWSA BRR (DALY) | 65+, moderate transmission", age = "65+", ar = "moderate"),
  list(label = "OWSA BRR (DALY) | 18-64, high transmission", age = "18-64", ar = "high"),
  list(label = "OWSA BRR (DALY) | 65+, high transmission", age = "65+", ar = "high")
)

# ------------------------------------------------------------
# 2. Run one-way sensitivity using PSA-like travel engine
# ------------------------------------------------------------
region_key <- c(
  "Ceará" = "ce", "Bahia" = "bh", "Paraíba" = "pa", "Pernambuco" = "pn",
  "Rio Grande do Norte" = "rg", "Piauí" = "pi", "Tocantins" = "tc",
  "Alagoas" = "ag", "Minas Gerais" = "mg", "Sergipe" = "se", "Goiás" = "go"
)

OWSA_ENTRY_SAMPLES <- 10
OWSA_AGES <- c("18-64", "65+")

run_engine_like_psa <- function(pars, entry_samples = OWSA_ENTRY_SAMPLES, seed = 12345) {
  set.seed(seed)
  states_local <- get("states_to_run", inherits = TRUE)
  foi_by_state <- get("foi_daily_by_state", inherits = TRUE)

  h0_by_state <- lapply(states_local, function(st) {
    foi0 <- foi_by_state[[st]]
    if (is.null(foi0)) stop("State not found in foi_daily_by_state: ", st)
    sum(foi0, na.rm = TRUE)
  })
  names(h0_by_state) <- states_local

  travel_days <- list(
    "7d" = as.numeric(pars$trav_7d),
    "14d" = as.numeric(pars$trav_14d),
    "30d" = as.numeric(pars$trav_30d),
    "90d" = as.numeric(pars$trav_90d)
  )
  if (any(!is.finite(unlist(travel_days)))) {
    stop("Travel duration parameters (trav_7d/trav_14d/trav_30d/trav_90d) are missing.")
  }

  ve_d <- as.numeric(pars$ve)
  symp_prop_d <- as.numeric(pars$symp_overall)
  out <- list()
  idx <- 1L

  for (st in states_local) {
    foi0 <- foi_by_state[[st]]
    key <- region_key[[st]]
    if (is.null(key)) stop("No region_key mapping for state: ", st)

    ar_col <- paste0("ar_", key)
    if (!ar_col %in% names(pars)) stop("Missing state AR parameter: ", ar_col)

    ar_target_raw <- as.numeric(pars[[ar_col]])
    eps_ar <- 1e-12
    ar_target <- pmin(pmax(ar_target_raw, eps_ar), 1 - eps_ar)

    h0 <- h0_by_state[[st]]
    if (!is.finite(h0) || h0 <= 0) stop("Non-positive FOI sum for state: ", st)

    ht <- -log(1 - ar_target)
    m <- ht / h0
    foi_daily_scaled <- foi0 * m
    ar_total_state <- 1 - exp(-sum(foi_daily_scaled, na.rm = TRUE))

    for (age in OWSA_AGES) {
      if (age == "65+") {
        p_sae_vacc <- as.numeric(pars$p_sae_vacc_65)
        p_death_vacc <- as.numeric(pars$p_death_vacc_65)
        p_sae_nat <- as.numeric(pars$p_sae_nat_65)
        p_death_nat <- as.numeric(pars$p_death_nat_65)
      } else {
        p_sae_vacc <- as.numeric(pars$p_sae_vacc_u65)
        p_death_vacc <- as.numeric(pars$p_death_vacc_u65)
        p_sae_nat <- as.numeric(pars$p_sae_nat_64)
        p_death_nat <- as.numeric(pars$p_death_nat_64)
      }

      for (days_label in names(travel_days)) {
        d_stay <- max(1L, round(as.numeric(travel_days[[days_label]])))
        l_days <- length(foi0)
        max_entry <- max(1L, l_days - d_stay + 1L)
        entry_days <- sample.int(max_entry, size = entry_samples, replace = TRUE)

        brr_vals <- vapply(entry_days, function(entry_day) {
          ar_travel <- compute_ar_travel_fn(foi_daily_scaled, entry_day, d_stay)

          res <- compute_outcome_fn(
            AR = ar_travel,
            p_hosp = symp_prop_d * p_sae_nat,
            p_death = symp_prop_d * p_death_nat,
            p_sae_vacc = p_sae_vacc,
            p_death_vacc = p_death_vacc,
            VE_hosp = ve_d,
            VE_death = ve_d
          )

          symp_nv_10k <- 1e4 * ar_travel * symp_prop_d
          symp_v_10k <- symp_nv_10k * (1 - ve_d)
          nonhosp_symp_nv_10k <- pmax(0, symp_nv_10k - res$risk_nv_hosp)
          nonhosp_symp_v_10k <- pmax(0, symp_v_10k - res$risk_v_hosp)

          dz_nv <- compute_daly_one_fn(
            age_group = age,
            deaths_10k = res$risk_nv_death,
            hosp_10k = res$risk_nv_hosp,
            nonhosp_symp_10k = nonhosp_symp_nv_10k,
            symp_10k = symp_nv_10k,
            sae_10k = 0, deaths_sae_10k = 0,
            draw_pars = pars
          )
          dz_v <- compute_daly_one_fn(
            age_group = age,
            deaths_10k = res$risk_v_death,
            hosp_10k = res$risk_v_hosp,
            nonhosp_symp_10k = nonhosp_symp_v_10k,
            symp_10k = symp_v_10k,
            sae_10k = 0, deaths_sae_10k = 0,
            draw_pars = pars
          )
          sae <- compute_daly_one_fn(
            age_group = age,
            deaths_10k = 0, hosp_10k = 0, nonhosp_symp_10k = 0, symp_10k = 0,
            sae_10k = 1e4 * p_sae_vacc,
            deaths_sae_10k = 1e4 * p_death_vacc,
            draw_pars = pars
          )

          daly_averted <- dz_nv$daly_dz - dz_v$daly_dz
          daly_sae <- sae$daly_sae
          ifelse(daly_sae > 0, daly_averted / daly_sae, NA_real_)
        }, numeric(1))

        out[[idx]] <- data.frame(
          state = st,
          age_group = age,
          days = days_label,
          AR_total = ar_total_state,
          brr_daly = median(brr_vals, na.rm = TRUE)
        )
        idx <- idx + 1L
      }
    }
  }

  bind_rows(out) %>%
    mutate(
      ar_cat = case_when(
        AR_total < 0.05 ~ "low",
        AR_total < 0.10 ~ "moderate",
        AR_total >= 0.10 ~ "high"
      ),
      ar_cat = factor(.data$ar_cat, levels = c("low", "moderate", "high"))
    )
}

scenario_brr <- function(engine_df, age, ar_cat_level) {
  sub <- engine_df %>% filter(.data$age_group == age, .data$ar_cat == ar_cat_level)
  if (nrow(sub) == 0) return(NA_real_)
  median(sub$brr_daly, na.rm = TRUE)
}

calc_scenario_vector <- function(pars, scenarios_input, seed = 12345) {
  engine_df <- run_engine_like_psa(pars, entry_samples = OWSA_ENTRY_SAMPLES, seed = seed)
  values <- vapply(
    scenarios_input,
    function(s) scenario_brr(engine_df, s$age, s$ar),
    numeric(1)
  )
  names(values) <- vapply(scenarios_input, `[[`, character(1), "label")
  values
}

target_params <- intersect(params_of_interest, names(median_pars))
base_vec <- calc_scenario_vector(median_pars, scenarios, seed = 4001)

rows <- lapply(seq_along(target_params), function(i) {
  p <- target_params[[i]]
  cat("OWSA parameter", i, "of", length(target_params), ":", p, "\n")

  pars_minus <- vary_param_10pct(median_pars, p, 0.90)
  pars_plus <- vary_param_10pct(median_pars, p, 1.10)

  minus_vec <- calc_scenario_vector(pars_minus, scenarios, seed = 4001)
  plus_vec <- calc_scenario_vector(pars_plus, scenarios, seed = 4001)

  data.frame(
    scenario = names(base_vec),
    param = p,
    brr_base = as.numeric(base_vec),
    brr_minus10 = as.numeric(minus_vec[names(base_vec)]),
    brr_plus10 = as.numeric(plus_vec[names(base_vec)])
  )
})

owsa_results <- bind_rows(rows) %>%
  mutate(
    param = as.character(param),
    param_label = dplyr::coalesce(unname(param_labels[param]), param),
    low = pmin(brr_minus10, brr_plus10, na.rm = TRUE),
    high = pmax(brr_minus10, brr_plus10, na.rm = TRUE),
    range = high - low
  ) %>%
  filter(is.finite(brr_base), is.finite(low), is.finite(high))

TOP_N <- 12
SHARE_X_ACROSS_PANELS <- TRUE

# Global x-axis range for direct panel-to-panel comparison
x_all <- c(
  owsa_results$low,
  owsa_results$high,
  owsa_results$brr_base,
  owsa_results$brr_minus10,
  owsa_results$brr_plus10
)
x_all <- x_all[is.finite(x_all)]
if (length(x_all) == 0) stop("No finite BRR values available for plotting.")

x_min <- min(x_all, na.rm = TRUE)
x_max <- max(x_all, na.rm = TRUE)
x_lower <- min(0, 1, x_min)
x_upper <- max(1, x_max)
x_pad <- 0.03 * (x_upper - x_lower)
x_limits_shared <- c(x_lower - x_pad, x_upper + x_pad)

plot_one_owsa <- function(scenario_label, top_n = TOP_N) {
  x_limits_local <- get("x_limits_shared", inherits = TRUE)

  df <- owsa_results %>%
    filter(.data$scenario == scenario_label) %>%
    arrange(desc(.data$range)) %>%
    slice_head(n = top_n) %>%
    arrange(.data$range) %>%
    mutate(param_label = factor(.data$param_label, levels = .data$param_label))

  ggplot(df, aes(y = .data$param_label)) +
    geom_segment(
      aes(x = .data$low, xend = .data$high, yend = .data$param_label),
      linewidth = 4, colour = "grey85"
    ) +
    geom_point(aes(x = .data$brr_minus10, colour = "-10%"), size = 2) +
    geom_point(aes(x = .data$brr_base, colour = "Baseline median"), size = 2.2) +
    geom_point(aes(x = .data$brr_plus10, colour = "+10%"), size = 2) +
    geom_vline(xintercept = 1, linewidth = 0.35, colour = "grey25", linetype = "dashed") +
    {
      if (SHARE_X_ACROSS_PANELS) {
        scale_x_continuous(
          limits = x_limits_local,
          labels = scales::label_number(accuracy = 0.1, big.mark = ",")
        )
      } else {
        scale_x_continuous(
          labels = scales::label_number(accuracy = 0.1, big.mark = ",")
        )
      }
    } +
    scale_colour_manual(
      values = c(
        "-10%" = "#E69F00",
        "Baseline median" = "#222222",
        "+10%" = "#56B4E9"
      ),
      breaks = c("-10%", "Baseline median", "+10%")
    ) +
    labs(
      title = scenario_label,
      x = "Benefit-risk ratio (BRR, DALY)",
      y = NULL,
      colour = NULL
    ) +
    theme_minimal(base_size = 6, base_family = "Calibri") +
    theme(
      plot.background = element_rect(fill = "white", colour = "grey72", linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
      axis.text.y = element_text(size = 10, hjust = 1),
      axis.text.x = element_text(size = 10),
      legend.position = "none"
    ) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 4))
    )
}

scenario_labels <- vapply(scenarios, `[[`, character(1), "label")
owsa_plots <- lapply(scenario_labels, plot_one_owsa)
names(owsa_plots) <- scenario_labels

# Shared legend
legend_plot <- ggplot(
  data.frame(
    x = c(1, 2, 3),
    y = c(1, 1, 1),
    change = factor(c("-10%", "Baseline median", "+10%"),
                    levels = c("-10%", "Baseline median", "+10%"))
  ),
  aes(x = x, y = y, colour = change)
) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c(
      "-10%" = "#E69F00",
      "Baseline median" = "#222222",
      "+10%" = "#56B4E9"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 11, family = "Calibri"),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_blank()
  )

shared_legend <- cowplot::get_legend(legend_plot)

combined <- wrap_plots(owsa_plots, ncol = 2) +
  plot_annotation(
    title = "One-way sensitivity analysis for BRR (DALY)",
    subtitle = paste0(
      "Center point uses all-parameter medians. ",
      "Each row shows BRR when one parameter is changed by +/-10%."
    ),
    caption = "Dashed vertical line at BRR = 1 (net-benefit threshold).",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", family = "Calibri"),
      plot.subtitle = element_text(size = 12, colour = "grey20", family = "Calibri", margin = margin(b = 6)),
      plot.caption = element_text(size = 12, colour = "grey20", family = "Calibri", hjust = 0)
    )
  )

final_plot <- wrap_plots(
  combined,
  patchwork::wrap_elements(shared_legend),
  ncol = 1,
  heights = c(18, 1.2)
)

print(final_plot)

output_dir <- "03_Manuscript/Submission/second_submission_lancet/suppl_figs/fig3_owsa"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(output_dir, "owsa.pdf"),
  combined,
  width = 12,
  height = 8,
  device = cairo_pdf
)
