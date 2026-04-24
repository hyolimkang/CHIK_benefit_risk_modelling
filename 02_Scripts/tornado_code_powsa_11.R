## ============================================================
## Probabilistic One-way Sensitivity Tornado (p-OWSA, BRR DALY)
## - Baseline: PSA draws (no parameter perturbation)
## - One-way: target parameter only +/-10%, others remain probabilistic
## - Engine: state FOI scaling + entry-day sampling (PSA-like)
## ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(showtext)
library(sysfonts)

# Optional font setup
if (file.exists("calibri.ttf") && file.exists("calibrib.ttf")) {
  font_add("Calibri", regular = "calibri.ttf", bold = "calibrib.ttf")
}
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

utils::globalVariables(c(
  "scenario", "param_label", "range_med",
  "ui_low", "ui_high", "low_med", "high_med",
  "brr_minus_med", "brr_base_med", "brr_plus_med"
))

bound_prob <- function(x) pmin(pmax(x, 0), 1)

vary_param_10pct <- function(pars, param, multiplier) {
  out <- pars
  if (!param %in% names(out)) return(out)
  val <- suppressWarnings(as.numeric(out[[param]]))
  if (!is.finite(val)) return(out)

  out[[param]] <- val * multiplier

  # Keep probabilistic quantities in [0, 1]
  prob_like <- grepl("^(p_|ve$|symp_|hosp$|fatal_|dw_|acute$|subac$|chr)", param)
  if (prob_like) out[[param]] <- bound_prob(as.numeric(out[[param]]))

  out
}

region_key <- c(
  "Ceará" = "ce", "Bahia" = "bh", "Paraíba" = "pa", "Pernambuco" = "pn",
  "Rio Grande do Norte" = "rg", "Piauí" = "pi", "Tocantins" = "tc",
  "Alagoas" = "ag", "Minas Gerais" = "mg", "Sergipe" = "se", "Goiás" = "go"
)

scenarios <- list(
  list(label = "p-OWSA BRR (DALY) | 18-64, low transmission", age = "18-64", ar = "low"),
  list(label = "p-OWSA BRR (DALY) | 65+, low transmission", age = "65+", ar = "low"),
  list(label = "p-OWSA BRR (DALY) | 18-64, moderate transmission", age = "18-64", ar = "moderate"),
  list(label = "p-OWSA BRR (DALY) | 65+, moderate transmission", age = "65+", ar = "moderate"),
  list(label = "p-OWSA BRR (DALY) | 18-64, high transmission", age = "18-64", ar = "high"),
  list(label = "p-OWSA BRR (DALY) | 65+, high transmission", age = "65+", ar = "high")
)

# ------------------------------------------------------------
# 1. p-OWSA runtime controls
# ------------------------------------------------------------
# For strict/full p-OWSA, set POWSA_N_DRAWS <- nrow(lhs_sample)
POWSA_N_DRAWS <- min(100L, nrow(lhs_sample))
POWSA_ENTRY_SAMPLES <- 20L
POWSA_SEED <- 2026
TOP_N <- 12
SHARE_X_ACROSS_PANELS <- FALSE

draw_ids <- if (POWSA_N_DRAWS < nrow(lhs_sample)) {
  set.seed(POWSA_SEED)
  sort(sample.int(nrow(lhs_sample), size = POWSA_N_DRAWS, replace = FALSE))
} else {
  seq_len(nrow(lhs_sample))
}

states_local <- get("states_to_run", inherits = TRUE)
foi_by_state <- get("foi_daily_by_state", inherits = TRUE)
lhs_sample_local <- get("lhs_sample", inherits = TRUE)

h0_by_state <- lapply(states_local, function(st) {
  foi0 <- foi_by_state[[st]]
  if (is.null(foi0)) stop("State not found in foi_daily_by_state: ", st)
  sum(foi0, na.rm = TRUE)
})
names(h0_by_state) <- states_local

scenario_from_engine <- function(engine_df, age, ar_cat_level) {
  sub <- engine_df %>% filter(.data$age_group == age, .data$ar_cat == ar_cat_level)
  if (nrow(sub) == 0) return(NA_real_)
  median(sub$brr_daly, na.rm = TRUE)
}

run_engine_one_draw <- function(pars, entry_samples = POWSA_ENTRY_SAMPLES) {
  travel_days <- list(
    "7d" = as.numeric(pars$trav_7d),
    "14d" = as.numeric(pars$trav_14d),
    "30d" = as.numeric(pars$trav_30d),
    "90d" = as.numeric(pars$trav_90d)
  )
  if (any(!is.finite(unlist(travel_days)))) {
    stop("Travel duration parameters are missing for this draw.")
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

    for (age in c("18-64", "65+")) {
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

calc_scenario_vector_for_draw <- function(draw_idx, param_name = NULL, multiplier = 1) {
  pars <- as.list(lhs_sample_local[draw_idx, , drop = FALSE])
  if (!is.null(param_name)) {
    pars <- vary_param_10pct(pars, param_name, multiplier)
  }

  engine_df <- run_engine_one_draw(pars, entry_samples = POWSA_ENTRY_SAMPLES)
  vals <- vapply(scenarios, function(s) {
    scenario_from_engine(engine_df, s$age, s$ar)
  }, numeric(1))
  names(vals) <- vapply(scenarios, `[[`, character(1), "label")
  vals
}

run_psa_for_case <- function(param_name = NULL, multiplier = 1) {
  out <- matrix(
    NA_real_,
    nrow = length(draw_ids),
    ncol = length(scenarios),
    dimnames = list(NULL, vapply(scenarios, `[[`, character(1), "label"))
  )

  for (i in seq_along(draw_ids)) {
    if (i %% 10 == 0) {
      cat("p-OWSA case", ifelse(is.null(param_name), "baseline", param_name),
          "draw", i, "of", length(draw_ids), "\n")
    }
    out[i, ] <- calc_scenario_vector_for_draw(
      draw_idx = draw_ids[[i]],
      param_name = param_name,
      multiplier = multiplier
    )
  }
  out
}

summarise_case <- function(mat_case) {
  data.frame(
    scenario = colnames(mat_case),
    med = apply(mat_case, 2, median, na.rm = TRUE),
    lo = apply(mat_case, 2, quantile, probs = 0.025, na.rm = TRUE),
    hi = apply(mat_case, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

# ------------------------------------------------------------
# 2. Baseline and +/-10% for each parameter
# ------------------------------------------------------------
baseline_mat <- run_psa_for_case(param_name = NULL, multiplier = 1)
baseline_sum <- summarise_case(baseline_mat) %>%
  rename(brr_base_med = med, brr_base_lo = lo, brr_base_hi = hi)

target_params <- intersect(params_of_interest, names(lhs_sample_local))
rows <- lapply(seq_along(target_params), function(i) {
  p <- target_params[[i]]
  cat("Running p-OWSA parameter", i, "of", length(target_params), ":", p, "\n")

  minus_mat <- run_psa_for_case(param_name = p, multiplier = 0.90)
  plus_mat <- run_psa_for_case(param_name = p, multiplier = 1.10)

  minus_sum <- summarise_case(minus_mat) %>%
    rename(brr_minus_med = med, brr_minus_lo = lo, brr_minus_hi = hi)
  plus_sum <- summarise_case(plus_mat) %>%
    rename(brr_plus_med = med, brr_plus_lo = lo, brr_plus_hi = hi)

  baseline_sum %>%
    left_join(minus_sum, by = "scenario") %>%
    left_join(plus_sum, by = "scenario") %>%
    mutate(param = p)
})

powsa_results <- bind_rows(rows) %>%
  mutate(
    param = as.character(.data$param),
    param_label = dplyr::coalesce(unname(param_labels[.data$param]), .data$param),
    low_med = pmin(.data$brr_minus_med, .data$brr_plus_med, na.rm = TRUE),
    high_med = pmax(.data$brr_minus_med, .data$brr_plus_med, na.rm = TRUE),
    ui_low = pmin(.data$brr_minus_lo, .data$brr_plus_lo, na.rm = TRUE),
    ui_high = pmax(.data$brr_minus_hi, .data$brr_plus_hi, na.rm = TRUE),
    range_med = .data$high_med - .data$low_med
  ) %>%
  filter(
    is.finite(.data$brr_base_med),
    is.finite(.data$low_med),
    is.finite(.data$high_med),
    is.finite(.data$ui_low),
    is.finite(.data$ui_high)
  )

# ------------------------------------------------------------
# 3. Plot
# ------------------------------------------------------------
x_all <- c(
  powsa_results$ui_low, powsa_results$ui_high,
  powsa_results$brr_minus_med, powsa_results$brr_base_med, powsa_results$brr_plus_med
)
x_all <- x_all[is.finite(x_all)]
if (length(x_all) == 0) stop("No finite BRR values available for plotting.")
x_min <- min(x_all, na.rm = TRUE)
x_max <- max(x_all, na.rm = TRUE)
x_lower <- min(0, 1, x_min)
x_upper <- max(1, x_max)
x_pad <- 0.03 * (x_upper - x_lower)
x_limits_shared <- c(x_lower - x_pad, x_upper + x_pad)

plot_one_powsa <- function(scenario_label, top_n = TOP_N) {
  df <- powsa_results %>%
    filter(.data$scenario == scenario_label) %>%
    arrange(desc(.data$range_med)) %>%
    slice_head(n = top_n) %>%
    arrange(.data$range_med) %>%
    mutate(param_label = factor(.data$param_label, levels = .data$param_label))

  ggplot(df, aes(y = .data$param_label)) +
    # outer uncertainty envelope across +/-10% (95% UI)
    geom_segment(
      aes(x = .data$ui_low, xend = .data$ui_high, yend = .data$param_label),
      linewidth = 2.2, colour = "grey85"
    ) +
    # median impact range between -10% and +10%
    geom_segment(
      aes(x = .data$low_med, xend = .data$high_med, yend = .data$param_label),
      linewidth = 4.3, colour = "grey70"
    ) +
    geom_point(aes(x = .data$brr_minus_med, colour = "-10% median"), size = 2) +
    geom_point(aes(x = .data$brr_base_med, colour = "Baseline median"), size = 2.2) +
    geom_point(aes(x = .data$brr_plus_med, colour = "+10% median"), size = 2) +
    geom_vline(xintercept = 1, linewidth = 0.35, colour = "grey25", linetype = "dashed") +
    scale_x_continuous(
      labels = scales::label_number(accuracy = 0.1, big.mark = ",")
    ) +
    scale_colour_manual(
      values = c(
        "-10% median" = "#E69F00",
        "Baseline median" = "#222222",
        "+10% median" = "#56B4E9"
      ),
      breaks = c("-10% median", "Baseline median", "+10% median")
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
powsa_plots <- lapply(scenario_labels, plot_one_powsa)
names(powsa_plots) <- scenario_labels

legend_plot <- ggplot(
  data.frame(
    x = c(1, 2, 3),
    y = c(1, 1, 1),
    change = factor(c("-10% median", "Baseline median", "+10% median"),
                    levels = c("-10% median", "Baseline median", "+10% median"))
  ),
  aes(x = x, y = y, colour = change)
) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c(
      "-10% median" = "#E69F00",
      "Baseline median" = "#222222",
      "+10% median" = "#56B4E9"
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

combined <- wrap_plots(powsa_plots, ncol = 2) +
  plot_annotation(
    title = "Probabilistic one-way sensitivity analysis for BRR (DALY)",
    subtitle = paste0(
      "Baseline and +/-10% effects use PSA draws (n=", length(draw_ids),
      "). Grey thin bar: pooled 95% UI across +/-10%; grey thick bar: median impact."
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

# Save figure and raw table
output_dir <- "03_Manuscript/Submission/second_submission_lancet/suppl_figs/fig3_powsa"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(output_dir, "powsa.pdf"), combined,
       width = 12, height = 8, device = cairo_pdf)

write.csv(
  powsa_results,
  file = file.path(output_dir, "powsa_results.csv"),
  row.names = FALSE
)
