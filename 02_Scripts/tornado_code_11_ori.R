# ============================================================
# tornado_code_11_ori.R
#
# Spearman rank-correlation tornado plots for outbreak response
# immunisation (ORI). Mirrors the structure of tornado_code_11.R
# (traveller vaccination) but uses ORI draw-level BRR output.
#
# Required objects in the workspace before sourcing:
#   - draw_level_xy_serostatus_filtered (from brazil_travel_map.R
#     or brazil_all_draws_ori_v2.R; RR_seropos == 0 subset)
#   - lhs_sample (LHS parameter draws used downstream for DALY)
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(showtext)
library(sysfonts)

font_add("Calibri", regular    = "calibri.ttf",
                    bold       = "calibrib.ttf",
                    italic     = "calibrii.ttf",
                    bolditalic = "calibriz.ttf")
showtext_auto()
showtext_opts(dpi = 300)

# ------------------------------------------------------------
# 0. Parameters, labels, and ORI-specific choices
# ------------------------------------------------------------
# Pick which VE layer to profile. ORI has two VE levels:
#   "VE0"    = Disease blocking only
#   "VE98.9" = Disease and infection blocking
target_ve     <- "VE98.9"
target_ve_lab <- dplyr::case_when(
  target_ve == "VE0"    ~ "Disease blocking only",
  target_ve == "VE98.9" ~ "Disease and infection blocking",
  TRUE ~ target_ve
)

# BRR column to correlate against. Either "brr_base" (no
# serostatus adjustment) or "brr_adj" (serostatus adjusted).
target_brr_col <- "brr_base"

# Parameters of interest. Attack-rate, symptomatic-fraction, and
# per-infection SAE/death probabilities are not ORI-relevant LHS
# inputs (transmission is handled by the dynamic model + setting;
# disease harms come from hosp/CFR draws). VE is fixed at the VE
# column level, so the LHS ve draw is also excluded.
#
# Age-specific vaccine SAE / life-expectancy parameters are only
# included in the tornado for their own age scenario so that 65+
# inputs (e.g. p_death_vacc_65) are not shown in the 18-64 panel.
params_common <- c(
  "hosp", "fatal_hosp", "fatal_nonhosp",
  "dw_hosp", "dw_nonhosp", "dw_chronic",
  "dur_acute", "dur_subac", "dur_6m", "dur_12m", "dur_30m",
  "acute", "chr6m", "chr12m", "chr30m"
)

params_age_u65 <- c(
  "p_sae_vacc_u65", "p_death_vacc_u65",
  "le_lost_18_64"
)

params_age_65p <- c(
  "p_sae_vacc_65", "p_death_vacc_65",
  "le_lost_65"
)

get_params_for_age <- function(age) {
  age_key <- as.character(age)
  age_params <- switch(age_key,
    "18-64" = params_age_u65,
    "65+"   = params_age_65p,
    character(0)
  )
  c(age_params, params_common)
}

params_of_interest <- unique(c(params_common, params_age_u65, params_age_65p))

param_labels <- c(
  p_sae_vacc_u65   = "Vaccine SAE rate (18\u201364)",
  p_sae_vacc_65    = "Vaccine SAE rate (65+)",
  p_death_vacc_u65 = "Vaccine death rate (18\u201364)",
  p_death_vacc_65  = "Vaccine death rate (65+)",
  ve               = "Vaccine efficacy",
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
  le_lost_18_64    = "Life expectancy lost (18\u201364)",
  le_lost_65       = "Life expectancy lost (65+)"
)

# ------------------------------------------------------------
# 1. Merge LHS draws into ORI draw-level table by draw_id
# ------------------------------------------------------------
if (!exists("draw_level_xy_serostatus_filtered", inherits = TRUE)) {
  stop("draw_level_xy_serostatus_filtered not found. Run ",
       "brazil_all_draws_ori_v2.R and brazil_travel_map.R first.")
}
if (!exists("lhs_sample", inherits = TRUE)) {
  stop("lhs_sample not found. Source lhsm_param_dists.R / setup.R first.")
}

lhs_df         <- as.data.frame(lhs_sample)
lhs_df$draw_id <- seq_len(nrow(lhs_df))

ori_merged <- draw_level_xy_serostatus_filtered %>%
  dplyr::filter(
    outcome == "DALY",
    VE      == target_ve
  ) %>%
  dplyr::left_join(lhs_df, by = "draw_id")

# Normalise setting factor so Low / Moderate / High is consistent.
ori_merged <- ori_merged %>%
  dplyr::mutate(
    setting = factor(as.character(setting),
                     levels = c("Low", "Moderate", "High"))
  )

cat("ORI rows per setting:\n"); print(table(ori_merged$setting))
cat("ORI rows per AgeCat:\n"); print(table(ori_merged$AgeCat))

# ------------------------------------------------------------
# 2. Spearman correlation function (identical to travel version)
# ------------------------------------------------------------
compute_spearman <- function(data, brr_col, params) {
  data_clean <- data %>%
    dplyr::filter(is.finite(.data[[brr_col]]), !is.na(.data[[brr_col]]))

  results <- lapply(params, function(p) {
    if (!p %in% colnames(data_clean)) return(NULL)
    x <- data_clean[[p]]
    y <- data_clean[[brr_col]]
    v <- var(x, na.rm = TRUE)
    if (is.na(v) || v == 0) return(NULL)
    ct <- cor.test(x, y, method = "spearman", exact = FALSE)
    data.frame(param = p, rho = ct$estimate, pval = ct$p.value)
  })
  dplyr::bind_rows(results)
}

# ------------------------------------------------------------
# 3. Define scenarios: (age) x (transmission setting)
# ------------------------------------------------------------
scenarios <- list(
  list(label = "ORI BRR (DALY) | 18\u201364, low transmission",
       age = "18-64", setting = "Low",      brr_col = target_brr_col),
  list(label = "ORI BRR (DALY) | 65+, low transmission",
       age = "65+",   setting = "Low",      brr_col = target_brr_col),
  list(label = "ORI BRR (DALY) | 18\u201364, moderate transmission",
       age = "18-64", setting = "Moderate", brr_col = target_brr_col),
  list(label = "ORI BRR (DALY) | 65+, moderate transmission",
       age = "65+",   setting = "Moderate", brr_col = target_brr_col),
  list(label = "ORI BRR (DALY) | 18\u201364, high transmission",
       age = "18-64", setting = "High",     brr_col = target_brr_col),
  list(label = "ORI BRR (DALY) | 65+, high transmission",
       age = "65+",   setting = "High",     brr_col = target_brr_col)
)

# ------------------------------------------------------------
# 4. Run Spearman for all scenarios
# ------------------------------------------------------------
spearman_results <- lapply(scenarios, function(s) {
  sub <- ori_merged %>%
    dplyr::filter(as.character(AgeCat) == s$age,
                  as.character(setting) == s$setting)
  if (nrow(sub) == 0) {
    warning(paste("Zero rows for:", s$label)); return(NULL)
  }
  res <- compute_spearman(sub, s$brr_col, get_params_for_age(s$age))
  if (is.null(res) || nrow(res) == 0) {
    warning(paste("No Spearman results for:", s$label)); return(NULL)
  }
  res$scenario <- s$label
  res$brr_col  <- s$brr_col
  res
}) %>% dplyr::bind_rows()

cat("\nScenarios computed:", length(unique(spearman_results$scenario)), "\n")

# ------------------------------------------------------------
# 5. Per-scenario tornado plot function (same visual grammar
#    as the traveller version)
# ------------------------------------------------------------
TOP_N <- 12

plot_one_tornado <- function(scenario_label, top_n = TOP_N) {

  df <- spearman_results %>%
    dplyr::filter(scenario == scenario_label) %>%
    dplyr::mutate(
      param_label = dplyr::recode(param, !!!param_labels),
      abs_rho     = abs(rho),
      direction   = factor(
        ifelse(rho > 0, "Positive (favourable)", "Negative (unfavourable)"),
        levels = c("Positive (favourable)", "Negative (unfavourable)")
      )
    ) %>%
    dplyr::arrange(dplyr::desc(abs_rho)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::arrange(rho) %>%
    dplyr::mutate(param_label = factor(param_label, levels = param_label))

  ggplot(df, aes(x = rho, y = param_label, fill = direction)) +
    geom_col(width = 0.70, colour = NA) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey25") +
    geom_text(
      aes(
        x     = ifelse(rho < 0, rho - 0.015, rho + 0.015),
        label = sprintf("%.2f", rho),
        hjust = ifelse(rho < 0, 1, 0)
      ),
      size = 3, colour = "grey20"
    ) +
    scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(-0.75, 0.75, 0.25),
      labels = function(x) sprintf("%.2f", x)
    ) +
    scale_fill_manual(
      values = c(
        "Positive (favourable)"   = "#56B4E9",
        "Negative (unfavourable)" = "#E69F00"
      ),
      drop = FALSE
    ) +
    labs(
      title = scenario_label,
      x     = NULL,
      y     = NULL,
      fill  = NULL
    ) +
    theme_minimal(base_size = 6, base_family = "Calibri") +
    theme(
      plot.background    = element_rect(fill = "white", colour = "grey72",
                                        linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
      axis.text.y        = element_text(size = 10, hjust = 1),
      axis.text.x        = element_text(size = 10),
      legend.position    = "none"
    ) +
    theme(
      plot.title = element_text(
        size   = 13,
        face   = "bold",
        hjust  = 0.5,
        margin = margin(b = 4)
      )
    )
}

# ------------------------------------------------------------
# 6. Build all 6 panels
# ------------------------------------------------------------
scenario_labels      <- sapply(scenarios, `[[`, "label")
tornado_plots        <- lapply(scenario_labels, plot_one_tornado)
names(tornado_plots) <- scenario_labels

# ------------------------------------------------------------
# 7. Shared legend
# ------------------------------------------------------------
legend_plot <- ggplot(
  data.frame(
    x = c(0.1, -0.1),
    y = c("a", "b"),
    direction = factor(c("Positive (favourable)", "Negative (unfavourable)"),
                       levels = c("Positive (favourable)",
                                  "Negative (unfavourable)"))
  ),
  aes(x = x, y = y, fill = direction)
) +
  geom_col() +
  scale_fill_manual(
    values = c("Positive (favourable)"   = "#56B4E9",
               "Negative (unfavourable)" = "#E69F00")
  ) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 11, family = "Calibri"),
    legend.key.size  = unit(0.5, "cm"),
    legend.title     = element_blank()
  )

shared_legend <- cowplot::get_legend(legend_plot)

# ------------------------------------------------------------
# 8. Assemble 6 panels in 2 columns
# ------------------------------------------------------------
combined <- wrap_plots(tornado_plots, ncol = 2) +
  plot_annotation(
    title    = paste0(
      "Spearman rank correlations between input parameters and the benefit-risk ratio (DALY) ",
      "of outbreak response immunisation"
    ),
    subtitle = paste0(
      "Top ", TOP_N, " parameters by |\u03c1| per scenario (low \u2192 moderate \u2192 high transmission).  ",
      "VE layer: ", target_ve_lab, "; BRR type: ",
      ifelse(target_brr_col == "brr_base", "base (no serostatus adjustment)", "serostatus-adjusted"), ".\n",
      "Positive \u03c1 \u2192 BRR increases (vaccination more favourable).  ",
      "Negative \u03c1 \u2192 BRR decreases (vaccination less favourable)."
    ),
    caption  = paste0(
      "Transmission categories are defined by state-level setting assignment (Low / Moderate / High).\n",
      "BRR (DALY) = DALY averted / DALY from vaccine-associated SAE per 10,000.  BRR >1 = net benefit."
    ),
    theme = theme(
      plot.title    = element_text(size = 16, face = "bold", family = "Calibri",
                                   margin = margin(t = 2, b = 4)),
      plot.subtitle = element_text(size = 12, colour = "grey20", family = "Calibri",
                                   margin = margin(t = 2, b = 10), lineheight = 1.25),
      plot.caption  = element_text(size = 12, colour = "grey20", family = "Calibri",
                                   hjust = 0, margin = margin(t = 10, b = 4),
                                   lineheight = 1.25),
      plot.margin   = margin(t = 10, r = 10, b = 10, l = 10)
    )
  )

print(combined)

# ------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------
output_dir <- "03_Manuscript/Submission/second_submission_lancet/suppl_figs/fig2_spearman"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(output_dir, paste0("tornado_spearman_ori_", target_ve, "_", target_brr_col, ".pdf")),
  combined,
  width = 12, height = 9.5, device = cairo_pdf
)
