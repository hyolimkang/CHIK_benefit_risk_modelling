# ============================================================
# PSA: Spearman Tornado + P(BRR > 1)
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(showtext)
library(sysfonts)

# Load Calibri font
# Option 1: if Calibri is installed on the system (Windows/Mac with Office)
font_add("Calibri", regular    = "calibri.ttf",
                    bold       = "calibrib.ttf",
                    italic     = "calibrii.ttf",
                    bolditalic = "calibriz.ttf")
showtext_auto()
showtext_opts(dpi = 300)

# ------------------------------------------------------------
# 0. Parameters and labels
# ------------------------------------------------------------
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
  p_sae_vacc_u65   = "Vaccine SAE rate (18\u201364)",
  p_sae_vacc_65    = "Vaccine SAE rate (65+)",
  p_death_vacc_u65 = "Vaccine death rate (18\u201364)",
  p_death_vacc_65  = "Vaccine death rate (65+)",
  p_sae_nat_64     = "Natural infection SAE rate (18\u201364)",
  p_sae_nat_65     = "Natural infection SAE rate (65+)",
  p_death_nat_64   = "Natural infection death rate (18\u201364)",
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
  le_lost_18_64    = "Life expectancy lost (18\u201364)",
  le_lost_65       = "Life expectancy lost (65+)",
  symp_overall     = "Symptomatic fraction"
)

# ------------------------------------------------------------
# 1. Merge LHS draws into psa_df by draw index
# ------------------------------------------------------------
lhs_df      <- as.data.frame(lhs_sample)
lhs_df$draw <- seq_len(nrow(lhs_df))

psa_merged <- psa_df %>%
  left_join(lhs_df, by = "draw")

# ------------------------------------------------------------
# 2. AR categorisation: manuscript-defined cutoffs
# ------------------------------------------------------------
psa_merged <- psa_merged %>%
  mutate(
    ar_cat = case_when(
      AR_total <  0.05 ~ "low",
      AR_total <  0.10 ~ "moderate",
      AR_total >= 0.10 ~ "high"
    ),
    ar_cat = factor(ar_cat, levels = c("low", "moderate", "high"))
  )

cat("AR cutoffs: low <5%, moderate 5-10%, high >=10%\n")
cat("Rows per ar_cat:\n"); print(table(psa_merged$ar_cat))

# ------------------------------------------------------------
# 3. Spearman correlation function
# ------------------------------------------------------------
compute_spearman <- function(data, brr_col, params) {
  data_clean <- data %>%
    filter(is.finite(.data[[brr_col]]), !is.na(.data[[brr_col]]))

  results <- lapply(params, function(p) {
    if (!p %in% colnames(data_clean)) return(NULL)
    x <- data_clean[[p]]
    y <- data_clean[[brr_col]]
    v <- var(x, na.rm = TRUE)
    if (is.na(v) || v == 0) return(NULL)
    ct <- cor.test(x, y, method = "spearman", exact = FALSE)
    data.frame(param = p, rho = ct$estimate, pval = ct$p.value)
  })
  bind_rows(results)
}

# ------------------------------------------------------------
# 4. Define scenarios
# ------------------------------------------------------------
scenarios <- list(
  list(label = "BRR (DALY) | 18\u201364, low transmission",
       age = "18-64", ar = "low",      brr_col = "brr_daly"),
  list(label = "BRR (DALY) | 65+, low transmission",
       age = "65+",   ar = "low",      brr_col = "brr_daly"),
  list(label = "BRR (DALY) | 18\u201364, moderate transmission",
       age = "18-64", ar = "moderate", brr_col = "brr_daly"),
  list(label = "BRR (DALY) | 65+, moderate transmission",
       age = "65+",   ar = "moderate", brr_col = "brr_daly"),
  list(label = "BRR (DALY) | 18\u201364, high transmission",
       age = "18-64", ar = "high",     brr_col = "brr_daly"),
  list(label = "BRR (DALY) | 65+, high transmission",
       age = "65+",   ar = "high",     brr_col = "brr_daly")
)

# ------------------------------------------------------------
# 5. Run Spearman for all scenarios
# ------------------------------------------------------------
spearman_results <- lapply(scenarios, function(s) {
  sub <- psa_merged %>%
    filter(age_group == s$age, ar_cat == s$ar)
  if (nrow(sub) == 0) {
    warning(paste("Zero rows for:", s$label)); return(NULL)
  }
  res <- compute_spearman(sub, s$brr_col, params_of_interest)
  if (is.null(res) || nrow(res) == 0) {
    warning(paste("No Spearman results for:", s$label)); return(NULL)
  }
  res$scenario <- s$label
  res$brr_col  <- s$brr_col
  res
}) %>% bind_rows()

cat("\nScenarios computed:", length(unique(spearman_results$scenario)), "\n")

# ------------------------------------------------------------
# 6. Per-scenario tornado plot function
#    Each plot has its OWN y factor — no cross-contamination
# ------------------------------------------------------------
TOP_N <- 12

plot_one_tornado <- function(scenario_label, top_n = TOP_N) {

  df <- spearman_results %>%
    filter(scenario == scenario_label) %>%
    mutate(
      param_label = dplyr::recode(param, !!!param_labels),
      abs_rho     = abs(rho),
      direction   = factor(
        ifelse(rho > 0, "Positive (favourable)", "Negative (unfavourable)"),
        levels = c("Positive (favourable)", "Negative (unfavourable)")
      )
    ) %>%
    arrange(desc(abs_rho)) %>%
    slice_head(n = top_n) %>%
    arrange(rho) %>%
    # Build a CLEAN factor for this scenario only
    mutate(param_label = factor(param_label, levels = param_label))

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
        "Positive (favourable)"   = "#56B4E9",  # Okabe-Ito sky blue
        "Negative (unfavourable)" = "#E69F00"   # Okabe-Ito orange
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
    # Separate theme() call to force bold title AFTER theme_minimal override
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
# 7. Build all 6 panels
# ------------------------------------------------------------
scenario_labels <- sapply(scenarios, `[[`, "label")
tornado_plots   <- lapply(scenario_labels, plot_one_tornado)
names(tornado_plots) <- scenario_labels

# ------------------------------------------------------------
# 8. Shared legend
# ------------------------------------------------------------
# Build a dummy plot just for the legend
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

# Assemble 6 panels in 2 columns
# NOTE: do NOT use & theme() here — it overrides plot.title face="bold"
# in individual panels. Legend is already suppressed in each panel.
combined <- wrap_plots(tornado_plots, ncol = 2) +
  plot_annotation(
    title    = "Spearman rank correlations between input parameters and BRR (DALY)",
    subtitle = paste0(
      "Top ", TOP_N, " parameters by |\u03c1| per scenario (low \u2192 moderate \u2192 high transmission).  ",
      "Positive \u03c1 \u2192 BRR increases (vaccination more favourable).  ",
      "Negative \u03c1 \u2192 BRR decreases \n(vaccination less favourable)."
    ),
    caption  = paste0(
      "Transmission categories: low AR <5%, moderate 5\u201310%, high \u226510%.  ",
      "BRR (DALY) = DALY averted / DALY from vaccine-associated SAE per 10,000.  BRR >1 = net benefit."
    ),
    theme = theme(
      plot.title    = element_text(size = 16, face = "bold",   family = "Calibri"),
      plot.subtitle = element_text(size = 12, colour = "grey20", family = "Calibri",
                                   margin = margin(b = 6)),
      plot.caption  = element_text(size = 12, colour = "grey20", family = "Calibri", hjust = 0)
    )
  )

# Display
print(combined)

# Save
output_dir <- "03_Manuscript/Submission/second_submission_lancet/suppl_figs/fig2_spearman"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(output_dir, "tornado_spearman.pdf"), combined,
       width = 12, height = 8, device = cairo_pdf)

