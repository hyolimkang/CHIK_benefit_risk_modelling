shp_dir <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/CHIK_vaccine_impact/CHIK_ORV_impact/00_Data/0_1_Raw/country_shape/gadm41_BRA_shp"
br_states <- st_read(file.path(shp_dir, "gadm41_BRA_1.shp"), quiet = TRUE)

############################################################
## SETTINGS
############################################################

days_levels <- c("7d","14d","30d","90d")
target_outcome <- "brr_daly"

# Threshold is no longer used for "minimum days" map
# because we will show categorised Pr(BRR>1) and let readers decide.
thr <- 0.8

map_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = ggplot2::margin(2, 2, 2, 2),
    legend.position = "bottom",
    legend.box = "horizontal"
  )
############################################################
## HELPERS
############################################################

norm_state <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_trim() %>%
    stringi::stri_trans_general("Latin-ASCII") %>%
    str_replace_all("[^a-z0-9]+", "")
}

fmt_pct <- function(x, digits = 1){
  if (max(x, na.rm = TRUE) <= 1) scales::percent(x, accuracy = 10^(-digits))
  else paste0(round(x, digits), "%")
}

# Keep only age groups that exist in PSA input (e.g., >=18 if under-18 removed).
age_levels_plot <- psa_df %>%
  dplyr::transmute(age_group = dplyr::if_else(age_group == "65", "65+", as.character(age_group))) %>%
  dplyr::filter(!is.na(age_group), age_group %in% c("18-64", "65+")) %>%
  dplyr::distinct(age_group) %>%
  dplyr::pull(age_group)


############################################################
## ATTACK RATE 
############################################################

ar_state <- psa_df %>%
  group_by(state, draw) %>%
  summarise(AR_total_draw = mean(AR_total, na.rm=TRUE), .groups="drop") %>%
  group_by(state) %>%
  summarise(AR_total_med = median(AR_total_draw, na.rm=TRUE), .groups="drop")

# Keep AR as continuous % (recommended for interpretation)
states_base <- br_states %>%
  mutate(state_key = norm_state(NAME_1))

states_ar <- states_base %>%
  left_join(ar_state %>% mutate(state_key = norm_state(state)), by="state_key") %>%
  mutate(AR_pct = 100 * AR_total_med)

ar_range <- range(states_ar$AR_pct[is.finite(states_ar$AR_pct)], na.rm = TRUE)
if (!all(is.finite(ar_range))) ar_range <- c(0, 1)
if (diff(ar_range) == 0) ar_range <- ar_range + c(-0.5, 0.5)


############################################################
## Pr(BRR>1) CALCULATION 
############################################################

p_ar2 <- ggplot(states_ar) +
  geom_sf(aes(fill = AR_pct), colour="white", linewidth=0.08) +  # thinner borders
  labs(fill = "Attack rate (%)") +
  coord_sf(datum = NA) +
  map_theme +
  scale_fill_distiller(
    palette = "Reds",
    direction = 1,
    limits = ar_range,
    oob = scales::squish,
    na.value = "grey80",
    breaks = scales::pretty_breaks(n = 5),
    labels = function(x) paste0(round(x), "%")
  ) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  theme(plot.margin =  ggplot2::margin(t = 2, r = -25, b = 2, l = 2, unit = "mm"))

pr_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
pr_labels <- c("0–20%", "20–40%", "40–60%", "60–80%", ">80%")
pr_palette <- setNames(
  viridisLite::viridis(length(pr_labels), option = "C", direction = 1),
  pr_labels
)

build_travel_probability_map <- function(target_outcome = "brr_daly",
                                         save_path = NULL,
                                         width = 14, height = 8) {
  travel_outcome_label <- dplyr::case_when(
    identical(target_outcome, "brr_sae") ~ "SAE",
    identical(target_outcome, "brr_death") ~ "Death",
    TRUE ~ "DALY"
  )

  brr_state_days_age <- psa_df %>%
    dplyr::mutate(age_group = dplyr::if_else(age_group == "65", "65+", as.character(age_group))) %>%
    dplyr::filter(!is.na(age_group), age_group %in% age_levels_plot) %>%
    dplyr::mutate(days = factor(days, levels = days_levels)) %>%
    dplyr::group_by(state, age_group, days, draw) %>%
    dplyr::summarise(brr_draw = mean(.data[[target_outcome]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(state, age_group, days) %>%
    dplyr::summarise(pr_gt1 = mean(brr_draw > 1, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      age_group = factor(age_group, levels = age_levels_plot) %>% forcats::fct_drop(),
      state_key = norm_state(state)
    ) %>%
    dplyr::select(state_key, age_group, days, pr_gt1)

  age_levels_outcome <- brr_state_days_age %>%
    dplyr::group_by(age_group) %>%
    dplyr::summarise(has_data = any(is.finite(pr_gt1) & !is.na(pr_gt1)), .groups = "drop") %>%
    dplyr::filter(has_data) %>%
    dplyr::pull(age_group) %>%
    as.character()
  if (length(age_levels_outcome) == 0) age_levels_outcome <- age_levels_plot

  brr_state_days_age <- brr_state_days_age %>%
    dplyr::filter(as.character(age_group) %in% age_levels_outcome) %>%
    dplyr::mutate(age_group = factor(as.character(age_group), levels = age_levels_outcome))

  states_pr_full <- states_base %>%
    tidyr::crossing(
      age_group = factor(age_levels_outcome, levels = age_levels_outcome),
      days = factor(days_levels, levels = days_levels)
    ) %>%
    dplyr::left_join(brr_state_days_age, by = c("state_key", "age_group", "days")) %>%
    dplyr::mutate(
      pr_cat = cut(pr_gt1, breaks = pr_breaks, include.lowest = TRUE, labels = pr_labels),
      pr_cat = factor(as.character(forcats::fct_relevel(pr_cat, pr_labels)), levels = pr_labels)
    ) %>%
    sf::st_as_sf()

  p_pr_cat <- ggplot(states_pr_full) +
    geom_sf(aes(fill = pr_cat), colour = "white", linewidth = 0.01) +
    facet_grid(age_group ~ days) +
    labs(fill = paste0("Pr(BRR(", travel_outcome_label, ")>1)")) +
    coord_sf(datum = NA) +
    map_theme +
    scale_fill_manual(
      values = pr_palette,
      na.value = "grey80",
      breaks = pr_labels,
      limits = pr_labels,
      drop = FALSE,
      guide = ggplot2::guide_legend(
        override.aes = list(fill = unname(pr_palette[pr_labels]), colour = NA)
      )
    ) +
    theme(legend.position = "bottom", legend.box = "horizontal") +
    theme(plot.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = -25, unit = "mm"))

  fig_out <- (p_ar2 + p_pr_cat) + patchwork::plot_layout(widths = c(1.1, 4.5), guides = "keep")
  fig_out <- fig_out & theme(legend.position = "bottom", legend.box = "horizontal") +
    theme(text = element_text(family = "Calibri", size = 13))

  if (!is.null(save_path)) {
    ggsave(save_path, plot = fig_out, width = width, height = height, device = cairo_pdf)
  }

  list(fig = fig_out, map = p_pr_cat, pr_table = brr_state_days_age)
}

travel_map_daly <- build_travel_probability_map(
  target_outcome = "brr_daly",
  save_path = "06_Results/brr_brazil_map_fig_travel_daly.pdf"
)
travel_map_sae <- build_travel_probability_map(
  target_outcome = "brr_sae",
  save_path = "06_Results/brr_brazil_map_fig_travel_sae.pdf"
)
travel_map_death <- build_travel_probability_map(
  target_outcome = "brr_death",
  save_path = "06_Results/brr_brazil_map_fig_travel_death.pdf"
)
travel_map_daly$fig


## this is the end of traveller scenario--------------------------------------------------------------

### ORI ------------------------------------------------------------------------

age_levels_ori <- c("18-64", "65+")
pr_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
pr_labels <- c("0–20%", "20–40%", "40–60%", "60–80%", ">80%")
pr_palette <- setNames(
  viridisLite::viridis(length(pr_labels), option = "C", direction = 1),
  pr_labels
)

if (!exists("draw_level_xy_serostatus_filtered", inherits = FALSE)) {
  if (!exists("draw_level_xy_serostatus", inherits = FALSE)) {
    stop("Run brazil_all_draws_ori_v2.R first (need draw_level_xy_serostatus_filtered).")
  }
  draw_level_xy_serostatus_filtered <- draw_level_xy_serostatus %>%
    dplyr::filter(.data$RR_seropos == 0)
}

build_ori_probability_map <- function(target_outcome = c("SAE", "DALY", "Death"),
                                      save_path = NULL,
                                      width = 12, height = 6) {
  target_outcome <- match.arg(target_outcome)
  src <- draw_level_xy_serostatus_filtered

  ori_dat <- src %>%
    dplyr::filter(.data$outcome == target_outcome, .data$AgeCat %in% age_levels_ori) %>%
    dplyr::mutate(
      state_key = norm_state(.data$Region),
      AgeCat = as.character(.data$AgeCat),
      VE_label = factor(.data$VE_label)
    )

  ori_dat <- ori_dat %>%
    dplyr::mutate(AgeCat = factor(.data$AgeCat, levels = age_levels_ori))

  # Build BRR series for both base and adjusted risks when available.
  if (all(c("brr_base", "brr_adj") %in% names(ori_dat))) {
    ori_dat <- ori_dat %>%
      dplyr::transmute(
        state_key, VE_label, AgeCat, draw_id,
        brr_base = as.numeric(.data$brr_base),
        brr_adj = as.numeric(.data$brr_adj)
      ) %>%
      tidyr::pivot_longer(
        cols = c("brr_base", "brr_adj"),
        names_to = "brr_type",
        values_to = "brr_value"
      ) %>%
      dplyr::mutate(
        brr_type = dplyr::recode(.data$brr_type,
                                 "brr_base" = "Base",
                                 "brr_adj" = "Adjusted"),
        brr_type = factor(.data$brr_type, levels = c("Base", "Adjusted"))
      )
  } else if ("brr_adj" %in% names(ori_dat)) {
    ori_dat <- ori_dat %>%
      dplyr::mutate(
        brr_type = factor("Adjusted", levels = c("Base", "Adjusted")),
        brr_value = as.numeric(.data$brr_adj)
      )
  } else if ("x_10k_adj" %in% names(ori_dat) && "y_10k" %in% names(ori_dat)) {
    ori_dat <- ori_dat %>%
      dplyr::mutate(
        brr_type = factor("Adjusted", levels = c("Base", "Adjusted")),
        brr_value = dplyr::if_else(
          is.na(.data$x_10k_adj) | .data$x_10k_adj == 0,
          NA_real_,
          .data$y_10k / .data$x_10k_adj
        )
      )
  } else if ("x_10k" %in% names(ori_dat) && "y_10k" %in% names(ori_dat)) {
    ori_dat <- ori_dat %>%
      dplyr::mutate(
        brr_type = factor("Adjusted", levels = c("Base", "Adjusted")),
        brr_value = dplyr::if_else(
          is.na(.data$x_10k) | .data$x_10k == 0,
          NA_real_,
          .data$y_10k / .data$x_10k
        )
      )
  } else {
    stop("No compatible BRR columns found in draw_level_xy_serostatus_filtered (expected brr_adj, or y_10k with x_10k_adj/x_10k).")
  }

  brr_pr_ve_age <- ori_dat %>%
    dplyr::group_by(.data$state_key, .data$brr_type, .data$VE_label, .data$AgeCat, .data$draw_id) %>%
    dplyr::summarise(brr_draw = mean(.data$brr_value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(.data$state_key, .data$brr_type, .data$VE_label, .data$AgeCat) %>%
    dplyr::summarise(pr_gt1 = mean(.data$brr_draw > 1, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      pr_cat = cut(.data$pr_gt1, breaks = pr_breaks, include.lowest = TRUE, labels = pr_labels),
      pr_cat = forcats::fct_relevel(.data$pr_cat, pr_labels)
    )

  age_levels_ori_outcome <- brr_pr_ve_age %>%
    dplyr::group_by(.data$AgeCat) %>%
    dplyr::summarise(has_data = any(is.finite(.data$pr_gt1) & !is.na(.data$pr_gt1)), .groups = "drop") %>%
    dplyr::filter(.data$has_data) %>%
    dplyr::pull(.data$AgeCat) %>%
    as.character()
  if (length(age_levels_ori_outcome) == 0) age_levels_ori_outcome <- age_levels_ori

  brr_pr_ve_age <- brr_pr_ve_age %>%
    dplyr::filter(as.character(.data$AgeCat) %in% age_levels_ori_outcome) %>%
    dplyr::mutate(AgeCat = factor(as.character(.data$AgeCat), levels = age_levels_ori_outcome))
  
  single_age_row <- length(age_levels_ori_outcome) == 1

  available_types <- base::intersect(
    c("Base", "Adjusted"),
    unique(as.character(brr_pr_ve_age$brr_type))
  )
  if (length(available_types) == 0) {
    stop("No BRR map panels available after filtering by brr_type.")
  }

  ve_levels_all <- sort(unique(as.character(brr_pr_ve_age$VE_label)))

  ori_plot_data <- states_base %>%
    tidyr::crossing(
      brr_type = factor(available_types, levels = available_types),
      AgeCat   = factor(age_levels_ori_outcome, levels = age_levels_ori_outcome),
      VE_label = factor(ve_levels_all, levels = ve_levels_all)
    ) %>%
    dplyr::left_join(
      brr_pr_ve_age %>% dplyr::mutate(brr_type = factor(brr_type, levels = available_types)),
      by = c("state_key", "brr_type", "VE_label", "AgeCat")
    ) %>%
    dplyr::mutate(pr_cat = factor(as.character(.data$pr_cat), levels = pr_labels)) %>%
    sf::st_as_sf()

  nested_facet <- if (single_age_row) {
    ggh4x::facet_nested(. ~ brr_type + VE_label)
  } else {
    ggh4x::facet_nested(AgeCat ~ brr_type + VE_label)
  }

  p_brr_cat <- ggplot(ori_plot_data) +
    geom_sf(aes(fill = pr_cat), colour = "white", linewidth = 0.01) +
    nested_facet +
    labs(fill = paste0("Pr(BRR(", target_outcome, ")>1)")) +
    coord_sf(datum = NA) +
    map_theme +
    scale_fill_manual(
      values = pr_palette,
      na.value = "grey80",
      breaks = pr_labels,
      limits = pr_labels,
      drop = FALSE,
      guide = ggplot2::guide_legend(
        title.position = "left",
        title.vjust = 0.9,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(fill = unname(pr_palette[pr_labels]), colour = NA)
      )
    ) +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      plot.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 0, unit = "mm")
    )

  fig_out <- (p_ar2 + p_brr_cat) +
    patchwork::plot_layout(widths = c(1.1, 4.5), guides = "keep")

  fig_out <- fig_out &
    theme(
      text = element_text(family = "Calibri"),
      legend.box = "horizontal",
      legend.title = element_text(margin = ggplot2::margin(r = 8)),
      legend.spacing.x = grid::unit(2, "mm"),
      legend.key.width = grid::unit(4, "mm"),
      legend.key.height = grid::unit(4, "mm")
    )

  if (!is.null(save_path)) {
    ggsave(save_path, plot = fig_out, width = width, height = height, device = cairo_pdf)
  }

  list(fig = fig_out, map = p_brr_cat, pr_table = brr_pr_ve_age)
}

# Default outputs for all three outcomes.
ori_map_sae <- build_ori_probability_map("SAE", save_path = "06_Results/brr_brazil_map_ori_sae.pdf")
ori_map_daly <- build_ori_probability_map("DALY", save_path = "06_Results/brr_brazil_map_ori_daly.pdf")
ori_map_death <- build_ori_probability_map("Death", save_path = "06_Results/brr_brazil_map_ori_death.pdf")

ori_map_sae$fig

 

#### END-------------------------------------------------------------------------