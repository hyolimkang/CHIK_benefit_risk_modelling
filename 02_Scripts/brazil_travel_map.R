shp_dir <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/CHIK_vaccine_impact/CHIK_ORV_impact/00_Data/0_1_Raw/country_shape/gadm41_BRA_shp"
br_states <- st_read(file.path(shp_dir, "gadm41_BRA_1.shp"), quiet = TRUE)

############################################################
## SETTINGS
############################################################

days_levels <- c("7d","14d","30d","90d")
target_outcome <- "brr_sae"

# Threshold is no longer used for "minimum days" map
# because we will show categorised Pr(BRR>1) and let readers decide.
thr <- 0.8

map_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
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

ar_max <- max(states_ar$AR_pct, na.rm = TRUE)
ar_breaks <- seq(0, ceiling(ar_max/5)*5, by = 5)


############################################################
## Pr(BRR>1) CALCULATION 
############################################################

brr_state_days_age <- psa_df %>%
  filter(!is.na(age_group)) %>%
  mutate(days = factor(days, levels = days_levels)) %>%
  group_by(state, age_group, days, draw) %>%
  summarise(brr_draw = mean(.data[[target_outcome]], na.rm = TRUE), .groups="drop") %>%
  group_by(state, age_group, days) %>%
  summarise(
    pr_gt1 = mean(brr_draw > 1, na.rm = TRUE),
    .groups="drop"
  ) %>%
  mutate(
    age_group = factor(age_group, levels = c("1-11","12-17","18-64","65+")) %>% fct_drop(),
    state_key = norm_state(state)
  ) %>%
  dplyr::select(state_key, age_group, days, pr_gt1)


############################################################
## IMPORTANT UPDATE:
## COMPLETE GRID so NA states remain on the map as grey
############################################################

# Create a full grid of (all states) x (all ages) x (all day levels)
# so states with no data are still drawn (as NA -> grey).
states_pr_full <- states_base %>%
  tidyr::crossing(
    age_group = factor(c("1-11","12-17","18-64","65+"),
                       levels = c("1-11","12-17","18-64","65+")),
    days      = factor(days_levels, levels = days_levels)
  ) %>%
  left_join(brr_state_days_age, by = c("state_key","age_group","days"))

############################################################
## IMPORTANT UPDATE:
## Categorise Pr(BRR>1) so readers can pick their own threshold
############################################################

pr_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
pr_labels <- c("0–20%", "20–40%", "40–60%", "60–80%", "≥80%")

states_pr_full <- states_pr_full %>%
  mutate(
    pr_cat = cut(
      pr_gt1,
      breaks = pr_breaks,
      include.lowest = TRUE,
      labels = pr_labels
    ),
    pr_cat = forcats::fct_relevel(pr_cat, pr_labels)
  )%>%
  sf::st_as_sf()

class(states_pr_full)
names(states_pr_full)
sf::st_geometry(states_pr_full) %>% head()
inherits(states_pr_full, "sf")
############################################################
## MAPS
############################################################

# Attack rate map (continuous)
p_ar2 <- ggplot(states_ar) +
  geom_sf(aes(fill = AR_pct), colour="white", linewidth=0.08) +  # thinner borders
  labs(fill = "Attack rate (%)") +
  coord_sf(datum = NA) +
  map_theme +
  scale_fill_distiller(
    palette = "Reds",
    direction = 1,
    limits = c(0, max(ar_breaks)),
    oob = scales::squish,
    na.value = "grey80",
    breaks = scales::pretty_breaks(n = 5),
    labels = function(x) paste0(round(x), "%")
  ) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  theme(plot.margin = margin(t = 2, r = -25, b = 2, l = 2, unit = "mm"))


# UPDATED: Pr(BRR>1) map as categorical bins
p_pr_cat <- ggplot(states_pr_full) +
  geom_sf(aes(fill = pr_cat), colour="white", linewidth=0.01) +  # thinner borders
  facet_grid(days ~ age_group) +
  labs(fill = "Pr(BRR>1)") +
  coord_sf(datum = NA) +
  map_theme +
  scale_fill_viridis_d(
    option = "C",
    direction = 1,
    na.value = "grey80",   # NA states remain visible in grey for geographic context
    #drop = FALSE,
    drop = TRUE
  ) +
  theme(legend.position = "bottom", legend.box = "horizontal")+
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = -25, unit = "mm"))
############################################################
## FINAL: SINGLE FIGURE
## (Remove the "minimum travel days" figure and keep one combined figure)
############################################################

fig_one <- (p_ar2 + p_pr_cat) + 
  plot_layout(widths = c(1.1, 4.5), guides = "keep") + 
  plot_annotation(title = "Benefit-risk ratio of traveller vaccination: Pr(BRR(SAE)>1) by travel duration and age group", 
                  theme = theme(plot.title = element_text(size = 10)) 
                  )
fig_one <- fig_one & theme(legend.position = "bottom", legend.box = "horizontal") +
           theme(text = element_text(family = "Calibri"))

fig_one

ggsave("06_Results/brr_brazil_map_fig_travel_sae.pdf", plot = fig_one, width = 11, height = 8, device = cairo_pdf)


## this is the end of traveller scenario--------------------------------------------------------------

### ORI ------------------------------------------------------------------------

age_levels <- c("1-11","12-17","18-64","65+")

daly_dat <- draw_level_xy_true %>%
  filter(outcome == "SAE") %>%   # IMPORTANT: subset to DALY only
  mutate(
    state_key = norm_state(Region),
    AgeCat    = factor(AgeCat, levels = age_levels),
    setting   = factor(setting)
  ) %>%
  # Compute BRR safely (avoid division by 0)
  mutate(
    brr_sae = dplyr::if_else(
      is.na(sae_10k) | sae_10k == 0,
      NA_real_,
      y_10k / sae_10k
    )
  )

############################################################
## 2) Summarise Pr(BRR>1) across draws for each state/setting/age
############################################################
pr_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
pr_labels <- c("0–20%", "20–40%", "40–60%", "60–80%", "≥80%")

brr_pr_ve_age <- daly_dat %>%
  group_by(state_key, VE_label, AgeCat, draw_id) %>%
  summarise(brr_draw = mean(brr_sae, na.rm = TRUE), .groups = "drop") %>%
  group_by(state_key, VE_label, AgeCat) %>%
  summarise(pr_gt1 = mean(brr_draw > 1, na.rm = TRUE), .groups = "drop")

brr_pr_ve_age <- brr_pr_ve_age %>%
  mutate(
    pr_cat = cut(pr_gt1, breaks = pr_breaks, include.lowest = TRUE, labels = pr_labels),
    pr_cat = forcats::fct_relevel(pr_cat, pr_labels),
    pr_cat = forcats::fct_na_value_to_level(pr_cat, level = "No data")
  )
############################################################
## 3) Categorise probability so readers can choose threshold
############################################################
ve_levels  <- sort(unique(brr_pr_ve_age$VE_label))
age_levels <- c("1-11","12-17","18-64","65+")

states_pr_full_outbreak <- states_base %>%
  tidyr::crossing(
    AgeCat   = factor(age_levels, levels = age_levels),
    VE_label = factor(ve_levels, levels = ve_levels)
  ) %>%
  left_join(brr_pr_ve_age, by = c("state_key","VE_label","AgeCat")) %>%
  sf::st_as_sf()

bb <- sf::st_bbox(states_base)

states_pr_full_outbreak$pr_cat <- droplevels(states_pr_full_outbreak$pr_cat)


# 5) Plot
p_brr_cat <- ggplot(states_pr_full_outbreak) +
  geom_sf(aes(fill = pr_cat), colour="white", linewidth=0.01) +
  facet_grid(VE_label ~ AgeCat) +
  labs(fill = "Pr(BRR(SAE)>1)") +
  coord_sf(datum = NA) +
  map_theme +
  scale_fill_viridis_d(
    option = "C",
    direction = 1,
    na.value = "grey80",
    drop = FALSE
  ) +
  theme(legend.position = "bottom", legend.box = "horizontal")+
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 0, unit = "mm"))

p_brr_cat

fig_one <- (p_ar2 + p_brr_cat) + 
  plot_layout(widths = c(1.1, 4.5), guides = "keep") + 
  plot_annotation(title = "Benefit-risk ratio of outbreak response immunisation: Pr(BRR(SAE)>1) by vaccine protection mechanism and age group", 
                  theme = theme(plot.title = element_text(size = 10)) 
  )

fig_one <- fig_one & theme(legend.position = "bottom", legend.box = "horizontal") + 
           theme(text = element_text(family = "Calibri"))
fig_one 

ggsave("06_Results/brr_brazil_map_ori_sae.pdf", plot = fig_one, width = 12, height = 6, device = cairo_pdf)

 

#### END-------------------------------------------------------------------------