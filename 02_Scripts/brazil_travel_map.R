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

ar_state <- psa_data_mid %>%
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

brr_state_days_age <- psa_data_mid %>%
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
  theme(legend.position = "bottom", legend.box = "horizontal") 


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
    drop = FALSE
  ) +
  theme(legend.position = "bottom", legend.box = "horizontal") 
############################################################
## FINAL: SINGLE FIGURE
## (Remove the "minimum travel days" figure and keep one combined figure)
############################################################

fig_one <- (p_ar2 + p_pr_cat) + 
  plot_layout(widths = c(1.1, 4.5), guides = "collect") + 
  plot_annotation(title = "Benefit-risk ratio of traveller vaccination: Pr(BRR(DALY)>1) by travel duration and age group", 
                  theme = theme(plot.title = element_text(size = 10)) 
                  )
fig_one <- fig_one & theme(legend.position = "bottom", legend.box = "horizontal")
fig_one

ggsave("06_Results/brr_brazil_map_fig_one.pdf", plot = fig_one, width = 8, height = 6)


## this is the end of traveller scenario--------------------------------------------------------------

### ORI ------------------------------------------------------------------------

age_levels <- c("1-11","12-17","18-64","65+")

daly_dat <- draw_level_xy_true %>%
  filter(outcome == "DALY") %>%   # IMPORTANT: subset to DALY only
  mutate(
    state_key = norm_state(Region),
    AgeCat    = factor(AgeCat, levels = age_levels),
    setting   = factor(setting)
  ) %>%
  # Compute BRR safely (avoid division by 0)
  mutate(
    brr_daly = dplyr::if_else(
      is.na(x_daly_10k) | x_daly_10k == 0,
      NA_real_,
      y_10k / x_daly_10k
    )
  )

############################################################
## 2) Summarise Pr(BRR>1) across draws for each state/setting/age
############################################################

brr_pr <- daly_dat %>%
  group_by(state_key, setting, AgeCat, draw_id) %>%
  summarise(brr_draw = mean(brr_daly, na.rm = TRUE), .groups = "drop") %>%
  group_by(state_key, setting, AgeCat) %>%
  summarise(
    pr_gt1 = mean(brr_draw > 1, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
## 3) Categorise probability so readers can choose threshold
############################################################

pr_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
pr_labels <- c("0–20%", "20–40%", "40–60%", "60–80%", "≥80%")

brr_pr <- brr_pr %>%
  mutate(
    pr_cat = cut(
      pr_gt1,
      breaks = pr_breaks,
      include.lowest = TRUE,
      labels = pr_labels
    ),
    pr_cat = forcats::fct_explicit_na(pr_cat, na_level = "No data")
  )

states_brr_full <- states_base %>%
  tidyr::crossing(
    setting = sort(unique(brr_pr$setting)),
    AgeCat  = factor(age_levels, levels = age_levels)
  ) %>%
  left_join(brr_pr, by = c("state_key","setting","AgeCat"))%>%
  sf::st_as_sf()

bb <- sf::st_bbox(states_base)

p_brr_cat <- g# 1) DALY만 필터 + BRR(DALY)= y_10k/x_daly_10k
age_levels <- c("1-11","12-17","18-64","65+")

daly_dat <- draw_level_xy_true %>%
  filter(outcome == "DALY") %>%
  mutate(
    state_key = norm_state(Region),
    AgeCat    = factor(AgeCat, levels = age_levels),
    setting   = factor(setting),
    brr_daly  = dplyr::if_else(is.na(x_daly_10k) | x_daly_10k == 0, NA_real_, y_10k / x_daly_10k)
  )

# 2) Pr(BRR>1)
brr_pr <- daly_dat %>%
  group_by(state_key, setting, AgeCat, draw_id) %>%
  summarise(brr_draw = mean(brr_daly, na.rm = TRUE), .groups = "drop") %>%
  group_by(state_key, setting, AgeCat) %>%
  summarise(pr_gt1 = mean(brr_draw > 1, na.rm = TRUE), .groups = "drop")

# 3) Categorise Pr, keep NA as NA (IMPORTANT)
pr_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
pr_labels <- c("0–20%", "20–40%", "40–60%", "60–80%", "≥80%")

brr_pr <- brr_pr %>%
  mutate(
    pr_cat = cut(pr_gt1, breaks = pr_breaks, include.lowest = TRUE, labels = pr_labels)
  )

# 4) sf-safe expansion (sf on the left)
states_brr_full <- states_base %>%
  sf::st_drop_geometry() %>%
  tidyr::crossing(
    setting = sort(unique(brr_pr$setting)),
    AgeCat  = factor(age_levels, levels = age_levels)
  ) %>%
  left_join(brr_pr, by = c("state_key","setting","AgeCat")) %>%
  # 2) Re-attach geometry from the sf base
  left_join(
    states_base %>% dplyr::select(state_key, geometry),
    by = "state_key"
  ) %>%
  sf::st_as_sf()

# (optional but recommended) lock bbox to full Brazil
bb <- sf::st_bbox(states_base)

# 5) Plot: EXACTLY the same style as your previous code
p_brr_cat <- ggplot(states_brr_full) +
  geom_sf(aes(fill = pr_cat), colour="white", linewidth=0.01) +
  facet_grid(setting ~ AgeCat) +
  labs(fill = "Pr(BRR(DALY)>1)") +
  coord_sf(
    datum = NA,
    xlim = c(bb["xmin"], bb["xmax"]),
    ylim = c(bb["ymin"], bb["ymax"])
  ) +
  map_theme +
  scale_fill_viridis_d(
    option = "C",
    direction = 1,
    na.value = "grey80",  # NA states visible in grey
    drop = FALSE
  ) +
  theme(legend.position = "bottom")

p_brr_cat

fig_one <- (p_ar2 + p_brr_cat) + 
  plot_layout(widths = c(1.1, 4.5), guides = "collect") + 
  plot_annotation(title = "Benefit-risk ratio of outbreak response immunisation: Pr(BRR(DALY)>1) by attack rate and age group", 
                  theme = theme(plot.title = element_text(size = 10)) 
  )

fig_one <- fig_one & theme(legend.position = "bottom", legend.box = "horizontal")
fig_one

ggsave("06_Results/brr_brazil_map_ori.pdf", plot = fig_one, width = 8, height = 6)

right_column <- (p_pr_cat / p_brr_cat)

## combined figs
p_brr_cat_no_legend <- p_brr_cat + theme(legend.position = "none")

design <- "
  ABB
  ACC
"
combined_fig <- wrap_plots(
  A = p_ar2, 
  B = p_pr_cat, 
  C = p_brr_cat_no_legend, 
  design = design
) +
  plot_layout(
    widths = c(1, 3.5), 
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = 'A',
    title = "Comparative benefit-risk assessment of Ixchiq vaccination strategies",
    subtitle = "A: Regional attack rate; B: Traveller vaccination; C: Outbreak response vaccination.",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  ) & 
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )

combined_fig
#### old code -----------------------------------------------------------------
min_days_state_age <- brr_state_days_age %>%
  arrange(state, age_group, days) %>%
  group_by(state, age_group) %>%
  summarise(
    min_days = {
      ok <- as.character(days[pr_gt1 >= thr])
      if (length(ok) == 0) "Not favourable" else ok[1]
    },
    .groups = "drop"
  ) %>%
  mutate(
    min_days_cat = factor(min_days, levels = c(days_levels, "Not favourable")),
    age_group = factor(age_group, levels = c("1-11","12-17","18-64","65+")) %>% fct_drop()
  )

states_base <- br_states %>%
  mutate(state_key = norm_state(NAME_1))

states_ar <- states_base %>%
  left_join(ar_state %>% mutate(state_key = norm_state(state)), by="state_key")

states_min <- states_base %>%
  left_join(min_days_state_age %>% mutate(state_key = norm_state(state)), by="state_key")

states_pr <- states_base %>%
  left_join(brr_state_days_age %>% mutate(state_key = norm_state(state)), by="state_key")

map_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(2, 2, 2, 2),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

##################

p_min <- ggplot(states_min %>% filter(!is.na(age_group))) +
  geom_sf(aes(fill = min_days_cat), colour="white", linewidth=0.15) +
  facet_wrap(~ age_group, nrow = 1) +
  labs(fill = "Minimum travel days") +
  coord_sf(datum = NA) +
  map_theme +
  scale_fill_discrete(na.value = "grey80") +
  theme(legend.position = "bottom")  

p_ar2 <- ggplot(states_ar) +
  geom_sf(aes(fill = AR_bin), colour="white", linewidth=0.15) +
  labs(fill = "Attack rate (median)") +
  coord_sf(datum = NA) +
  map_theme +
  scale_fill_discrete(na.value = "grey80") +
  theme(legend.position = "bottom")

p_pr2 <- ggplot(states_pr %>% filter(!is.na(age_group), !is.na(days))) +
  geom_sf(aes(fill = pr_gt1), colour="white", linewidth=0.15) +
  facet_grid(days ~ age_group) +
  labs(fill = "Pr(BRR>1)") +
  coord_sf(datum = NA) +
  map_theme +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    oob = scales::squish,
    na.value = "grey80",
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme(legend.position = "bottom")


states_ar <- states_ar %>%
  mutate(AR_pct = 100 * AR_total_med)

ar_max <- max(states_ar$AR_pct, na.rm = TRUE)
ar_breaks <- seq(0, ceiling(ar_max/5)*5, by = 5)

p_ar2 <- ggplot(states_ar) +
  geom_sf(aes(fill = AR_pct), colour="white", linewidth=0.15) +
  labs(fill = "Attack rate (%)") +
  coord_sf(datum = NA) +
  map_theme +
  scale_fill_distiller(
    palette = "Reds",
    direction = 1,
    #breaks = ar_breaks,
    limits = c(0, max(ar_breaks)),
    oob = scales::squish,
    na.value = "grey80",
    breaks = scales::pretty_breaks(n = 5),
    labels = function(x) paste0(round(x), "%")
  ) +
  theme(legend.position = "bottom")

fig1 <- (p_ar2 + p_min) + 
  plot_layout(
    widths = c(1.5, 5),    
    guides = "collect"
  ) + 
  plot_annotation(
    title = "Minimum recommended travel duration = earliest day where Pr(BRR>1) ≥ 0.8"
  )

fig1 <- fig1 & theme(legend.position = "bottom", legend.box = "horizontal")

fig1

fig2 <- (p_ar2 + p_pr2) +
  plot_layout(widths = c(1.1, 4.5), guides = "collect") +
  plot_annotation(
    title = paste0("Pr(BRR>1) by travel duration and age group (BRR for DALY)"),
    theme = theme(plot.title = element_text(size = 10))
  )

fig2 <- fig2 & theme(legend.position = "bottom", legend.box = "horizontal")
fig2

ggsave("06_Results/brr_brazil_map_ori.pdf", plot = fig2, width = 7, height = 6)
