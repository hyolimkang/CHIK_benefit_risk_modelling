shp_dir <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/CHIK_vaccine_impact/CHIK_ORV_impact/00_Data/0_1_Raw/country_shape/gadm41_BRA_shp"
br_states <- st_read(file.path(shp_dir, "gadm41_BRA_1.shp"), quiet = TRUE)

days_levels <- c("7d","14d","30d","90d")
target_outcome <- "brr_daly"

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

ar_state <- psa_data_mid %>%
  group_by(state, draw) %>%
  summarise(AR_total_draw = mean(AR_total, na.rm=TRUE), .groups="drop") %>%
  group_by(state) %>%
  summarise(AR_total_med = median(AR_total_draw, na.rm=TRUE), .groups="drop")

ar_breaks <- as.numeric(quantile(ar_state$AR_total_med, probs = c(0, .25, .5, .75, 1), na.rm = TRUE))
ar_breaks <- unique(ar_breaks)

ar_labels <- paste0(
  fmt_pct(ar_breaks[-length(ar_breaks)]), "–", fmt_pct(ar_breaks[-1])
)

ar_state <- ar_state %>%
  mutate(
    AR_bin = cut(
      AR_total_med,
      breaks = ar_breaks,
      include.lowest = TRUE,
      labels = ar_labels
    )
  )

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
    age_group = factor(age_group, levels = c("1-11","12-17","18-64","65+")) %>% fct_drop()
  )

thr <- 0.8  # Pr(BRR>1) >= 0.8

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

ggsave("06_Results/brr_brazil_map_travelday.pdf", plot = fig1, width = 8, height = 6)
ggsave("06_Results/brr_brazil_map.pdf", plot = fig2, width = 7, height = 6)
