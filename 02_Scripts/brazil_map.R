shp_dir <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/CHIK_vaccine_impact/CHIK_ORV_impact/00_Data/0_1_Raw/country_shape/gadm41_BRA_shp"
br_states <- st_read(file.path(shp_dir, "gadm41_BRA_1.shp"), quiet = TRUE)

library(dplyr)


ar_state <- psa_df %>%
  group_by(state, draw) %>%
  summarise(AR_total_draw = mean(AR_total, na.rm=TRUE), .groups="drop") %>%
  group_by(state) %>%
  summarise(
    AR_total_med = median(AR_total_draw, na.rm=TRUE),
    .groups="drop"
  )

target_days <- "90d"
target_age  <- unique(psa_df$age_group)[1]  

brr14_state <- psa_df %>%
  filter(days == target_days, age_group == target_age) %>%
  group_by(state, draw) %>%
  summarise(brr_draw = mean(brr_daly, na.rm=TRUE), .groups="drop") %>%
  group_by(state) %>%
  summarise(
    brr_med = median(brr_draw, na.rm=TRUE),
    pr_gt1  = mean(brr_draw > 1, na.rm=TRUE),
    .groups="drop"
  )


states_join <- br_states %>%
  mutate(state_key = norm(NAME_1)) %>%
  left_join(ar_state  %>% mutate(state_key = norm(state)), by="state_key") %>%
  left_join(brr14_state %>% mutate(state_key = norm(state)), by="state_key")

map_long <- states_join %>%
  select(geometry, AR_total_med, pr_gt1) %>%
  pivot_longer(c(AR_total_med, pr_gt1), names_to="metric", values_to="value") %>%
  mutate(metric = recode(metric,
                         AR_total_med="Attack rate (median)",
                         pr_gt1="Pr(BRR>1) for 90d"))

ggplot(map_long) +
  geom_sf(aes(fill=value), colour="white", linewidth=0.15) +
  facet_wrap(~metric, ncol=2) +
  labs(fill=NULL) +
  theme_minimal()