# setwd
setwd("C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/CHIK_benefit_risk")

## hosp and death
# hosp and death 
hosp  <- readRDS("01_Data/hosp_chikvna.RDS")
death <- readRDS("01_Data/death_chikvna.RDS")

total_n <- hosp %>% group_by(age_group)%>%
  mutate(total_n = sum(n)) %>% 
  mutate(hosp_rate = n / total_n)

hosp$total_n <- total_n$total_n

death_age <- death %>%
  # One row per age_group: get deaths (deaht == TRUE) and the denominator (total_n)
  group_by(age_group) %>%
  summarise(
    death_n = sum(n[deaht], na.rm = TRUE), # deaht == TRUE rows
    denom_n = max(total_n, na.rm = TRUE), # total_n repeated across TRUE/FALSE rows
    death_rate = death_n / denom_n, # age-specific death rate
    .groups = "drop"
  )

death_band <- death_age %>%
  mutate(
    age_band = case_when(
      age_group >= 1 & age_group <= 11 ~ "1-11",
      age_group >= 12 & age_group <= 17 ~ "12-17",
      age_group >= 18 & age_group <= 64 ~ "18-64",
      age_group >= 65 ~ "65+",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_band)) %>%
  group_by(age_band) %>%
  summarise(
    death_n = sum(death_n),
    denom_n = sum(denom_n),
    death_rate = death_n / denom_n, # weighted/overall rate within band
    .groups = "drop"
  )

hosp <- hosp %>%
  # Collapse to one row per age_group (to avoid duplicated total_n)
  group_by(age_group) %>%
  summarise(
    hosp_n = sum(n[hosp], na.rm = TRUE),
    denom_n = max(total_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Create age bands
  mutate(
    age_band = case_when(
      age_group >= 1 & age_group <= 11 ~ "1-11",
      age_group >= 12 & age_group <= 17 ~ "12-17",
      age_group >= 18 & age_group <= 64 ~ "18-64",
      age_group >= 65 ~ "65+",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_band)) %>%
  group_by(age_band) %>%
  summarise(
    # Weighted/overall hospitalisation rate for the band
    hosp_rate = sum(hosp_n) / sum(denom_n),
    hosp_n = sum(hosp_n),
    denom_n = sum(denom_n),
    .groups = "drop"
  )


################################################################################
# 20 age group summary
################################################################################
age_breaks <- c(0, 1, 5, 10, 12, 18, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, Inf)
age_labels <- c("<1","1-4","5-9","10-11","12-17","18-19",
                "20-24","25-29","30-34","35-39","40-44","45-49",
                "50-54","55-59","60-64","65-69","70-74","75-79",
                "80-84","85+")


hosp_20 <- hosp %>%
  mutate(
    age_label = cut(
      age_group, # 
      breaks = age_breaks,
      labels = age_labels,
      right = FALSE, # [0,1), [1,5) ...
      include.lowest = TRUE
    )
  ) %>%
  filter(!is.na(age_label)) %>%
  group_by(age_label) %>%
  summarise(
    hosp_n = sum(n[hosp == TRUE], na.rm = TRUE),
    nonhosp_n = sum(n[hosp == FALSE], na.rm = TRUE),
    total_n = hosp_n + nonhosp_n,
    p_hosp = hosp_n / total_n,
    .groups = "drop"
  )

death_20 <- death %>%
  mutate(
    age_label = cut(
      age_group, # 
      breaks = age_breaks,
      labels = age_labels,
      right = FALSE, # [0,1), [1,5) ...
      include.lowest = TRUE
    )
  ) %>%
  filter(!is.na(age_label)) %>%
  group_by(age_label) %>%
  summarise(
    death_n = sum(n[deaht == TRUE], na.rm = TRUE),
    nondeath_n = sum(n[deaht == FALSE], na.rm = TRUE),
    total_n = death_n + nondeath_n,
    p_death = death_n / total_n,
    .groups = "drop"
  )

hosp <- hosp_20$p_hosp
fatal <- death_20$p_death

save(hosp, fatal, nh_fatal, file = "00_Data/0_2_Processed/chikv_fatal_hosp_rate.RData")
