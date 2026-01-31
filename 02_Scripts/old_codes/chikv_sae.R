# input data 
death_plot <- overall %>%
  filter(death_chikv == TRUE)
death_plot <- death_plot %>% filter(age_group != "0-99")

ggplot(death_plot, aes(x = age_group, y = death_prop, group = 1)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(size = 2, color = "steelblue") +
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "grey40") +
  labs(
    x = "Age group",
    y = "Proportion of deaths among total N",
    title = "CHIKV CFR by age group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )

cases_by_age <- overall %>%                  
  filter(age_group != "0-99") %>%             
  group_by(age_group) %>%
  summarise(total_cases = sum(n), .groups = "drop")

hosp_by_age <- hosp_chikv %>%                   
  filter(age_group != "0-99") %>%         
  group_by(age_group) %>%
  summarise(hosp_n = sum(n), .groups = "drop")

hosp_rate <- hosp_by_age %>%
  left_join(cases_by_age, by = "age_group") %>%
  mutate(hosp_rate = hosp_n / total_cases)

hosp_rate <- hosp_rate %>% filter(age_group != "0-99")

ggplot(hosp_rate, aes(x = age_group, y = hosp_rate, group = 1)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(size = 2, color = "steelblue") +
  labs(
    x = "Age group",
    y = "Proportion of hospitalisation among total N",
    title = "CHIKV hospitalisation by age group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )


hosp_rate <- hosp_rate %>%
  mutate(age_cat = ifelse(age_group %in% c("[0,10)", "[10,20)", "[20,30)", "[30,40)", "[40,50)", "[50,60)"),
                          "<60", "60+")) %>%
  group_by(age_cat) %>%
  summarise(
    hosp_n_total   = sum(hosp_n),
    total_cases_sum= sum(total_cases),
    pooled_rate    = hosp_n_total / total_cases_sum
  )
