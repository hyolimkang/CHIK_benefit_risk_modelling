### traveller 
conflicted::conflict_prefer("align", "flextable")

ft <- flextable(ar_table_wide2) %>%
  set_header_labels(
    outcome     = "Outcome",
    ar_category = "Setting",
    age_group   = "Age Group",
    `7d`  = "7 days",
    `14d` = "14 days",
    `30d` = "30 days",
    `90d` = "90 days"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1:3, align = "left", part = "all") %>%
  bold(part = "header") %>%
  merge_v(j = c("outcome", "setting")) %>%
  valign(j = c("outcome", "setting"), valign = "top") %>%
  border_remove()

doc <- read_docx() %>%
  body_add_par("Benefit–Risk Ratio (BRR) by Outcome, Attack Rate, Age Group, and Travel Duration", style = "heading 2") %>%
  body_add_flextable(ft)

print(doc, target = "06_Results/BRR_table_travel.docx")

ar_summary_all_df <- as_tibble(ar_summary_all)

# final ver.
ar_table_long <- ar_summary_all_df %>%
  left_join(pr_ge1_tbl, by = c("outcome", "setting", "age_group", "days")) %>%
  mutate(
    benefit = sprintf("%.2f [%.2f–%.2f]", av_med, av_lo, av_hi),
    risk    = sprintf("%.2f [%.2f–%.2f]", ca_med, ca_lo, ca_hi),
    brr     = sprintf("%.2f [%.2f–%.2f]", brr_med, brr_lo, brr_hi),
    prob    = sprintf("%.1f%%", 100 * pr_brr_ge_1)
  ) %>%
  dplyr::select(outcome, setting, age_group, days, benefit, risk, brr, prob)

ft_final <- flextable(ar_table_long) %>%
  set_header_labels(
    outcome = "Outcome", setting = "Setting", age_group = "Age group",
    days = "Travel\nduration",
    benefit = "Benefit:\nOutcomes averted...",
    risk = "Risk:\nOutcomes attributable...",
    brr = "Benefit-risk ratio...",
    prob = "Probability\n(BRR > 1)\n(%)"
  ) %>%
  merge_v(j = c("outcome", "setting", "age_group")) %>%
  theme_booktabs() %>%
  autofit()

doc <- read_docx() %>%
  body_add_par("Benefit–Risk Ratio (BRR) by Outcome, Attack Rate, Age Group, and Travel Duration", style = "heading 2") %>%
  body_add_flextable(ft_final)

print(doc, target = "06_Results/BRR_table_travel.docx")


## ori by setting
### ori
part1 <- brr_table_wide_setting2 %>%
  dplyr::select(
    Outcome, Setting, `Age group`,
    Benefit = `Disease blocking only (Benefit)`,
    Risk    = `Disease blocking only (Risk)`,
    BRR     = `Disease blocking only (BRR)`,
    prob    = `Disease blocking only Pr(BRR>1)`
  ) %>%
  mutate(mechanism = "Disease blocking only")

# 2. Extract and rename columns for 'Disease and infection blocking'
part2 <- brr_table_wide_setting2 %>%
  dplyr::select(
    Outcome, Setting, `Age group`,
    Benefit = `Disease and infection blocking (Benefit)`,
    Risk    = `Disease and infection blocking (Risk)`,
    BRR     = `Disease and infection blocking (BRR)`,
    prob    = `Disease and infection blocking Pr(BRR>1)`
  ) %>%
  mutate(mechanism = "Disease and infection blocking")

# 3. Stack them into one final long-form table
brr_table_final_long <- bind_rows(part1, part2) %>%
  arrange(Outcome, Setting, `Age group`) %>%
  dplyr::select(Outcome, Setting, `Age group`, mechanism, Benefit, Risk, BRR, prob)

ft_brr <- flextable(brr_table_final_long) %>%
  # Set professional multi-line headers (English)
  set_header_labels(
    Outcome   = "Outcome",
    Setting   = "Setting",
    `Age group` = "Age group",
    mechanism = "Vaccine protection\nmechanism",
    Benefit   = "Benefit:\nOutcomes averted\n(per 10,000)",
    Risk      = "Risk:\nOutcomes attributable\n(per 10,000)",
    BRR       = "Benefit-risk ratio:\n(Prevented per 1 caused)",
    prob      = "Probability\n(BRR > 1)\n(%)"
  ) %>%
  # Styling
  theme_booktabs() %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = 1:4, align = "left", part = "all") %>% 
  
  # Vertical merging for Outcome, Setting, and Age group
  merge_v(j = c("Outcome", "Setting", "Age group")) %>%
  valign(j = c("Outcome", "Setting", "Age group"), valign = "top") %>%
  
  # Fine-tuning
  fontsize(size = 9, part = "all") %>%
  autofit()

ft_brr

doc <- read_docx() %>%
  body_add_par(
    "Benefit–Risk Ratio (BRR) by Outcome, Age group, and VE",
    style = "heading 2"
  ) %>%
  body_add_flextable(ft_brr)

print(doc, target = "06_Results/BRR_table_ori_setting.docx")

