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

### ori
ft_brr <- flextable(brr_table_wide) %>%
  theme_booktabs() %>%
  autofit() %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = 1:2, align = "left", part = "all") %>%
  merge_v(j = "Outcome") %>%
  valign(j = "Outcome", valign = "top")

ft_brr

doc <- read_docx() %>%
  body_add_par(
    "Benefit–Risk Ratio (BRR) by Outcome, Age group, and VE",
    style = "heading 2"
  ) %>%
  body_add_flextable(ft_brr)

print(doc, target = "06_Results/BRR_table_ori.docx")


## ori by setting
### ori
brr_table_wide_setting2 <- brr_table_wide_setting2 %>%
  arrange(Outcome, Setting, `Age group`)

ft_brr <- flextable(brr_table_wide_setting2) %>%
  set_header_labels(
    Outcome = "Outcome",
    Setting = "Setting",
    `Age group` = "Age Group",
    `Disease blocking only` = "Disease blocking only",
    `Disease and infection blocking` = "Disease and infection blocking"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = 1:3, align = "left", part = "all") %>%
  merge_v(j = c("Outcome", "Setting")) %>%          
  valign(j = c("Outcome", "Setting"), valign = "top") %>%
  border_remove()

ft_brr
doc <- read_docx() %>%
  body_add_par(
    "Benefit–Risk Ratio (BRR) by Outcome, Age group, and VE",
    style = "heading 2"
  ) %>%
  body_add_flextable(ft_brr)

print(doc, target = "06_Results/BRR_table_ori_setting.docx")

## ori national
### ori
ft_brr <- flextable(brr_table_wide_true) %>%
  theme_booktabs() %>%
  autofit() %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "all") %>%
  flextable::merge_v(j = "Outcome") %>%
  flextable::valign(j = "Outcome", valign = "top")

ft_brr

doc <- read_docx() %>%
  body_add_par(
    "Benefit–Risk Ratio (BRR) by Outcome, Age group, and VE",
    style = "heading 2"
  ) %>%
  body_add_flextable(ft_brr)

print(doc, target = "06_Results/BRR_table_ori_national.docx")
