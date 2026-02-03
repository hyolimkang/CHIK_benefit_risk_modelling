ft <- flextable(ar_table_wide) %>%
  set_header_labels(
    outcome     = "Outcome",
    ar_category = "Attack Rate",
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
  merge_v(j = c("outcome", "ar_category")) %>%
  valign(j = c("outcome", "ar_category"), valign = "top") %>%
  border_remove()

doc <- read_docx() %>%
  body_add_par("Benefit–Risk Ratio (BRR) by Outcome, Attack Rate, Age Group, and Travel Duration", style = "heading 2") %>%
  body_add_flextable(ft)

print(doc, target = "06/Results/BRR_table.docx")