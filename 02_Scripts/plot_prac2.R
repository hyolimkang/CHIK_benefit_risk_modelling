## ------------------------------------------------------------
## FINAL: Risk–Benefit plane where slope y/x = BRR
## x = excess adverse outcomes attributable to vaccination
## y = infection-related adverse outcomes averted by vaccination
## Background = log10(BRR) with midpoint at BRR=1 (y=x)
## ------------------------------------------------------------


# show BRR in [0.1, 10] most clearly (tune if needed)
log_lim <- 1  # -1..1 corresponds to BRR 0.1..10

# 3) Errorbar thickness helpers
y_range_global <- y_max_global - 0

# 1) Global ranges (from your representative summary)

x_max_global <- max(br_representative_benefit$x_hi, na.rm = TRUE) * 1.2
y_max_global <- max(br_representative_benefit$y_hi, na.rm = TRUE) * 1.2

cat("  X: [0,", round(x_max_global, 2), "]\n")
cat("  Y: [0,", round(y_max_global, 2), "]\n\n")

# 2) Background grid: log10(BRR) = log10(y/x)
eps_x <- x_max_global * 1e-6
bg_grid_base <- expand.grid(
  x = seq(eps_x, x_max_global, length.out = 250),
  y = seq(0, y_max_global, length.out = 250)
) %>%
  mutate(
    brr = y / x,
    log10_brr = log10(brr)
  )

facet_combinations <- br_representative_benefit %>%
  distinct(ar_category, days)

bg_grid_by_facet <- facet_combinations %>%
  rowwise() %>%
  mutate(grid = list(bg_grid_base %>% 
                       mutate(ar_category = ar_category,
                              days = days))) %>%
  ungroup() %>%
  select(grid) %>%
  tidyr::unnest(cols = grid)

actual_brr_range <- br_representative_benefit %>%
  mutate(brr = y_med / x_med) %>%
  summarise(
    min_brr = min(brr, na.rm = TRUE),
    max_brr = max(brr, na.rm = TRUE)
  )
cat("  Min BRR:", round(actual_brr_range$min_brr, 2), "\n")
cat("  Max BRR:", round(actual_brr_range$max_brr, 2), "\n")

log_min <- floor(log10(actual_brr_range$min_brr))
log_max <- ceiling(log10(actual_brr_range$max_brr))
log_breaks <- seq(log_min, log_max, by = 1)
brr_breaks <- 10^log_breaks

plot_risk_travel_improved <- 
  ggplot() +
  
  geom_raster(
    data = bg_grid_by_facet,
    aes(x = x, y = y, fill = log10_brr),
    alpha = 0.85,
    interpolate = TRUE
  ) +
  
  scale_fill_gradient2(
    name = "Benefit–risk ratio\n(log10 scale; BRR = averted / excess)",
    low = "#ca0020",      
    mid = "#f7f7f7",      #  BRR = 1
    high = "#0571b0",     #  BRR > 1
    midpoint = 0,         # log10(BRR) = 0 → BRR = 1
    limits = c(log_min, log_max),
    breaks = log_breaks,
    labels = paste0("BRR = ", brr_breaks),
    oob = scales::squish
  ) +
  
  # BRR = 1  (y = x)
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", linewidth = 0.7, color = "black"
  ) +
  
  # Y error bars
  geom_errorbar(
    data = br_representative_benefit,
    aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome),
    width = 0,  
    linewidth = 0.8,
    alpha = 0.8
  ) +
  
  # X error ars
  geom_errorbarh(
    data = br_representative_benefit,
    aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome),
    height = 0,  
    linewidth = 0.8,
    alpha = 0.8
  ) +
  
  # points
  geom_point(
    data = br_representative_benefit,
    aes(x = x_med, y = y_med, shape = AgeCat, color = outcome),
    size = 3, stroke = 0.9, fill = "white"
  ) +
  
  # Facet
  facet_grid(
    rows = vars(ar_category),
    cols = vars(days),
    scales = "fixed"  
  ) +
  
  scale_color_manual(
    values = c("SAE" = "#1B7F1B", "Death" = "#B8860B"),
    name = "Outcome type"
  ) +
  
  scale_shape_manual(
    values = c("18-64" = 22, "65+" = 24),
    name = "Age group"
  ) +
  
  coord_cartesian(
    xlim = c(0, x_max_global),
    ylim = c(0, y_max_global),
    expand = FALSE
  ) +
  guides(
    fill  = guide_colorbar(order = 1, title.position = "top", barwidth = 1.5, barheight = 12),
    color = guide_legend(order = 2, title.position = "top"),
    shape = guide_legend(order = 3, title.position = "top")
  ) +
  
  labs(
    x = "Excess adverse outcomes attributable to vaccination\n(per 10,000 vaccinated)",
    y = "Infection-related adverse outcomes averted by vaccination\n(per 10,000 vaccinated)",
    title = "Benefit-Risk Space: Vaccination outcomes across travel duration and attack rates"
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "right",
    legend.box = "vertical",
    legend.key.height = unit(0.8, "lines"),
    panel.spacing = unit(0.6, "lines"),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
  )

print(plot_risk_travel_improved)
