# ----------------------------------------------------------
# Gantt Chart Generator using ganttrify
# ----------------------------------------------------------
# This script creates a Gantt chart using the ganttrify package
# A Gantt chart shows the timeline of tasks or activities

# Load required packages
# Install ganttrify if not already installed:
# install.packages("ganttrify")

library(ganttrify)
library(dplyr)

# ----------------------------------------------------------
# Example 1: Basic Gantt Chart
# ----------------------------------------------------------

# Create example data
# Required columns: wp (work package), activity, start_date, end_date
# Dates should be in format "YYYY-MM-DD"
tasks_basic <- data.frame(
  wp = c("WP1", "WP1", "WP2", "WP2", "WP3"),
  activity = c("Task A", "Task B", "Task C", "Task D", "Task E"),
  start_date = c("2024-01-01", "2024-01-15", "2024-02-01", "2024-02-10", "2024-03-01"),
  end_date = c("2024-01-20", "2024-02-10", "2024-02-28", "2024-03-15", "2024-03-30"),
  stringsAsFactors = FALSE
)

# Create Gantt chart
gantt_basic <- ganttrify(project = tasks_basic,
                         spots = NULL,
                         by_date = FALSE,
                         exact_date = TRUE,
                         project_start_date = NULL,
                         colour_palette = wesanderson::wes_palette("Darjeeling1"),
                         font_family = "sans",
                         mark_quarters = FALSE,
                         mark_years = FALSE,
                         size_wp = 6,
                         size_activity = 4,
                         size_text_relative = 1,
                         label_wp = NULL,
                         hide_wp = FALSE,
                         month_breaks = TRUE,
                         alpha_wp = 0.5,
                         alpha_activity = 0.8,
                         line_end = "round",
                         month_number_label = TRUE,
                         show_vertical_lines = TRUE,
                         axis_text_align = "right",
                         colour_stripe = "lightgray")

print(gantt_basic)

# Save plot
# ggsave("gantt_basic.pdf", plot = gantt_basic, width = 12, height = 6)


# ----------------------------------------------------------
# Example 2: Gantt Chart with Spots (Milestones/Deadlines)
# ----------------------------------------------------------

# Create example data with tasks
tasks_with_spots <- data.frame(
  wp = c("WP1", "WP1", "WP2", "WP2", "WP3"),
  activity = c("Planning", "Data Collection", "Analysis", "Writing", "Review"),
  start_date = c("2024-01-01", "2024-01-10", "2024-02-01", "2024-03-01", "2024-04-01"),
  end_date = c("2024-01-09", "2024-01-31", "2024-02-29", "2024-03-31", "2024-04-14"),
  stringsAsFactors = FALSE
)

# Create spots (milestones) data frame
# Required columns: spot (milestone name), spot_date, activity (must match activity in project)
spots <- data.frame(
  spot = c("Kick-off", "Data Ready", "First Draft", "Submission"),
  spot_date = c("2024-01-01", "2024-01-31", "2024-03-31", "2024-04-14"),
  activity = c("Planning", "Data Collection", "Writing", "Review"),
  stringsAsFactors = FALSE
)

# Create Gantt chart with spots
gantt_with_spots <- ganttrify(project = tasks_with_spots,
                              spots = spots,
                              by_date = FALSE,
                              exact_date = TRUE,
                              project_start_date = NULL,
                              colour_palette = wesanderson::wes_palette("Darjeeling1"),
                              font_family = "sans",
                              mark_quarters = FALSE,
                              mark_years = FALSE,
                              size_wp = 6,
                              size_activity = 4,
                              size_text_relative = 1,
                              hide_wp = FALSE,
                              month_breaks = TRUE,
                              alpha_wp = 0.5,
                              alpha_activity = 0.8)

print(gantt_with_spots)

# Save plot
# ggsave("gantt_with_spots.pdf", plot = gantt_with_spots, width = 12, height = 7)


# ----------------------------------------------------------
# Example 3: Gantt Chart with Quarters and Years
# ----------------------------------------------------------

tasks_quarters <- data.frame(
  wp = c("Phase 1", "Phase 1", "Phase 2", "Phase 2", "Phase 3"),
  activity = c("Task A", "Task B", "Task C", "Task D", "Task E"),
  start_date = c("2024-01-01", "2024-02-15", "2024-04-01", "2024-06-01", "2024-09-01"),
  end_date = c("2024-02-14", "2024-03-31", "2024-05-31", "2024-08-31", "2024-12-31"),
  stringsAsFactors = FALSE
)

# Create Gantt chart with quarter and year markers
gantt_quarters <- ganttrify(project = tasks_quarters,
                            spots = NULL,
                            by_date = FALSE,
                            exact_date = TRUE,
                            project_start_date = NULL,
                            colour_palette = wesanderson::wes_palette("Darjeeling1"),
                            font_family = "sans",
                            mark_quarters = TRUE,
                            mark_years = TRUE,
                            size_wp = 6,
                            size_activity = 4,
                            size_text_relative = 1,
                            hide_wp = FALSE,
                            month_breaks = TRUE,
                            alpha_wp = 0.5,
                            alpha_activity = 0.8)

print(gantt_quarters)

# Save plot
# ggsave("gantt_quarters.pdf", plot = gantt_quarters, width = 14, height = 7)


# ----------------------------------------------------------
# Example 4: Customized Gantt Chart
# ----------------------------------------------------------

# Create project data
project_data <- data.frame(
  wp = c("WP1", "WP1", "WP2", "WP2", "WP2", "WP3"),
  activity = c("Planning", "Preparation", "Data Collection", "Analysis", "Writing", "Submission"),
  start_date = c("2024-01-01", "2024-01-10", "2024-02-01", "2024-03-01", "2024-04-01", "2024-05-01"),
  end_date = c("2024-01-09", "2024-01-31", "2024-02-29", "2024-03-31", "2024-04-30", "2024-05-15"),
  stringsAsFactors = FALSE
)

# Create milestones
milestones <- data.frame(
  spot = c("Start", "Data Complete", "Analysis Complete", "Final Submission"),
  spot_date = c("2024-01-01", "2024-02-29", "2024-03-31", "2024-05-15"),
  activity = c("Planning", "Data Collection", "Analysis", "Submission"),
  stringsAsFactors = FALSE
)

# Create customized Gantt chart
gantt_custom <- ganttrify(project = project_data,
                          spots = milestones,
                          by_date = FALSE,
                          exact_date = TRUE,
                          project_start_date = "2024-01-01",
                          colour_palette = wesanderson::wes_palette("GrandBudapest1"),
                          font_family = "sans",
                          mark_quarters = TRUE,
                          mark_years = FALSE,
                          size_wp = 7,
                          size_activity = 5,
                          size_text_relative = 1.2,
                          label_wp = NULL,
                          hide_wp = FALSE,
                          month_breaks = TRUE,
                          alpha_wp = 0.6,
                          alpha_activity = 0.9,
                          line_end = "round",
                          month_number_label = TRUE,
                          show_vertical_lines = TRUE,
                          axis_text_align = "right",
                          colour_stripe = "lightblue")

print(gantt_custom)

# Save plot
# ggsave("gantt_custom.pdf", plot = gantt_custom, width = 14, height = 8)


# ----------------------------------------------------------
# Example 5: Simple Usage (Minimal Parameters)
# ----------------------------------------------------------

# Simple data frame
simple_project <- data.frame(
  wp = c("Phase 1", "Phase 2", "Phase 3"),
  activity = c("Planning & Setup", "Execution", "Finalization"),
  start_date = c("2024-01-01", "2024-02-01", "2024-03-01"),
  end_date = c("2024-01-31", "2024-02-29", "2024-03-31"),
  stringsAsFactors = FALSE
)

# Create simple Gantt chart (using defaults)
gantt_simple <- ganttrify(project = simple_project)

print(gantt_simple)

# Save plot
# ggsave("gantt_simple.pdf", plot = gantt_simple, width = 10, height = 5)


# ----------------------------------------------------------
# Notes:
# ----------------------------------------------------------
# 
# Required columns in project data frame:
# - wp: Work package (grouping variable)
# - activity: Activity/task name
# - start_date: Start date (format: "YYYY-MM-DD")
# - end_date: End date (format: "YYYY-MM-DD")
#
# Optional spots data frame (for milestones):
# - spot: Milestone name
# - spot_date: Milestone date (format: "YYYY-MM-DD")
# - activity: Activity name (must match activity in project)
#
# Note: You may need to install the 'wesanderson' package for color palettes:
# install.packages("wesanderson")
#
# Alternative: You can use custom color palettes or hex color codes
# Example: colour_palette = c("#3498db", "#e74c3c", "#2ecc71", "#f39c12")
