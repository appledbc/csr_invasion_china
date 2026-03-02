# -------------------------------------------------------------------------
# @Description: Radar plot and bar plot visualization
# -------------------------------------------------------------------------
# Clear memory
cat("\014")
rm(list = ls())
gc()

# Load packges
library(broom)
library(car)
library(ggpubr)
library(ggraph)
library(gridExtra)
library(igraph)
library(lmtest)
library(patchwork)
library(purrr)
library(readr)
library(stringr)
library(tidygraph)
library(tidyverse)

# Set working directory (Please change this to your own data path)
# Set file path
file_path <- "D:/CSRdata"

# Check if file exists
if (!file.exists(file_path)) {
  file.create(file_path)
  cat("File created:", file_path, "\n")
  cat("Note: Please copy your code and data file to this location before running\n")
} else {
  cat("File exists:", file_path, "\n")

  # Check if file is empty
  if (file.info(file_path)$size == 0) {
    cat("Warning: File is empty, please add code and data file first\n")
  } else {
    cat("File is ready to use\n")
  }
}

setwd(file_path)
getwd()

# Import data (Please adjust the paths to point to your local data files)
Matti.gdp <- read_csv("./data/region_mean_Matti.gdp.csv")
Matti.hdi <- read_csv("./data/region_mean_Matti.hdi.csv")
chelsa.bio <- read_csv("./data/region_mean_chelsa.bio1_12.csv")

region.data <- Matti.gdp %>%
  left_join(Matti.hdi, by = "Geographic_region7") %>%
  left_join(chelsa.bio, by = "Geographic_region7") %>%
  arrange(desc(bio1))

region.data

# Extract data
Geo.regions <- c(
    "Northeast",
    "Northwest",
    "North",
    "Southwest",
    "Central",
    "East",
    "South")

# New order: N, NE, E, S, SW, C, NW (N at top)
Geo.regions <- c("North", "Northeast", "East", "South", "Southwest", "Central", "Northwest")

# ==========================================================================
# Load invasion data
# ==========================================================================
model.ias <- vector()

# Loop through regions
for (i in  1:length(Geo.regions)){

  # Extract single region
  region.i <- Geo.regions[i]

  print(paste(i, region.i, sep = "|"))

  # Load region-specific file
  mit02_file_name <- paste0("./results/phylo.sem.model_ias.tidy02_", region.i, ".csv")

  model_ias <- read.csv(file = mit02_file_name)

  # Preprocess data
  model_ias <- model_ias %>%
    filter(effect == "fixed" & !str_detect(term, "Intercept")) %>%
    mutate(region = region.i) %>%
    mutate(term = str_replace_all(term, "[.]", ""),
           response = str_replace_all(response, "invasionstatus", "invasion")) %>%
    mutate(Estimate = round(Estimate, 2),
           std.error = round(std.error, 3),
           conf.low = round(conf.low, 3),
           conf.high = round(conf.high, 3),
           Post.Prob = round(Post.Prob, 4)) %>%
    mutate(Star = case_when(Post.Prob > 0.999 ~ "***",
                          Post.Prob > 0.99 ~ "**",
                          Post.Prob > 0.95 ~ " * ",
                          Post.Prob < 0.95 ~ "",
                          )) %>%
    mutate(color = case_when(Estimate > 0 & Post.Prob > 0.95 ~ "red",
                           Estimate < 0 & Post.Prob > 0.95 ~ "blue",
                           Post.Prob < 0.95 ~ "grey",
                          ))

  # Select data
  model_ias01 <- model_ias %>%
    transmute(region = region,
              model = model,
              from = term,
              to = response,
              Estimate = Estimate,
              std.error = std.error,
              conf.low = conf.low,
              conf.high = conf.high,
              Post.Prob = Post.Prob,
              color = color,
              Star = Star,
              label = paste0(Estimate, Star)
              )

  # Combine ias data
  model.ias <- rbind(model.ias, model_ias01)

}

model.ias01 <- model.ias %>%
  mutate(from = str_replace_all(from, "[csr][_]scorescaled", "scorescaled")) %>%
  mutate(status = "ias")

# ==========================================================================
# Load naturalization data
# ==========================================================================
model.nat <- vector()

# Loop through regions
for (i in  1:length(Geo.regions)){

  # Extract single region
  region.i <- Geo.regions[i]

  print(paste(i, region.i, sep = "|"))

  # Load region-specific file
  mit02_file_name <- paste0("./results/phylo.sem.model_nat.tidy02_", region.i, ".csv")

  model_nat <- read.csv(file = mit02_file_name)

  # Preprocess data
  model_nat <- model_nat %>%
    filter(effect == "fixed" & !str_detect(term, "Intercept")) %>%
    mutate(region = region.i) %>%
    mutate(term = str_replace_all(term, "[.]", ""),
           response = str_replace_all(response, "invasionstatus", "invasion")) %>%
    mutate(Estimate = round(Estimate, 2),
           std.error = round(std.error, 3),
           conf.low = round(conf.low, 3),
           conf.high = round(conf.high, 3),
           Post.Prob = round(Post.Prob, 4)) %>%
    mutate(Star = case_when(Post.Prob > 0.999 ~ "***",
                          Post.Prob > 0.99 ~ "**",
                          Post.Prob > 0.95 ~ "*",
                          Post.Prob < 0.95 ~ "",
                          )) %>%
    mutate(color = case_when(Estimate > 0 & Post.Prob > 0.95 ~ "red",
                           Estimate < 0 & Post.Prob > 0.95 ~ "blue",
                           Post.Prob < 0.95 ~ "grey",
                          ))

  # Select data
  model_nat01 <- model_nat %>%
    transmute(region = region,
              model = model,
              from = term,
              to = response,
              Estimate = Estimate,
              std.error = std.error,
              conf.low = conf.low,
              conf.high = conf.high,
              Post.Prob = Post.Prob,
              color = color,
              Star = Star,
              label = paste0(Estimate, Star)
              )

  # Combine nat data
  model.nat <- rbind(model.nat, model_nat01)

}

model.nat01 <- model.nat %>%
  mutate(from = str_replace_all(from, "[csr][_]scorescaled", "scorescaled")) %>%
  mutate(status = "nat")

# Combine naturalization and invasion data
model.ias_nat01 <- rbind(model.ias01, model.nat01)
model.ias_nat02 <- model.ias_nat01 %>%
  left_join(region.data, by = c("region" = "Geographic_region7")) %>%
  mutate(log.gdp = log(gdp/1e+09),
         log.bio12 = log(bio12))

str(model.ias_nat02)

# ==========================================================================
# Group and apply regression models
# ==========================================================================
regression_models <- model.ias_nat02 %>%
  group_by(status, model, from, to) %>%
  nest() %>%
  mutate(
    models = map(data, ~ list(
      log.gdp = lm(Estimate ~ log.gdp, data = .x),
      hdi = lm(Estimate ~ hdi, data = .x),
      bio1 = lm(Estimate ~ bio1, data = .x),
      log.bio12 = lm(Estimate ~ log.bio12, data = .x)
    )),
    models = set_names(models, paste(status, model, from, to, sep = "_")),
    bp_test = map(models, ~ map(.x, bptest) %>% map(tidy)),
    shapiro_test = map(models, ~ map(.x, ~ shapiro.test(resid(.x))) %>% map(tidy)),
    model_summary = map(models, ~ map(.x, tidy)),
    model_summary = set_names(model_summary, paste(status, model, from, to, sep = "_"))
  )


# Extract and organize regression results
regression_results <- regression_models %>%
  mutate(coefs = map(model_summary, bind_rows, .id = 'predictor')) %>%
  select(status, model, from, to, coefs) %>%
  unnest(coefs)

# Save results to CSV
write_csv(regression_results, "./results/regression_results.csv")

# ==========================================================================
# Plotting
# ==========================================================================
# Define plotting function
plot_func <- function(data, group_vars){

  # Set title
  p.title <- paste(group_vars$status, group_vars$model, group_vars$from, group_vars$to, sep = " ~ ")

  # Extract unique color values
  color_values <- unique(data$color)

  # Standardize region and color
  data <- data %>%
    mutate(region = factor(region, levels = Geo.regions)) %>%
    mutate(color = factor(color, levels = color_values))

  # Create plot
  p <- ggplot(data, aes(x = region, y = Estimate)) +
      geom_bar(stat = "identity", aes(color = color, fill = color), size = 0) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      theme_classic() +
      theme(plot.title = element_text(family = "serif", colour = "black", size = 14),
            axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
            axis.text = element_text(family = "serif", colour = "black", size = 14),
            axis.text.x = element_text(angle = 45, vjust = 0.5),
            axis.title = element_text(family = "serif", colour = "black", size = 14),
            legend.title = element_text(family = "serif", colour = "black", size = 14),
            legend.text = element_text(family = "serif", colour = "black", size = 14),
            legend.position = "none"
      ) +
      scale_color_manual(values = color_values) +
      scale_fill_manual(values = color_values) +
      labs(title = "", x = "", y = "Estimate") +
      scale_y_continuous(limits = c(-1.5, 3), breaks = seq(-1,3, by = 1))

  return(p)

}

# ==========================================================================
# Batch generate data
# ==========================================================================
# Group data
grouped_data <- model.ias_nat01 %>%
  group_by(status, model, from, to)

# Extract group variables
group_vars <- grouped_data %>%
  group_keys()

# Batch generate new names
new_names <- map_chr(seq_len(nrow(group_vars)), ~ paste(group_vars$status[.x], group_vars$model[.x], group_vars$from[.x], group_vars$to[.x], sep = "_"))

# Batch plotting function
plot.ias_nat <- grouped_data %>%
  group_map(plot_func, .keep = TRUE)

# Rename elements in plot.ias_nat list
plot.ias_nat <- setNames(plot.ias_nat, new_names)

library(gridExtra)
do.call("grid.arrange", c(plot.ias_nat, ncol = 6))

# ==========================================================================
# Batch output plots by group
# ==========================================================================
# Define path lists and labels
path.lists <- list(
    score.to.mrt = c(
        "nat_c_model_nat_scorescaled_mrtscaled",
        "nat_s_model_nat_scorescaled_mrtscaled",
        "nat_r_model_nat_scorescaled_mrtscaled",
        "ias_c_model_ias_scorescaled_mrtscaled",
        "ias_s_model_ias_scorescaled_mrtscaled",
        "ias_r_model_ias_scorescaled_mrtscaled"),
    score.to.ias.nat = c(
        "nat_c_model_nat_scorescaled_invasion",
        "nat_s_model_nat_scorescaled_invasion",
        "nat_r_model_nat_scorescaled_invasion",
        "ias_c_model_ias_scorescaled_invasion",
        "ias_s_model_ias_scorescaled_invasion",
        "ias_r_model_ias_scorescaled_invasion"),
    score.to.eco = c(
        "nat_c_model_nat_scorescaled_ecousescaled",
        "nat_s_model_nat_scorescaled_ecousescaled",
        "nat_r_model_nat_scorescaled_ecousescaled",
        "ias_c_model_ias_scorescaled_ecousescaled",
        "ias_s_model_ias_scorescaled_ecousescaled",
        "ias_r_model_ias_scorescaled_ecousescaled"),
    eco.to.mrt = c(
        "nat_c_model_nat_ecousescaled_mrtscaled",
        "nat_s_model_nat_ecousescaled_mrtscaled",
        "nat_r_model_nat_ecousescaled_mrtscaled",
        "ias_c_model_ias_ecousescaled_mrtscaled",
        "ias_s_model_ias_ecousescaled_mrtscaled",
        "ias_r_model_ias_ecousescaled_mrtscaled"),
    mrt.to.ias.nat = c(
        "nat_c_model_nat_mrtscaled_invasion",
        "nat_s_model_nat_mrtscaled_invasion",
        "nat_r_model_nat_mrtscaled_invasion",
        "ias_c_model_ias_mrtscaled_invasion",
        "ias_s_model_ias_mrtscaled_invasion",
        "ias_r_model_ias_mrtscaled_invasion"),
    eco.to.ias.nat = c(
        "nat_c_model_nat_ecousescaled_invasion",
        "nat_s_model_nat_ecousescaled_invasion",
        "nat_r_model_nat_ecousescaled_invasion",
        "ias_c_model_ias_ecousescaled_invasion",
        "ias_s_model_ias_ecousescaled_invasion",
        "ias_r_model_ias_ecousescaled_invasion")
    )

path.labels <- list(
    score.to.mrt = c(
        "Naturalizaton: C-score -> MRT",
        "Naturalizaton: S-score -> MRT",
        "Naturalizaton: R-score -> MRT",
        "Invasion: C-score -> MRT",
        "Invasion: S-score -> MRT",
        "Invasion: R-score -> MRT"),
    score.to.ias.nat = c(
        "C-score -> Naturalization",
        "S-score -> Naturalization",
        "R-score -> Naturalization",
        "C-score -> Invasion",
        "S-score -> Invasion",
        "R-score -> Invasion"),
    score.to.eco = c(
        "Naturalizaton: C-score -> Eco.use",
        "Naturalizaton: S-score -> Eco.use",
        "Naturalizaton: R-score -> Eco.use",
        "Invasion: C-score -> Eco.use",
        "Invasion: S-score -> Eco.use",
        "Invasion: R-score -> Eco.use"),
    eco.to.mrt = c(
        "Naturalizaton: Eco.use -> MRT",
        "Naturalizaton: Eco.use -> MRT",
        "Naturalizaton: Eco.use -> MRT",
        "Invasion: Eco.use -> MRT",
        "Invasion: Eco.use -> MRT",
        "Invasion: Eco.use -> MRT"),
    mrt.to.ias.nat = c(
        "C-score: MRT -> Naturalization",
        "S-score: MRT -> Naturalization",
        "R-score: MRT -> Naturalization",
        "C-score: MRT -> Invasion",
        "S-score: MRT -> Invasion",
        "R-score: MRT -> Invasion"),
    eco.to.ias.nat = c(
        "C-score: Eco.use -> Naturalization",
        "S-score: Eco.use -> Naturalization",
        "R-score: Eco.use -> Naturalization",
        "C-score: Eco.use -> Invasion",
        "S-score: Eco.use -> Invasion",
        "R-score: Eco.use -> Invasion")
    )

# ==========================================================================
# 2. Ensure region is ordered factor in data preprocessing
# ==========================================================================
# Get Set2 colors
set2_colors <- RColorBrewer::brewer.pal(7, "Set2")

# Create named vector with abbreviations
custom_colors <- c("N" = set2_colors[1],
                   "NE" = set2_colors[2],
                   "E" = set2_colors[5],
                   "S" = set2_colors[6],
                   "SW" = set2_colors[7],
                   "C" = set2_colors[4],
                   "NW" = set2_colors[3])

# Calculate angle positions for labels
base_angle <- 360/7
custom_angles <- c("N" = 0,
                   "NE" = -base_angle*1,
                   "E" = -base_angle*2,
                   "S" = -base_angle*3,
                   "SW" = -base_angle*4,
                   "C" = -base_angle*5,
                   "NW" = -base_angle*6)

# ==========================================================================
# Part 2: Create and save radar plots
# ==========================================================================
# Create radar plots
radar_plot_list <- list()

for (path.list in names(path.lists)) {
  radar_plots <- list()

  folder_names <- paste0(file, "/results/", path.list)

  # Create folders if not exist
  for (folder in folder_names) {
    if (!dir.exists(folder)) {
      dir.create(folder, recursive = TRUE)
    }
  }

  for (i in seq_along(path.lists[[path.list]])) {
    p <- path.lists[[path.list]][i]
    l <- path.labels[[path.list]][i]

    # Create radar plot
    base_plot <- plot.ias_nat[[p]]
    plot_data <- base_plot$data
    plot_data$wrapped_region <- reorder(str_wrap(plot_data$region, 7), rep(1, length(plot_data$region)))

    radar_plots[[i]] <- base_plot +
      geom_hline(
        aes(yintercept = y),
        data.frame(y = seq(-0.5, 0.5, by = 0.5)),
        color = "lightgrey") +
      geom_segment(
        aes(
          x = region,
          y = -0.5,
          xend = region,
          yend = 0.5
        ),
        linetype = "dashed",
        color = "gray"
      ) +
      geom_bar(data = plot_data, aes(color = color, fill = color), stat = "identity") +
      coord_polar(start = -pi/7) +
      geom_hline(
        aes(yintercept = y),
        data.frame(y = 0),
        color = "black",
        size = 0.1) +
      annotate(
        x = 0.5,
        y = seq(-0.5, 0.5, by = 0.5),
        label = as.character(seq(-0.5, 0.5, by = 0.5)),
        geom = "text",
        color = "gray",
        size = 5,
        family = "serif"
      ) +
      ylim(-0.85, 0.85) +
      theme(
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(family = "serif", angle = unname(custom_angles), face = "bold", vjust = 0.5, hjust = 0.5, size = 16),
        legend.position = "none"
      ) +
      scale_x_discrete(labels = c("N", "NE", "E", "S", "SW", "C", "NW"))

    # Save plot with transparent background
    ggsave(filename = paste0(folder_names, "/20250401.radar_plots.v4.", p, ".png"),
           plot = radar_plots[[i]],
           width = 800 / 300,
           height = 800 / 300,
           dpi = 300,
           pointsize = 12)

  }

  # Save combined radar plots
  radar_plot_list[[path.list]] <- wrap_plots(radar_plots, ncol = 3) + plot_annotation(tag_levels = 'A')
  ggexport(radar_plot_list[[path.list]], filename = paste0("./results/20250401.radar_plots.v4.", path.list, ".png"),
    width = 2400,
    height = 1600,
    pointsize = 12,
    res = 300)
}
