# -------------------------------------------------------------------------
# @Description: Dot plot with linear regression fitting
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
library(ggrepel)

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

print(region.data)

# ==========================================================================
# Extract data
# ==========================================================================
Geo.regions <- c(
    "Northeast",
    "Northwest",
    "North",
    "Southwest",
    "Central",
    "East",
    "South")

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
    filter(effect == "fixed" & term != "Intercept") %>%
    mutate(region = region.i) %>%
    mutate(term = str_replace_all(term, "[.]", ""),
           response = str_replace_all(response, "invasionstatus", "invasion")) %>%
    mutate(Estimate = round(Estimate, 2),
           Est.Error = round(Est.Error, 3),
           CI.Lower = round(CI.Lower, 3),
           CI.Upper = round(CI.Upper, 3),
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
              Est.Error = Est.Error,
              CI.Lower = CI.Lower,
              CI.Upper = CI.Upper,
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
    filter(effect == "fixed" & term != "Intercept") %>%
    mutate(region = region.i) %>%
    mutate(term = str_replace_all(term, "[.]", ""),
           response = str_replace_all(response, "invasionstatus", "invasion")) %>%
    mutate(Estimate = round(Estimate, 2),
           Est.Error = round(Est.Error, 3),
           CI.Lower = round(CI.Lower, 3),
           CI.Upper = round(CI.Upper, 3),
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
              Est.Error = Est.Error,
              CI.Lower = CI.Lower,
              CI.Upper = CI.Upper,
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
# Batch regression plots
# ==========================================================================
# Data processing
all_plots <- list()
regression_results <- data.frame()

# Get unique combinations
unique_combinations <- unique(model.ias_nat02[, c("status", "model", "from", "to")])

# Environment variables
env_vars <- c("bio1", "log.bio12", "log.gdp", "hdi")
env_var_labels <- c("MAT (°C)", "Log(MAP; mm)", "Log(GDP; billion $)", "HDI")
names(env_var_labels) <- env_vars

# Value ranges
env_x_min <- c(0, 6, 5, 0.55)
env_x_max <- c(25, 8, 8, 0.7)
y_limits <- c(-1, 1)
y_breaks <- seq(-1, 1, by = 0.5)

env_df <- data.frame(env_vars, env_var_labels, env_x_min, env_x_max)

# Region abbreviation mapping
region_abbr <- c(
  "Northeast" = "NE",
  "Northwest" = "NW",
  "North" = "N",
  "Southwest" = "SW",
  "Central" = "C",
  "East" = "E",
  "South" = "S"
)

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
        "Naturalizaton: C-score -> NEU",
        "Naturalizaton: S-score -> NEU",
        "Naturalizaton: R-score -> NEU",
        "Invasion: C-score -> NEU",
        "Invasion: S-score -> NEU",
        "Invasion: R-score -> NEU"),
    eco.to.mrt = c(
        "Naturalizaton: NEU -> MRT",
        "Naturalizaton: NEU -> MRT",
        "Naturalizaton: NEU -> MRT",
        "Invasion: NEU -> MRT",
        "Invasion: NEU -> MRT",
        "Invasion: NEU -> MRT"),
    mrt.to.ias.nat = c(
        "C-score: MRT -> Naturalization",
        "S-score: MRT -> Naturalization",
        "R-score: MRT -> Naturalization",
        "C-score: MRT -> Invasion",
        "S-score: MRT -> Invasion",
        "R-score: MRT -> Invasion"),
    eco.to.ias.nat = c(
        "C-score: NEU -> Naturalization",
        "S-score: NEU -> Naturalization",
        "R-score: NEU -> Naturalization",
        "C-score: NEU -> Invasion",
        "S-score: NEU -> Invasion",
        "R-score: NEU -> Invasion")
    )

path.nums <- list(
    score.to.mrt = c(
        "C score \u2192 MRT",
        "S score \u2192 MRT",
        "R score \u2192 MRT",
        "C score \u2192 MRT",
        "S score \u2192 MRT",
        "R score \u2192 MRT"),
    score.to.ias.nat = c(
        "C score \u2192 Natur.",
        "S score \u2192 Natur.",
        "R score \u2192 Natur.",
        "C score \u2192 Invas.",
        "S score \u2192 Invas.",
        "R score \u2192 Invas."),
    score.to.eco = c(
        "C score \u2192 NEU",
        "S score \u2192 NEU",
        "R score \u2192 NEU",
        "C score \u2192 NEU",
        "S score \u2192 NEU",
        "R score \u2192 NEU"),
    eco.to.mrt = c(
        "NEU \u2192 MRT",
        "NEU \u2192 MRT",
        "NEU \u2192 MRT",
        "NEU \u2192 MRT",
        "NEU \u2192 MRT",
        "NEU \u2192 MRT"),
    mrt.to.ias.nat = c(
        "MRT \u2192 Natur.",
        "MRT \u2192 Natur.",
        "MRT \u2192 Natur.",
        "MRT \u2192 Invas.",
        "MRT \u2192 Invas.",
        "MRT \u2192 Invas."),
    eco.to.ias.nat = c(
        "NEU \u2192 Natur.",
        "NEU \u2192 Natur.",
        "NEU \u2192 Natur.",
        "NEU \u2192 Invas.",
        "NEU \u2192 Invas.",
        "NEU \u2192 Invas.")
    )


# Convert lists to vectors
paths <- unlist(path.lists)
labels <- unlist(path.labels)
nums <- unlist(path.nums)
path_df <- data.frame(path = paths, label = labels, num = nums)

# Loop through each combination
for (i in 1:nrow(unique_combinations)) {

  current_status <- unique_combinations$status[i]
  current_model <- unique_combinations$model[i]
  current_from <- unique_combinations$from[i]
  current_to <- unique_combinations$to[i]

  plot_name <- paste(current_status, current_model, current_from, current_to, sep = "_")
  print(plot_name)

  # Filter data
  plot_data <- model.ias_nat02 %>%
    filter(status == current_status,
           model == current_model,
           from == current_from,
           to == current_to)

  # Create list for regression results
  env_plots <- list()

  # Regression for each environment variable
  for (env_var in env_vars) {
    # Linear regression
    lm_formula <- as.formula(paste("Estimate ~", env_var))
    lm_model <- lm(lm_formula, data = plot_data)

    # Collect regression results
    intercept <- round(coef(lm_model)[1], 3)
    slope <- round(coef(lm_model)[2], 3)
    r_squared <- round(summary(lm_model)$r.squared, 3)
    p_value <- anova(lm_model)$`Pr(>F)`[1]

    # Save regression results
    result_row <- data.frame(
      status = current_status,
      model = current_model,
      from = current_from,
      to = current_to,
      env_var = env_var,
      intercept = intercept,
      slope = slope,
      p_value = p_value,
      r_squared = r_squared
    )

    regression_results <- rbind(regression_results, result_row)

    # Create p-value label with italic formatting
    if(p_value < 0.001) {
      p_text <- "italic(P) < 0.001"
    } else {
      p_text <- sprintf("italic(P) == %.3f", round(p_value, 3))
    }

    # Linear regression equation
    equation <- sprintf("atop(y == %.3f * x %+.3f, italic(R)^2 == %.3f~';'~%s)",
                        slope, intercept, r_squared, p_text)

    # Special handling for very small slope values
    if(abs(slope) < 0.001) {
      slope <- round(coef(lm_model)[2], 4)

      equation <- sprintf("atop(y == %.4f * x %+.3f, italic(R)^2 == %.3f~';'~%s)",
                          slope, intercept, r_squared, p_text)
    }

    # Extract key information
    x_label <- env_df %>%
      filter(env_vars == env_var) %>%
      pull(env_var_labels)

    x_min <- env_df %>%
      filter(env_vars == env_var) %>%
      pull(env_x_min)

    x_max <- env_df %>%
      filter(env_vars == env_var) %>%
      pull(env_x_max)

    y_label <- path_df %>%
      filter(path == plot_name) %>%
      pull(num)

    # Create base plot - add regression line only if p < 0.05
    if(p_value < 0.05){
      plot <- ggplot(plot_data, aes_string(x = env_var, y = "Estimate")) +
        geom_point(color = "black", alpha = 0.5, size = 5) +
        geom_text_repel(
          aes(label = region_abbr),
          size = 4,
          family = "serif",
          force = 1,
          max.overlaps = 15,
          box.padding = 0.5,
          point.padding = 0.3,
          segment.color = "grey50",
          min.segment.length = 0
        ) +
        geom_smooth(method = "lm", formula = y ~ x,
                    color = "red", fill = "red", alpha = 0.5)

    }else{
      plot <- ggplot(plot_data, aes_string(x = env_var, y = "Estimate")) +
        geom_point(color = "black", alpha = 0.5, size = 5)
    }

    # Add annotations and styling to the plot
    plot <- plot +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
        axis.text = element_text(family = "serif", colour = "black", size = 14),
        axis.title = element_text(family = "serif", colour = "black", size = 14),
        legend.title = element_text(family = "serif", colour = "black", size = 14),
        legend.text = element_text(family = "serif", colour = "black", size = 14)
      ) +
      scale_y_continuous(limits = y_limits, breaks = y_breaks) +
      scale_x_continuous(limits = c(x_min, x_max)) +
      labs(x = x_label, y = paste0("Path coeff. of\n", y_label))

    # Store combined plot
    plot_name <- paste(current_status, current_model, current_from, current_to, sep = "_")
    all_plots[[plot_name]][[env_var]] <- plot

  }
}

# View results
head(regression_results)

write_csv(regression_results, "./results/regression_detailed_results.csv")

# ==========================================================================
# Combine regression plots
# ==========================================================================
# Get unique combinations
unique_combinations <- unique(model.ias_nat02[, c("status", "model", "from", "to")])
unique_combinations02 <- unique(model.ias_nat02[, c("status", "model")])

# Environment variables
env_vars <- c("bio1", "log.bio12", "hdi")

# Create patchwork for each status-model combination
for (i in 1:nrow(unique_combinations02)) {

  current_status <- unique_combinations02$status[i]
  current_model <- unique_combinations02$model[i]

  # Filter current combination data
  current_combinations <- unique_combinations %>%
    filter(status == current_status, model == current_model)

  # Initialize plot list
  plot_combined <- list()

  # Create row for each from-to combination
  for (j in 1:nrow(current_combinations)) {
    current_from <- current_combinations$from[j]
    current_to <- current_combinations$to[j]

    # Build plot_name
    plot_name <- paste(current_status, current_model, current_from, current_to, sep = "_")

    # Create plot for each environment variable
    for (env_var in env_vars) {

      # If plot exists
      plot <- all_plots[[plot_name]][[env_var]]

      if(env_var != env_vars[1]){plot = plot + labs(y = "")}
      if(j != nrow(current_combinations)){plot = plot + labs(x = "")}

      # Add to row list
      plot_combined <- c(plot_combined, list(plot))

    }
    }

    # Combine all rows
    final_plot <- wrap_plots(plot_combined, ncol = 3) + plot_annotation(tag_levels = 'A')

    # Save combined plot
    file_name <- paste(current_status, current_model, "combined.noGDP.png", sep = "_")

    ggexport(final_plot, filename = paste0("./results/", file_name),
             width = 3600,
             height = 4200,
             pointsize = 12,
             res = 300)
}
