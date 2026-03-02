# -------------------------------------------------------------------------
# @Description: SEM diagram generation (national scale)
# -------------------------------------------------------------------------
# Clear memory
cat("\014")
rm(list = ls())
gc()

# Loading packages
library(readr)
library(stringr)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(purrr)
library(patchwork)
library(ggpubr)

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

# ==========================================================================
# Function
# ==========================================================================
sem_models <- function(data, model_name, response_var) {

  data <- data %>%
    filter(model == model_name) %>%
    mutate(color = factor(color))

  color <- levels(data$color)

  data.g <- as_tbl_graph(data) %>%
    activate(nodes) %>%
    mutate(name = case_when(
      name == "c_scorescaled" ~ "C-score",
      name == "s_scorescaled" ~ "S-score",
      name == "r_scorescaled" ~ "R-score",
      name == "ecousescaled" ~ "Economic\nUse",
      name == "mrtscaled" ~ "MRT",
      name == "invasion" ~ "Invasion\nincidence",
      name == "naturalization" ~ "Naturalization\nincidence"
    ))

  nodes <- data.g %>%
    activate(nodes) %>%
    as_tibble() %>%
    mutate(x = c(0, 1, -1, 0),
           y = c(0.5, 0, 0, -1))

  sem.plot <- ggraph(data.g, layout = nodes, x = x, y = y) +
      geom_edge_link(
        aes(
          label = label,
          color = factor(color),
          width = width,
          start_cap = circle(20, 'mm'),
          end_cap = circle(20, 'mm')
        ),
         angle_calc = "along",
         label_dodge = unit(8, "mm"),
         label_size = 4,
         arrow = arrow(length = unit(3, "mm"), type = "closed")) +
      geom_node_point(size = 0, color = "lightblue") +
      geom_node_text(aes(label = name),
                     size = 5,
                     color = 'black',
                     hjust = c(0.5, 0.8, -0.7, 0.5)) +
      scale_edge_color_manual(values = color) +
      theme_graph() +
      theme(plot.title = element_text(hjust = 0.5, size = 16),
            legend.position = "none"
      ) +
      labs(title = model_name)

  return(sem.plot)

}

# ==========================================================================
# Load invasion data
# ==========================================================================
model_ias <- read.csv("./results/phylo.sem.model_ias.tidy02.csv")
model_ias

model_ias <- model_ias %>%
  filter(effect == "fixed" & !str_detect(term, "Intercept")) %>%
  mutate(term = str_replace_all(term, "[.]", ""),
         response = str_replace_all(response, "invasionstatus", "invasion")) %>%
  mutate(estimate = round(estimate, 3),
         Post.Prob = round(Post.Prob, 4)) %>%
  mutate(Star = case_when(Post.Prob > 0.999 ~ "***",
                        Post.Prob > 0.99 ~ "**",
                        Post.Prob > 0.95 ~ "*",
                        Post.Prob < 0.95 ~ "",
                        )) %>%
  mutate(color = case_when(estimate > 0 & Post.Prob > 0.95 ~ "red",
                         estimate < 0 & Post.Prob > 0.95 ~ "blue",
                         Post.Prob < 0.95 ~ "grey",
                        ))

model_ias01 <- model_ias %>%
  transmute(from = term,
            to = response,
            model = model,
            estimate = estimate,
            Post.Prob = Post.Prob,
            width = abs(estimate) * 5,
            color = color,
            label = paste(estimate, Star)
            )

model_names <- c("c_model_ias", "s_model_ias", "r_model_ias")
response_vars <- c("C-score", "S-score", "R-score")

sem_models.ias <- purrr::map2(model_names, response_vars, ~sem_models(model_ias01, .x, .y))
sem_models.ias <- setNames(sem_models.ias, model_names)

sem_models.ias01c <- sem_models.ias$c_model_ias
sem_models.ias01s <- sem_models.ias$s_model_ias
sem_models.ias01r <- sem_models.ias$r_model_ias

# ==========================================================================
# Load naturalization data
# ==========================================================================
model_nat <- read.csv("./results/phylo.sem.model_nat.tidy02.csv")
model_nat

model_nat <- model_nat %>%
  filter(effect == "fixed" & !str_detect(term, "Intercept")) %>%
  mutate(term = str_replace_all(term, "[.]", ""),
         response = str_replace_all(response, "invasionstatus", "naturalization")) %>%
  mutate(estimate = round(estimate, 3),
         Post.Prob = round(Post.Prob, 4)) %>%
  mutate(Star = case_when(Post.Prob > 0.999 ~ "***",
                        Post.Prob > 0.99 ~ "**",
                        Post.Prob > 0.95 ~ "*",
                        Post.Prob < 0.95 ~ "",
                        )) %>%
  mutate(color = case_when(estimate > 0 & Post.Prob > 0.95 ~ "red",
                         estimate < 0 & Post.Prob > 0.95 ~ "blue",
                         Post.Prob < 0.95 ~ "grey",
                        ))

model_nat01 <- model_nat %>%
  transmute(from = term,
            to = response,
            model = model,
            estimate = estimate,
            Post.Prob = Post.Prob,
            width = abs(estimate) * 5,
            color = color,
            label = paste(estimate, Star)
            )

model_names <- c("c_model_nat", "s_model_nat", "r_model_nat")
response_vars <- c("C-score", "S-score", "R-score")

sem_models.nat <- purrr::map2(model_names, response_vars, ~sem_models(model_nat01, .x, .y))
sem_models.nat <- setNames(sem_models.nat, model_names)

sem_models.nat01c <- sem_models.nat$c_model_nat
sem_models.nat01s <- sem_models.nat$s_model_nat
sem_models.nat01r <- sem_models.nat$r_model_nat

# ==========================================================================
# Summarize data
# ==========================================================================
sem_models01 <-
  sem_models.nat01c +
  sem_models.nat01s +
  sem_models.nat01r +
  sem_models.ias01c +
  sem_models.ias01s +
  sem_models.ias01r +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")

sem_models01

ggexport(sem_models01, filename = "./results/sem_models01.png",
           width = 6000,
           height = 4000,
           pointsize = 12,
           res = 300)
