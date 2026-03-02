# -------------------------------------------------------------------------
# @Description: brms phylogenetic SEM construction (regional scale)
# -------------------------------------------------------------------------

# Tutorial: https://xiangyun.rbind.io/2024/05/stan-interfaces/

# Clear memory
cat("\014")
rm(list = ls())
gc()

# Loading packages
library(ape)
library(brms)
library(broom.mixed)
library(cmdstanr)
library(dplyr)
library(ggpubr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(tidytree)
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

# Set CmdStan path
set_cmdstan_path("C:/Users/dbc/Documents/.cmdstan/cmdstan-2.34.1")

# Loading CSR data
checklist.csr <- read_csv("./data/checklist_accepted_matches_CSR.WithCVH.mrt.v6simple.20251023.csv", show_col_types = FALSE)
str(checklist.csr)
View(checklist.csr)

# 31 provinces and 7 economic regions
table(checklist.csr$Provience_pinyin)
table(checklist.csr$Geo.region)

# Remove Taiwan, Macao, Hongkong
checklist.csr <- checklist.csr %>% filter(!(Provience_pinyin %in% c("Taiwan", "Macao", "Hongkong")))
str(checklist.csr)

checklist.csr.CHN.mrt <- checklist.csr %>%
  group_by(taxon_name, accepted_plant_name_id, invasion_status, TPL_names, c_score, s_score, r_score, wcup_eco_use_num, wcup_eco_use_num_2out, Geo.region) %>%
  summarise(fr.China = min(Checked_year), mrt.China = 2023 - fr.China) %>%
  ungroup()

checklist.csr.CHN.mrt01a <- checklist.csr.CHN.mrt %>%
  mutate_at(vars(c_score:wcup_eco_use_num_2out), ~ as.numeric(.)) %>%
  mutate(tip_label = str_replace(taxon_name, " ", "_"))

str(checklist.csr.CHN.mrt01a)

Geo.regions <- unique(checklist.csr.CHN.mrt01a$Geo.region)
Geo.regions

# Loop through regions for regional analysis
for (i in  1:length(Geo.regions)){

  # Set region
  region.i <- Geo.regions[i]

  checklist.csr.CHN.mrt01 <- checklist.csr.CHN.mrt01a %>%
    filter(Geo.region == region.i)

  # ==========================================================================
  # Get casual vs naturalized data (1 vs 2)
  # ==========================================================================
  checklist.csr.CHN.mrt.12 <- checklist.csr.CHN.mrt01 %>%
    mutate(invasion_status = case_when(invasion_status == 1 ~ 0,
                                       invasion_status %in% c(2, 3) ~ 1))

  table(checklist.csr.CHN.mrt01$invasion_status)
  table(checklist.csr.CHN.mrt.12$invasion_status)

  # Scale variables
  checklist.csr.CHN.mrt.12 <- checklist.csr.CHN.mrt.12 %>%
    mutate(c_score.scaled = scale(c_score),
           s_score.scaled = scale(s_score),
           r_score.scaled = scale(r_score),
           eco.use.scaled = scale(wcup_eco_use_num_2out),
           mrt.scaled = scale(mrt.China))

  # ==========================================================================
  # Load phylogenetic tree
  # ==========================================================================
  graft_status <- read_csv("./data/checklist_phylo.tree_graft_status_all.csr.sp.csv", show_col_types = FALSE)

  checklist.csr.CHN.mrt.12a <- checklist.csr.CHN.mrt.12 %>%
    left_join(graft_status, by = "tip_label")

  checklist.csr.CHN.mrt.12b <- checklist.csr.CHN.mrt.12a %>%
    arrange(tip_order) %>%
    tibble::column_to_rownames("tip_label")

  # Load phylogenetic database
  phylo_tree <- read.tree("./data/checklist_phylo.tree_all.csr.sp.tre")

  is.rooted(phylo_tree)

  # Important step
  sub.tree.12b <- geiger::treedata(phylo_tree, checklist.csr.CHN.mrt.12b, sort = TRUE, warnings = F)$phy
  is.rooted(sub.tree.12b)

  sub.tree.12b <- keep.tip(sub.tree.12b, checklist.csr.CHN.mrt.12b$species) %>% unroot()
  is.rooted(sub.tree.12b)

  # ==========================================================================
  # Build phylogenetic brms models
  # ==========================================================================
  num_cpu <- parallel::detectCores(logical = FALSE)
  num_cpu_logical <- parallel::detectCores(logical = TRUE)

  # Set matrix
  phy.m.corr.12b <- ape::vcv(sub.tree.12b, corr = TRUE)

  # Define formulas
  c_formula_eco.use <- brms::bf(eco.use.scaled ~ c_score.scaled + (1|gr(species, cov = phy.m.corr.12b)), family = gaussian())
  c_formula_mrt <- brms::bf(mrt.scaled ~ eco.use.scaled + c_score.scaled + (1|gr(species, cov = phy.m.corr.12b)), family = gaussian())
  c_formula_nat <- brms::bf(invasion_status ~ mrt.scaled + c_score.scaled + eco.use.scaled + (1|gr(species, cov = phy.m.corr.12b)), family = bernoulli(link = 'logit'))

  s_formula_eco.use <- brms::bf(eco.use.scaled ~ s_score.scaled + (1|gr(species, cov = phy.m.corr.12b)), family = gaussian())
  s_formula_mrt <- brms::bf(mrt.scaled ~ eco.use.scaled + s_score.scaled + (1|gr(species, cov = phy.m.corr.12b)), family = gaussian())
  s_formula_nat <- brms::bf(invasion_status ~ mrt.scaled + s_score.scaled + eco.use.scaled + (1|gr(species, cov = phy.m.corr.12b)), family = bernoulli(link = 'logit'))

  r_formula_eco.use <- brms::bf(eco.use.scaled ~ r_score.scaled + (1|gr(species, cov = phy.m.corr.12b)), family = gaussian())
  r_formula_mrt <- brms::bf(mrt.scaled ~ eco.use.scaled + r_score.scaled + (1|gr(species, cov = phy.m.corr.12b)), family = gaussian())
  r_formula_nat <- brms::bf(invasion_status ~ mrt.scaled + r_score.scaled + eco.use.scaled + (1|gr(species, cov = phy.m.corr.12b)), family = bernoulli(link = 'logit'))

  # Fit model
  c_model_nat <- brm(c_formula_eco.use + c_formula_mrt + c_formula_nat + set_rescor(FALSE),
    data = checklist.csr.CHN.mrt.12b,
    data2 = list(phy.m.corr.12b = phy.m.corr.12b),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    cores = 4,
    backend = "cmdstanr",
    threads = threading(7),
    iter = 2000,
    seed = 2024
  )

  # Fit model
  s_model_nat <- brm(s_formula_eco.use + s_formula_mrt + s_formula_nat + set_rescor(FALSE),
    data = checklist.csr.CHN.mrt.12b,
    data2 = list(phy.m.corr.12b = phy.m.corr.12b),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    cores = 4,
    backend = "cmdstanr",
    threads = threading(7),
    iter = 2000,
    seed = 2024
  )

  # Fit model
  r_model_nat <- brm(r_formula_eco.use + r_formula_mrt + r_formula_nat + set_rescor(FALSE),
    data = checklist.csr.CHN.mrt.12b,
    data2 = list(phy.m.corr.12b = phy.m.corr.12b),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    cores = 4,
    backend = "cmdstanr",
    threads = threading(7),
    iter = 2000,
    seed = 2024
  )

  # ==========================================================================
  # Get naturalized vs invasive data (2 vs 3)
  # ==========================================================================
  checklist.csr.CHN.mrt.23 <- checklist.csr.CHN.mrt01 %>%
    filter(invasion_status %in% c(2, 3)) %>%
    mutate(invasion_status = case_when(invasion_status == 2 ~ 0,
                                       invasion_status == 3 ~ 1))

  table(checklist.csr.CHN.mrt01$invasion_status)
  table(checklist.csr.CHN.mrt.23$invasion_status)

  # Scale variables
  checklist.csr.CHN.mrt.23 <- checklist.csr.CHN.mrt.23 %>%
    mutate(c_score.scaled = scale(c_score),
           s_score.scaled = scale(s_score),
           r_score.scaled = scale(r_score),
           eco.use.scaled = scale(wcup_eco_use_num_2out),
           mrt.scaled = scale(mrt.China))

  # ==========================================================================
  # Load phylogenetic tree
  # ==========================================================================
  graft_status <- read_csv("./data/checklist_phylo.tree_graft_status_all.csr.sp.csv", show_col_types = FALSE)

  checklist.csr.CHN.mrt.23a <- checklist.csr.CHN.mrt.23 %>%
    left_join(graft_status, by = "tip_label")

  checklist.csr.CHN.mrt.23b <- checklist.csr.CHN.mrt.23a %>%
    arrange(tip_order) %>%
    tibble::column_to_rownames("tip_label")

  # Load phylogenetic database
  phylo_tree <- read.tree("./data/checklist_phylo.tree_all.csr.sp.tre")

  is.rooted(phylo_tree)

  sub.tree.23b <- geiger::treedata(phylo_tree, checklist.csr.CHN.mrt.23b, sort = TRUE, warnings = F)$phy
  is.rooted(sub.tree.23b)

  sub.tree.23b <- keep.tip(sub.tree.23b, checklist.csr.CHN.mrt.23b$species) %>% unroot()
  is.rooted(sub.tree.23b)

  # ==========================================================================
  # Build brms models for invasion
  # ==========================================================================
  phy.m.corr.23b <- ape::vcv(sub.tree.23b, corr = TRUE)

  # Define formulas
  c_formula_eco.use <- brms::bf(eco.use.scaled ~ c_score.scaled + (1|gr(species, cov = phy.m.corr.23b)), family = gaussian())
  c_formula_mrt <- brms::bf(mrt.scaled ~ eco.use.scaled + c_score.scaled + (1|gr(species, cov = phy.m.corr.23b)), family = gaussian())
  c_formula_ias <- brms::bf(invasion_status ~ mrt.scaled + c_score.scaled + eco.use.scaled + (1|gr(species, cov = phy.m.corr.23b)), family = bernoulli(link = 'logit'))

  s_formula_eco.use <- brms::bf(eco.use.scaled ~ s_score.scaled + (1|gr(species, cov = phy.m.corr.23b)), family = gaussian())
  s_formula_mrt <- brms::bf(mrt.scaled ~ eco.use.scaled + s_score.scaled + (1|gr(species, cov = phy.m.corr.23b)), family = gaussian())
  s_formula_ias <- brms::bf(invasion_status ~ mrt.scaled + s_score.scaled + eco.use.scaled + (1|gr(species, cov = phy.m.corr.23b)), family = bernoulli(link = 'logit'))

  r_formula_eco.use <- brms::bf(eco.use.scaled ~ r_score.scaled + (1|gr(species, cov = phy.m.corr.23b)), family = gaussian())
  r_formula_mrt <- brms::bf(mrt.scaled ~ eco.use.scaled + r_score.scaled + (1|gr(species, cov = phy.m.corr.23b)), family = gaussian())
  r_formula_ias <- brms::bf(invasion_status ~ mrt.scaled + r_score.scaled + eco.use.scaled + (1|gr(species, cov = phy.m.corr.23b)), family = bernoulli(link = 'logit'))

  # Fit model
  c_model_ias <- brm(c_formula_eco.use + c_formula_mrt + c_formula_ias + set_rescor(FALSE),
    data = checklist.csr.CHN.mrt.23b,
    data2 = list(phy.m.corr.23b = phy.m.corr.23b),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    cores = 4,
    backend = "cmdstanr",
    threads = threading(7),
    iter = 2000,
    seed = 2024
  )

  # Fit model
  s_model_ias <- brm(s_formula_eco.use + s_formula_mrt + s_formula_ias + set_rescor(FALSE),
    data = checklist.csr.CHN.mrt.23b,
    data2 = list(phy.m.corr.23b = phy.m.corr.23b),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    cores = 4,
    backend = "cmdstanr",
    threads = threading(7),
    iter = 2000,
    seed = 2024
  )

  # Fit model
  r_model_ias <- brm(r_formula_eco.use + r_formula_mrt + r_formula_ias + set_rescor(FALSE),
    data = checklist.csr.CHN.mrt.23b,
    data2 = list(phy.m.corr.23b = phy.m.corr.23b),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    cores = 4,
    backend = "cmdstanr",
    threads = threading(7),
    iter = 2000,
    seed = 2024
  )

  # Save brm results
  brm.file_name <- paste0("./results/csr.brm_phylo.sem_models_", region.i, ".Rdata")

  save(c_model_nat,
       s_model_nat,
       r_model_nat,
       c_model_ias,
       s_model_ias,
       r_model_ias,
       file = brm.file_name)

  # ==========================================================================
  # Clean data
  # ==========================================================================
  c_model_nat.tidy <- c_model_nat %>%
    tidy() %>%
    mutate(
      dir = ifelse(estimate > 0, " > 0", " < 0"),
      term = str_replace_all(term, "[()]", ""),
      hyp = paste0(response, "_", term, dir),
      Hypothesis = paste0("(", response, "_", term, ")", dir)
    )

  c_model_nat.tidy01 <- c_model_nat.tidy %>%
    filter(effect == "fixed" & term != "Intercept") %>%
    mutate(hyp.test = purrr::map(hyp, ~ hypothesis(c_model_nat, .x, seed = 2024))) %>%
    mutate(hypothesis = purrr::map(hyp.test, "hypothesis"))

  c_model_nat.hyp <- c_model_nat.tidy01 %>%
    pull(hypothesis) %>%
    map_dfr(~.x)

  c_model_nat.tidy02 <- c_model_nat.tidy %>%
    left_join(c_model_nat.hyp, by = "Hypothesis") %>%
    mutate(model = "c_model_nat")

  s_model_nat.tidy <- s_model_nat %>%
    tidy() %>%
    mutate(
      dir = ifelse(estimate > 0, " > 0", " < 0"),
      term = str_replace_all(term, "[()]", ""),
      hyp = paste0(response, "_", term, dir),
      Hypothesis = paste0("(", response, "_", term, ")", dir)
    )

  s_model_nat.tidy01 <- s_model_nat.tidy %>%
    filter(effect == "fixed" & term != "Intercept") %>%
    mutate(hyp.test = purrr::map(hyp, ~ hypothesis(s_model_nat, .x, seed = 2024))) %>%
    mutate(hypothesis = purrr::map(hyp.test, "hypothesis"))

  s_model_nat.hyp <- s_model_nat.tidy01 %>%
    pull(hypothesis) %>%
    map_dfr(~.x)

  s_model_nat.tidy02 <- s_model_nat.tidy %>%
    left_join(s_model_nat.hyp, by = "Hypothesis") %>%
    mutate(model = "s_model_nat")

  r_model_nat.tidy <- r_model_nat %>%
    tidy() %>%
    mutate(
      dir = ifelse(estimate > 0, " > 0", " < 0"),
      term = str_replace_all(term, "[()]", ""),
      hyp = paste0(response, "_", term, dir),
      Hypothesis = paste0("(", response, "_", term, ")", dir)
    )

  r_model_nat.tidy01 <- r_model_nat.tidy %>%
    filter(effect == "fixed" & term != "Intercept") %>%
    mutate(hyp.test = purrr::map(hyp, ~ hypothesis(r_model_nat, .x, seed = 2024))) %>%
    mutate(hypothesis = purrr::map(hyp.test, "hypothesis"))

  r_model_nat.hyp <- r_model_nat.tidy01 %>%
    pull(hypothesis) %>%
    map_dfr(~.x)

  r_model_nat.tidy02 <- r_model_nat.tidy %>%
    left_join(r_model_nat.hyp, by = "Hypothesis") %>%
    mutate(model = "r_model_nat")

  model_nat.tidy02 <- rbind(c_model_nat.tidy02,
                            s_model_nat.tidy02,
                            r_model_nat.tidy02)

  # Save data
  mnt02_file_name <- paste0("./results/phylo.sem.model_nat.tidy02_", region.i, ".csv")
  write_csv(model_nat.tidy02, mnt02_file_name)

  # ==========================================================================
  c_model_ias.tidy <- c_model_ias %>%
    tidy() %>%
    mutate(
      dir = ifelse(estimate > 0, " > 0", " < 0"),
      term = str_replace_all(term, "[()]", ""),
      hyp = paste0(response, "_", term, dir),
      Hypothesis = paste0("(", response, "_", term, ")", dir)
    )

  c_model_ias.tidy01 <- c_model_ias.tidy %>%
    filter(effect == "fixed" & term != "Intercept") %>%
    mutate(hyp.test = purrr::map(hyp, ~ hypothesis(c_model_ias, .x, seed = 2024))) %>%
    mutate(hypothesis = purrr::map(hyp.test, "hypothesis"))

  c_model_ias.hyp <- c_model_ias.tidy01 %>%
    pull(hypothesis) %>%
    map_dfr(~.x)

  c_model_ias.tidy02 <- c_model_ias.tidy %>%
    left_join(c_model_ias.hyp, by = "Hypothesis") %>%
    mutate(model = "c_model_ias")

  s_model_ias.tidy <- s_model_ias %>%
    tidy() %>%
    mutate(
      dir = ifelse(estimate > 0, " > 0", " < 0"),
      term = str_replace_all(term, "[()]", ""),
      hyp = paste0(response, "_", term, dir),
      Hypothesis = paste0("(", response, "_", term, ")", dir)
    )

  s_model_ias.tidy01 <- s_model_ias.tidy %>%
    filter(effect == "fixed" & term != "Intercept") %>%
    mutate(hyp.test = purrr::map(hyp, ~ hypothesis(s_model_ias, .x, seed = 2024))) %>%
    mutate(hypothesis = purrr::map(hyp.test, "hypothesis"))

  s_model_ias.hyp <- s_model_ias.tidy01 %>%
    pull(hypothesis) %>%
    map_dfr(~.x)

  s_model_ias.tidy02 <- s_model_ias.tidy %>%
    left_join(s_model_ias.hyp, by = "Hypothesis") %>%
    mutate(model = "s_model_ias")

  r_model_ias.tidy <- r_model_ias %>%
    tidy() %>%
    mutate(
      dir = ifelse(estimate > 0, " > 0", " < 0"),
      term = str_replace_all(term, "[()]", ""),
      hyp = paste0(response, "_", term, dir),
      Hypothesis = paste0("(", response, "_", term, ")", dir)
    )

  r_model_ias.tidy01 <- r_model_ias.tidy %>%
    filter(effect == "fixed" & term != "Intercept") %>%
    mutate(hyp.test = purrr::map(hyp, ~ hypothesis(r_model_ias, .x, seed = 2024))) %>%
    mutate(hypothesis = purrr::map(hyp.test, "hypothesis"))

  r_model_ias.hyp <- r_model_ias.tidy01 %>%
    pull(hypothesis) %>%
    map_dfr(~.x)

  r_model_ias.tidy02 <- r_model_ias.tidy %>%
    left_join(r_model_ias.hyp, by = "Hypothesis") %>%
    mutate(model = "r_model_ias")

  model_ias.tidy02 <- rbind(c_model_ias.tidy02,
                            s_model_ias.tidy02,
                            r_model_ias.tidy02)

  # Save data
  mit02_file_name <- paste0("./results/phylo.sem.model_ias.tidy02_", region.i, ".csv")
  write_csv(model_ias.tidy02, mit02_file_name)

  print(region.i)

}
