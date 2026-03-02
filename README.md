### Data and R code for "Effects of plant adaptive strategies, economic use, and minimum residence time on plant invasions in China are both stage and region dependent"

This README file lists and describes the files used for the analyses involved in the manuscript titled "Effects of plant adaptive strategies, economic use, and minimum residence time on plant invasions in China are both stage and region dependent". Please see the manuscript itself and the Supporting Information for additional details regarding this analysis.

### Files included

### Dataset of Alien Plants in China

- **checklist_accepted_matches_CSR.WithCVH.mrt.v6simple.20251023.csv**: this is the main dataset used in all analyses at regional scale. Each row represents a species record with CSR strategy scores, invasion status, economic uses, and geographic information.
- **checklist_accepted_matches_CSR_v4.20251023.csv**: dataset without MRT (Mean Residence Time) data.

### Regional Environmental Data

- **region_mean_Matti.hdi.csv**: Human Development Index (HDI) by geographic region.
- **region_mean_chelsa.bio1_12.csv**: climate variables (MAT and MAP) by geographic region.

### Phylogenetic tree Data

- **checklist_phylo.tree_graft_status_all.csr.sp.csv**: phylogenetic order and graft status for species in phylogenetic analysis;
- **checklist_phylo.tree_all.csr.sp.tre**: the phylogenetic tree in Newick format for all CSR species.

### R code used to conduct analyses

- **001_csr_cvh_phylo_brm_sem_national.R**: Bayesian phylogenetic SEM at national scale. Prepares data, loads phylogenetic trees, fits brms models;
- **001_csr_cvh_phylo_brm_sem_national_noMRT.R**: same as above without MRT variable (sensitivity analysis);
- **001_csr_cvh_phylo_brm_sem_regional.R**: Bayesian phylogenetic SEM at regional scale (7 geographic regions);
- **002_csr_cvh_integration_brm_sem_diagram_national.R**: generates SEM pathway diagrams for national-scale results;
- **002_csr_cvh_integration_brm_sem_diagram_regional.R**: generates SEM pathway diagrams for regional-scale results;
- **002_csr_cvh_brm_sem_diagram_national_noMRT.R**: SEM diagrams without MRT variable;
- **003_estimate_plotting_regional_analysis_noGDP_dotplot_fitting.R**: dot plots with linear regression fitting curves;
- **004_estimate_plotting_regional_analysis_radar_v4_barplot.R**: radar plots and bar plots for regional variations.

### Metadata for the main dataset

- The columns in the main dataset (**checklist_accepted_matches_CSR.WithCVH.mrt.v6simple.20251023.csv**) are as follows:
  - **taxon_name**: name of taxon, standardized using rWCVP;
  - **taxon_authors**: authority of the taxon name;
  - **accepted_plant_name_id**: accepted plant name ID from rWCVP;
  - **keep**: whether matched to an accepted name;
  - **family**: plant family;
  - **genus**: plant genus;
  - **species**: species epithet;
  - **Standardized_names**: standardized name using TPL;
  - **TPL_names**: names from The Plant List;
  - **TPL_Author**: authority of the standardized names from TPL;
  - **TPL_names.withAuthor**: full name with authority from TPL;
  - **invasion_status**: invasion status in China (1 = non-native, 2 = naturalized, 3 = invasive);
  - **c_score**: C-score (Competitor) - CSR strategy component, 0-100 scale;
  - **s_score**: S-score (Stress-tolerator) - CSR strategy component, 0-100 scale;
  - **r_score**: R-score (Ruderal) - CSR strategy component, 0-100 scale;
  - **wcup_eco_use_AF/ EU/ FU/ GS/ HF/ IF/ MA/ ME/ PO/ SU**: economic use categories;
  - **wcup_eco_use_num**: total number of economic-use categories;
  - **wcup_eco_use_num_2out**: economic uses excluding GS and PO;
  - **wcup_eco_use_status**: status of economic use;
  - **Checked_year**: year first recorded in China (herbarium record);
  - **Provience_pinyin**: province name (Pinyin);
  - **Geo.region**: geographic region (7 regions);
  - **Geo.region4**: 4-region classification;
  - **provinces**: all provinces where the species is recorded;
  - **n.provinces**: number of provinces;
  - **BG_codes**: herbarium specimen codes.
