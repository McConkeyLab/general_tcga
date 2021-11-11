library(targets)
library(tarchetypes)
library(tibble)
library(future)
library(future.callr)

source("functions.R")

tar_option_set(
  packages = c(
    "tidyverse", "broom", "glue", "rvest", "GenomicDataCommons", "TCGAutils", "biomaRt",
    "SummarizedExperiment", "DESeq2", "GSVA", "future", "future.callr", "gt",
    "gtsummary", "survival", "survminer", "webshot", "ragg", "car", "tidymodels",
    "censored", "ggsci"
  )
)

plan(callr)

values = tibble(
  project = c("blca", "skcm", "luad", "lusc", "kirc", "coad", "stad", "paad",
              "lihc", "hnsc")
)

mapped <- tar_map(
  values = values,
  
  # File structure and data acquisition ----------------------------------------
  tar_target(dir, make_dir(paste0(data_dir, project, "/"))),
  tar_target(fig_dir, make_dir(paste0(dir, "figures/"))),
  tar_target(clin_dirty, get_clin(dir, project), format = "file"),
  
  # Wrangling and tidying ------------------------------------------------------
  tar_target(clin, tidy_clin(clin_dirty)),
  tar_target(man, make_man(rm_cases, project)),
  tar_target(man_w_paths, download_tcga_data(man)),
  tar_target(dds, man_to_dds(clin, man_w_paths)),
  
  # Normalization --------------------------------------------------------------
  tar_target(norm, normalize(dds, gene_ids)),
  tar_target(norm_file, write_norm(norm, project)),
  
  # GSVA -----------------------------------------------------------------------
  tar_target(gsva_scores, run_gsva(norm, gene_signatures)),
  tar_target(dds_w_scores, merge_gsva(norm, gsva_scores)),
  tar_target(dds_w_bin_scores, bin_gsva(dds_w_scores)),
  
  # Final colData tibble -------------------------------------------------------
  tar_target(coldata_tibble, get_coldata_tibble(dds_w_bin_scores, project)),
  tar_target(gsva_tibble, get_gsva_tibble(coldata_tibble)),
  
  # Survival analyses ----------------------------------------------------------
  tar_target(surv_tidy, tidy_for_survival(coldata_tibble, project)),
  tar_target(surv_tidy_select, tidy_surv_select(surv_tidy)),
  tar_target(univariate, run_all_uni_combos(surv_tidy, project)),
  tar_target(multivariable_names, get_multivariable_names(univariate, project)),
  tar_target(multivariable, run_all_multi_combos(surv_tidy, multivariable_names, project)),
  tar_target(interaction_test_multi, test_interaction_multi(surv_tidy, multivariable_names, project)),
  tar_target(interaction_test_uni, test_interaction_uni(surv_tidy, project)),
  
  # Tables ---------------------------------------------------------------------
  tar_target(clin_table, make_clin_table(dds_w_bin_scores, project), format = "file"),

  
  # Plots ----------------------------------------------------------------------
  tar_target(univariate_plots,
             make_univariate_plot(univariate, project)),
  
  tar_target(multivariable_plots,
             make_all_multivariable_plot_combos(multivariable, project)),

  tar_target(survival_plots,
             surv_ind_all(data = dds_w_bin_scores,
                          strata = c("b_bin", "cd8_bin", "imm_bin"),
                          project = project,
                          facet = c(NA, "sex")),
             format = "file")
)

list(
  # Init dirs ------------------------------------------------------------
  tar_target(data_dir, make_dir("./01_data/")),
  tar_target(common_dir, make_dir(paste0(data_dir, "00_common/"))),
  tar_target(gene_signatures, tidy_signatures()),
  tar_target(rm_cases, make_rm(common_dir), format = "file"),
  mapped, 
  tar_combine(combined_dds, mapped[["dds"]], command = list(!!!.x)),
  tar_target(gene_ids, get_hgnc(combined_dds)),
  tar_combine(combined_unis, mapped[["univariate"]], command = rbind(!!!.x)),
  tar_combine(combined_multis, mapped[["multivariable"]], command = rbind(!!!.x)),
  tar_combine(combined_survs, mapped[["surv_tidy_select"]], command = rbind(!!!.x)),
  tar_combine(combined_ints, mapped[["interaction_test_multi"]], mapped[["interaction_test_uni"]], command = rbind(!!!.x)),
  tar_combine(gsva_coldata, mapped[["gsva_tibble"]], command = rbind(!!!.x)),
  tar_target(density_plot, make_density_plot(gsva_coldata), format = "file"),
  tar_target(int_output, tidy_and_write_int(combined_ints), format = "file"),
  tar_target(hr_plot_unis, make_hr_plot(combined_unis), format = "file"),
  tar_target(hr_plot_multi, make_hr_plot_multi(combined_multis), format = "file")
)
