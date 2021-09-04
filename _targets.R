library(targets)
library(tarchetypes)
library(tibble)
library(future)
library(future.callr)

source("functions.R")

tar_option_set(
  packages = c(
    "tidyverse", "glue", "rvest", "GenomicDataCommons", "TCGAutils", "biomaRt",
    "SummarizedExperiment", "DESeq2", "GSVA", "future", "future.callr", "gt",
    "gtsummary", "survival", "survminer", "webshot", "ragg"
  )
)

plan(callr)

values = tibble(
  project = c("blca", "skcm", "luad", "lusc", "kirc", "coad", "stad", "paad",
              "lihc", "hnsc")
)

mapped <- tar_map(
  values = values,
  
  # File structure and data acquisition ----------------------------------
  tar_target(dir, make_dir(paste0(data_dir, project, "/"))),
  tar_target(fig_dir, make_dir(paste0(dir, "figures/"))),
  tar_target(clin_dirty, get_clin(dir, project), format = "file"),
  
  # Wrangling and tidying ------------------------------------------------
  tar_target(clin, tidy_clin(clin_dirty)),
  tar_target(man, make_man(rm_cases, project)),
  tar_target(dds, man_to_dds(clin, man)),
  
  # Analysis -------------------------------------------------------------
  tar_target(norm, normalize(dds, gene_ids)),
  tar_target(gsva_scores, run_gsva(norm, gene_signatures)),
  tar_target(dds_w_scores, merge_gsva(norm, gsva_scores)),
  tar_target(dds_w_bin_scores, bin_gsva(dds_w_scores)),
  
  # Plots ---------------------------------------------------------------- 
  tar_target(clin_table, make_clin_table(dds_w_bin_scores, project)),
  tar_target(
    dens_b, 
    dense_ind(dds_w_bin_scores, 
              b_cell, 
              file_name = "density_b.png", 
              project = project,
              color = sex), 
    format = "file"
  ),
  tar_target(
    surv_b_sex, 
    surv_ind(dds_w_bin_scores, 
             strata = "b_bin", 
             file_name =  "surv_b_sex.png", 
             facet = "sex", 
             project = project), 
    format = "file"
  ),
  tar_target(
    surv_b, 
    surv_ind(dds_w_bin_scores, 
             strata = "b_bin", 
             file_name = "surv_b.png", 
             project = project), 
    format = "file"
  )
)

list(
  # Init dirs ------------------------------------------------------------
  
  tar_target(data_dir, make_dir("./01_data/")),
  tar_target(common_dir, make_dir(paste0(data_dir, "00_common/"))),
  tar_target(gene_signatures, tidy_signatures()),
  tar_target(rm_cases, make_rm(common_dir), format = "file"),
  mapped, 
  tar_combine(combined_dds, mapped[[6]], command = list(!!!.x)),
  tar_target(gene_ids, get_hgnc(combined_dds))
)
