library(targets)
library(tarchetypes)

source("functions.R")

tar_option_set(
  packages = c(
    "tidyverse", "broom", "glue", "rvest", "GenomicDataCommons", "TCGAutils",
    "biomaRt", "SummarizedExperiment", "DESeq2", "tidymodels"
  ),
  iteration = "list"
)

list(
  # Init dirs ------------------------------------------------------------------
  tar_target(data_dir, make_dir("./01_data/")),
  tar_target(common_dir, make_dir(paste0(data_dir, "00_common/"))),
  tar_target(projects, c("blca", "skcm", "luad", "lusc", "kirc", "coad", "stad",
                         "paad", "lihc", "hnsc")),
  tar_target(dir, make_dir(paste0(data_dir, projects, "/")), pattern = map(projects), iteration = "vector"),

  # Create aggregate file of samples to be removed -----------------------------
  tar_file(rm_cases, make_rm(common_dir)),

  # Clinical data acquisition and tidying --------------------------------------
  tar_file(clin_dirty, get_clin(dir, projects), pattern = map(dir, projects), iteration = "vector"),
  tar_target(clin, tidy_clin(clin_dirty, projects), pattern = map(clin_dirty, projects)),

  # Wrangle and tidy count data ------------------------------------------------
  tar_target(man, make_man(rm_cases, projects), pattern = map(projects)),
  tar_target(man_w_paths, download_tcga_data(man), pattern = map(man)),

  # Join count and clinical data in dds ----------------------------------------
  tar_target(dds, man_to_dds(clin, man_w_paths), pattern = map(clin, man_w_paths)),

  # Get and add HGNC gene symbols to dds ---------------------------------------
  tar_target(gene_ids, get_hgnc(dds)),
  tar_target(dds_w_ids, join_ids(dds, gene_ids), pattern = map(dds)),

  # Normalize and transform counts with VST and add to assay slot --------------
  tar_target(norm, normalize(dds_w_ids), pattern = map(dds_w_ids)),

  # Output ---------------------------------------------------------------------
  tar_target(rds, write_rds(norm, file = paste0(dir, "dds.rds")), pattern = map(norm, dir))
)
