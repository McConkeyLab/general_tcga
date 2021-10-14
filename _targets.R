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
  tar_target(man_w_paths, download_tcga_data(man)),
  tar_target(dds, man_to_dds(clin, man_w_paths)),
  
  # Analysis -------------------------------------------------------------
  tar_target(norm, normalize(dds, gene_ids)),
  tar_target(norm_file, write_norm(norm, project)),
  tar_target(gsva_scores, run_gsva(norm, gene_signatures)),
  tar_target(dds_w_scores, merge_gsva(norm, gsva_scores)),
  tar_target(dds_w_bin_scores, bin_gsva(dds_w_scores)),
  tar_target(b_hrs, get_hr(dds_w_bin_scores, "b_bin + age + pathologic_tnm_t", project)),
  tar_target(t_hrs, get_hr(dds_w_bin_scores, "cd8_bin + age + pathologic_tnm_t", project)),
  tar_target(imm_hrs, get_hr(dds_w_bin_scores, "imm_bin + age + pathologic_tnm_t", project)),
  tar_target(surv_tidy, tidy_for_survival(dds_w_bin_scores, project)),
  tar_target(univariate_tidy, univariate(surv_tidy)),
  tar_target(univariate_glance, univariate(surv_tidy, TRUE)),
  tar_target(univariate_both, full_join(univariate_tidy, univariate_glance, by = c("name", "stratum"))),
  tar_target(multivariable_names, get_multivariable_names(univariate_glance)),
  tar_target(multivariable, run_all_multi_combos(surv_tidy, multivariable_names, project)),
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
  ),
  tar_target(
    dens_cd8, 
    dense_ind(dds_w_bin_scores, 
              cd8, 
              file_name = "density_cd8.png", 
              project = project,
              color = sex), 
    format = "file"
  ),
  tar_target(
    surv_cd8_sex, 
    surv_ind(dds_w_bin_scores, 
             strata = "cd8_bin", 
             file_name =  "surv_cd8_sex.png", 
             facet = "sex", 
             project = project), 
    format = "file"
  ),
  tar_target(
    surv_cd8, 
    surv_ind(dds_w_bin_scores, 
             strata = "cd8_bin", 
             file_name = "surv_cd8.png", 
             project = project), 
    format = "file"
  ),
  tar_target(
    dens_imm, 
    dense_ind(dds_w_bin_scores, 
              exp_immune, 
              file_name = "density_pan-imm.png", 
              project = project,
              color = sex), 
    format = "file"
  ),
  tar_target(
    surv_imm_sex, 
    surv_ind(dds_w_bin_scores, 
             strata = "imm_bin", 
             file_name =  "surv_pan-imm_sex.png", 
             facet = "sex", 
             project = project), 
    format = "file"
  ),
  tar_target(
    surv_imm, 
    surv_ind(dds_w_bin_scores, 
             strata = "imm_bin", 
             file_name = "surv_pan-imm.png", 
             project = project), 
    format = "file"
  ),
  tar_target(
    surv_sex,
    surv_ind(dds_w_bin_scores,
             strata = "sex",
             file_name = "surv_sex.png",
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
  tar_combine(combined_dds, mapped[[7]], command = list(!!!.x)),
  tar_target(gene_ids, get_hgnc(combined_dds)),
  tar_combine(combined_hrs, mapped[[13]], mapped[[14]], mapped[[15]], command = rbind(!!!.x)),
  tar_target(hr_plot, make_hr_plot(combined_hrs), format = "file"),
  tar_combine(combined_multis, mapped[[21]], command = rbind(!!!.x)),
  tar_target(hr_plot_multi, make_hr_plot_multi(combined_multis), format = "file")
)
