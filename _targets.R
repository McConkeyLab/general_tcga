library(targets)
library(tarchetypes)
library(tibble)

source("functions.R")

tar_option_set(
        packages = c(
                "tidyverse", "glue", "rvest", "GenomicDataCommons", "TCGAutils",
                "biomaRt", "SummarizedExperiment", "DESeq2", "GSVA", "future",
                "future.callr", "showtext", "gt", "gtsummary", "survival", 
                "survminer", "webshot", "ragg"
        )
)

library(future)
library(future.callr)
plan(callr)

values = tibble(
        project = c("blca", "skcm", "luad", "lusc", "kirc", "coad",
                    "stad", "paad", "lihc", "hnsc")
)


mapped <- tar_map(
        values = values,
        tar_target(dir, make_dir(paste0(data_dir, project, "/"))),
        tar_target(fig_dir, make_dir(paste0(dir, "figures/"))),
        tar_target(clin_file, get_clin_file(dir, project), format = "file"),
        tar_target(tidy_clin, tidy_clin(clin_file)),
        tar_target(man, make_man(rm_cases, project)),
        tar_target(dds, man_to_dds(tidy_clin, man)),
        tar_target(norm, norm(dds, gene_ids)),
        tar_target(gsva_scores, run_gsva(norm, gene_signatures)),
        tar_target(dds_w_scores, merge_gsva(norm, gsva_scores)),
        tar_target(dds_w_bin_scores, bin_gsva(dds_w_scores)),
        tar_target(clin_table, make_clin_table(dds_w_bin_scores, project)),
        # tar_target(dens_b, dense_ind(tcga_tidy_clin, b_cell, "density_b.png", color = sex), format = "file"),
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
        tar_target(signatures_dir, make_dir(paste0(data_dir, "signatures/"))),
        tar_target(gene_signatures, tidy_signatures(signatures_dir)),
        tar_target(rm_cases, make_rm(common_dir), format = "file"),
        mapped, 
        tar_combine(combined_dds, mapped[[6]], command = list(!!!.x)),
        tar_target(gene_ids, get_hgnc(combined_dds))
        
        
        # tar_combine(combined_gsva, mapped[[10]]),
        # 
        # tar_combine(combined_clin, mapped[[11]]),
        # tar_target(b_surv, survival_b_cell(combined_clin))
        # tar_target(b_dens, density_b_cell(combined_clin)),
        # tar_target(pan_dens, density_pan_immune_cell(combined_clin)),
        # 
        # tar_target(tp, test_plot(combined_gsva))
)


# TODO: 
# Density