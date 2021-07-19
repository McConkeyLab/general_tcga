library(targets)
library(tarchetypes)
library(tibble)

source("R/functions.R")

tar_option_set(
        packages = c(
                "tidyverse", "glue", "rvest", "GenomicDataCommons", "TCGAutils",
                "biomaRt", "SummarizedExperiment", "DESeq2", "GSVA", "future",
                "future.callr", "showtext", "gt", "gtsummary", "survival", 
                "survminer"
        )
)

library(future)
library(future.callr)
plan(callr)

values = tibble(
        project = c("blca", "brca", "ov", "skcm", "luad", "lusc", "kirc", "coad",
                    "stad", "paad", "lihc", "hnsc")
)


mapped <- tar_map(
        values = values,
        tar_target(dir, make_directory(paste0(data_dir, project, "/"))),
        tar_target(clin_file, download_tcga_clin_file(dir, project), format = "file"),
        tar_target(tidy_tcga_clin, tidy_tcga_clinical(clin_file), format = "file"),
        tar_target(tcga_man, create_tcga_manifest(rm_cases, project), format = "file"),
        tar_target(tcga_dds, tcga_manifest_to_dds(tidy_tcga_clin, tcga_man), format = "file"),
        tar_target(tcga_norm, norm_tcga_counts(tcga_dds, hgnc_symbols), format = "file"),
        tar_target(tcga_gsva_scores, run_gsva(tcga_norm, gene_signatures), format = "file"),
        tar_target(tcga_gsva_merged, merge_gsva(tcga_norm, tcga_gsva_scores), format = "file"),
        tar_target(tcga_gsva_tumor, rm_normal_tissue(tcga_gsva_merged), format = "file"),
        tar_target(tcga_gsva_unique_tumor, select_first_duplicate(tcga_gsva_tumor), format = "file"),
        tar_target(tcga_tidy_clin, prep_clin_data(tcga_gsva_unique_tumor), format = "file"),
        tar_target(tcga_clin_table, make_clin_table(tcga_tidy_clin), format = "file"),
        tar_target(tcga_dens_b, dense_ind(tcga_tidy_clin, b_cell, "density_b.png", color = sex), format = "file"),
        tar_target(tcga_surv_b_sex, surv_ind(tcga_tidy_clin, "b_bin", file_name =  "surv_b_sex.png", facet = "sex"), format = "file"),
        tar_target(tcga_surv_b, surv_ind(tcga_tidy_clin, strata = "b_bin", file_name = "surv_b.png"), format = "file")
)

list(
        # Init dirs ------------------------------------------------------------
        
        tar_target(
                data_dir,
                make_directory("./01_data/")
        ),
        tar_target(
                tcga_common_dir,
                make_directory(paste0(data_dir, "tcga-common/"))
        ),
        tar_target(
                signatures_dir,
                make_directory(paste0(data_dir, "signatures/"))
        ),
        tar_target(
                gene_signatures,
                tidy_signatures(signatures_dir),
                format = "file"
        ),
        tar_target(
                rm_cases,
                make_tcga_rm_file(tcga_common_dir),
                format = "file"
        ),
        
        mapped, 
        
        tar_combine(combined_dds, mapped[[5]]),
        
        tar_target(hgnc_symbols, get_hgnc(combined_dds)),
        
        tar_combine(combined_gsva, mapped[[10]]),
        
        tar_combine(combined_clin, mapped[[11]]),
        tar_target(b_surv, survival_b_cell(combined_clin)),
        tar_target(b_dens, density_b_cell(combined_clin)),
        tar_target(pan_dens, density_pan_immune_cell(combined_clin)),
        
        tar_target(tp, test_plot(combined_gsva))
)
