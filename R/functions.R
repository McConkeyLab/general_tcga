make_directory <- function(path) {
        if (!dir.exists(path)) {
                dir.create(path)
        }
        path
}

# TCGA common ------------------------------------------------------------------

make_tcga_rm_file <- function(tcga_common_dir) {
        
        base_path <- "http://gdac.broadinstitute.org/runs/stddata__latest/samples_report"
        download.file(
                glue("{base_path}/redactions_2016_01_28__00_00_08.tsv"),
                file.path(tcga_common_dir, "redacted.tsv")
        )
        download.file(
                glue("{base_path}/blacklist.tsv"),
                file.path(tcga_common_dir, "blacklisted.tsv")
        )
        download.file(
                glue("https://gdac.broadinstitute.org/runs/tmp/sample_report__2018_01_17/TCGA.filtered_samples.txt"),
                file.path(tcga_common_dir, "replicated.txt")
        )
        read_html(glue("{base_path}/FFPE_Cases.html")) |>
                html_node("table") |>
                html_table() |>
                write_tsv(file.path(tcga_common_dir, "ffpe.tsv"))
        
        redacted <- read_tsv(file.path(tcga_common_dir, "redacted.tsv"))
        blacklisted <- read_tsv(file.path(tcga_common_dir, "blacklisted.tsv"))
        ffpe <- read_tsv(file.path(tcga_common_dir, "ffpe.tsv"))
        replicated <- read_tsv(file.path(tcga_common_dir, "replicated.txt"))
        
        all <- tibble(id = c(redacted$Barcode, 
                             blacklisted$`TCGA ID`, 
                             ffpe$Barcode, 
                             replicated$`Removed Samples`)) |> 
                separate_rows(sep = ",")
        
        write_tsv(all, file.path(tcga_common_dir, "rm.tsv"))
        paste0(tcga_common_dir, "rm.tsv")
}

download_tcga_clin_file <- function(tcga_data_dir, tcga_proj) {
        
        file_to_download <- 
                glue("http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/{toupper(tcga_proj)}/20160128/gdac.broadinstitute.org_{toupper(tcga_proj)}.Merge_Clinical.Level_1.2016012800.0.0.tar.gz")
        
        dest_file_loc <- paste0(tcga_data_dir, tcga_proj, "_clin.tar.gz")
        
        download.file(file_to_download, dest_file_loc)
        untar(dest_file_loc, exdir = paste0("./01_data/", tcga_proj, "/"))
        
        clin_dir <- str_extract(file_to_download, "[^/]*(?=.tar.gz)")
        
        clin_file <- list.files(paste0(tcga_data_dir, clin_dir), pattern = "clin.merged.txt")
        
        paste0(tcga_data_dir, clin_dir, "/", clin_file)
}

tidy_tcga_clinical <- function(tcga_clin_path) {
        clinical <-
                read_tsv(tcga_clin_path,
                         col_names = F
                ) |>
                t()
        
        colnames(clinical) <- clinical[1, ]
        
        clinical <- clinical |>
                as_tibble()
        clinical <- clinical[-1, ]
        
        # Remove Columns with Only One Unique Entry
        
        is_just_one <- function(x) {
                length(unique(x)) == 1
        }
        
        res <- apply(clinical, 2, is_just_one)
        
        clin_filtered <- clinical[, !res] |>
                as_tibble()
        
        
        # Create Survival Time Column
        
        get_latest <- function(data) {
                rowwise(data) |>
                        mutate(across(everything(), as.numeric)) |>
                        mutate(collapsed = if_else(!is.infinite(max(c_across(), na.rm = T)),
                                                   max(c_across(), na.rm = T),
                                                   NA_real_
                        )) |>
                        dplyr::select(collapsed)
        }
        
        death <-
                clin_filtered |>
                dplyr::select(contains("days_to_death")) |>
                get_latest()
        
        fl <-
                clin_filtered |>
                dplyr::select(contains("days_to_last_followup")) |>
                get_latest()
        
        all_clin <-
                tibble(death, fl, .name_repair = "unique") |>
                `colnames<-`(c("death_days", "followUp_days")) |>
                mutate(
                        new_death = if_else(is.na(death_days), followUp_days, death_days),
                        death_event = if_else(clinical$patient.vital_status == "dead" | !is.na(death_days), 1, 0)
                ) |>
                cbind(clin_filtered)
        
        path <- str_remove(tcga_clin_path, "[^/]*$")
        write_tsv(all_clin, paste0(path, "clin.tsv"))
        paste0(path, "clin.tsv")
}

create_tcga_manifest <- function(rm_cases_file, tcga_project) {
        
        tcga_project_gdc <- paste0("TCGA-", toupper(tcga_project))

        rm_cases <- read_tsv(rm_cases_file)
        
        manifest <- files() |>
                GenomicDataCommons::filter(cases.project.project_id == tcga_project_gdc) |>
                GenomicDataCommons::filter(type == "gene_expression") |>
                GenomicDataCommons::filter(analysis.workflow_type == "HTSeq - Counts") |>
                manifest() |>
                # Creates a 'path' feature that DESeq uses for its sampleTable argument
                mutate(path = paste(id, filename, sep = "/"))
        
        uuids <- manifest$id |>
                UUIDtoBarcode(from_type = "file_id") |>
                dplyr::rename(submitter_id = associated_entities.entity_submitter_id)
        
        # Match Counts to Clinical
        manifest <- manifest |>
                full_join(uuids, by = c("id" = "file_id")) |>
                mutate(short_id = str_sub(submitter_id, 1, 12))
        
        # Keep just tumor cases
        #manifest <- manifest |>
        #        dplyr::filter(str_detect(submitter_id, "^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-0"))
        
        # Remove FFPE, blacklisted, and redacted samples
        manifest <- manifest[!((manifest$submitter_id %in% rm_cases$id)), ]
        file_path <- paste0("./01_data/", tcga_project, "/manifest.tsv")
        write_tsv(manifest, file_path)
        file_path
}

tcga_manifest_to_dds <- function(clin_path, manifest) {
        clin <- read_tsv(clin_path)
        manifest <- read_tsv(manifest) |>
                mutate(short_id = tolower(short_id)) |>
                dplyr::select(id, path, submitter_id, short_id)
        gdc_set_cache("./01_data/gdcdata", create_without_asking = T)
        lapply(manifest$id, gdcdata)
        # Attach Tumor Calls and Clinical Data
        tumor_data <- manifest |>
                left_join(clin, by = c("short_id" = "patient.bcr_patient_barcode")) |>
                mutate(
                        short_id = toupper(short_id),
                        case_id = make.names("short_id", unique = T)
                )
        # Generate Sample Table
        sample_table <- data.frame(
                sampleName = tumor_data$id,
                fileName = tumor_data$path
        )
        dds <- DESeqDataSetFromHTSeqCount(
                sampleTable = sample_table,
                directory = gdc_cache(),
                design = ~1
        )
        samples <- tumor_data[match(colnames(dds), tumor_data$id), ]
        colData(dds) <- cbind(colData(dds), samples)
        colnames(dds) <- colData(dds)$short_id
        
        # Strip Ensembl Version Number
        rownames(dds) <- str_replace(rownames(dds), ".[0-9]+$", "")
        write_rds(dds, paste0(str_remove(clin_path, "[^/]*/clin.tsv$"), "dds.Rds"))
        paste0(str_remove(clin_path, "[^/]*/clin.tsv$"), "dds.Rds")
}

get_hgnc <- function(tcga_dds = list()) {
        genes <- tcga_dds |> 
                map(\(x) {
                        x |> 
                                read_rds() |> 
                                rownames()
                }) |> 
                unlist() |> 
                as_tibble() |> 
                distinct()

        mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        ids <- getBM(
                attributes = c("hgnc_symbol", "entrezgene_id", "ensembl_gene_id"),
                filters = "ensembl_gene_id",
                values = genes$value,
                mart = mart,
                useCache = FALSE
        )

        write_tsv(ids, "./01_data/tcga-common/all_ids.tsv")
        "./01_data/tcga-common/all_ids.tsv"
}

norm_tcga_counts <- function(dds_path, hgnc_ids) {
        dds <- read_rds(dds_path)
        
        ids <- read_tsv(hgnc_ids)
        
        # Many (~5%, or 3000) transcripts mapped to now deprecated genes
        # and thus return NA
        
        ids <- ids[match(rownames(dds), ids$ensembl_gene_id), ]
        rowData(dds) <- cbind(rowData(dds), ids)
        
        rownames(dds) <- rowData(dds)$hgnc_symbol
        
        # Remove NA and Blank HUGO Symbols
        dds <- dds[!is.na(rownames(dds)), ]
        dds <- dds[(rownames(dds) != ""), ]
        
        # Normalize and Scale Counts
        norm_counts <- dds |>
                estimateSizeFactors() |>
                vst() |>
                assay()
        
        assay(dds, 2) <- norm_counts
        assayNames(dds)[[2]] <- "vst"
        
        write_rds(dds, paste0(str_remove(dds_path, "dds.Rds"), "norm-counts.Rds"))
        paste0(str_remove(dds_path, "dds.Rds"), "norm-counts.Rds")
}

tidy_signatures <- function(signatures_file) {
        
        # Pan B-Cell Signature: https://jitc.bmj.com/content/5/1/18.long
        
        # CD8+ T Effector Signature (cd8_rosen): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5480242/#SD1
        # Supplementary figure 6
        
        # CD8+ T (cd8_prat): https://cancerres.aacrjournals.org/content/canres/77/13/3540/F1.large.jpg?width=800&height=600&carousel=1
        
        # CD8+ T (cd8_sade): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6641984/bin/NIHMS1510803-supplement-11.xlsx
        
        # CD8+ T (cd8_fehr): https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(16)00587-0/fulltext
        
        # IFNg Signature: https://www.jci.org/articles/view/91190/table/2
        
        # Myeloid Inflammation Signature: https://www.nature.com/articles/s41591-018-0053-3/figures/2
        
        b_cell <- c(
                "BLK", "CD19", "FCRL2", "MS4A1", "KIAA0125", "TNFRSF17", "TCL1A",
                "SPIB", "PNOC"
        )
        cd8_rosen <- c(
                "CD8A", "GZMA", "GZMB", "IFNG", "CXCL9", "CXCL10", "PRF1", "TBX21"
        )
        cd8_prat <- c("PRF1", "CD8A", "CD8B", "GZMM", "FLT3LG")
        cd8_sade <- c(
                "IL7R", "GPR183", "LMNA", "NR4A3", "TCF7", "MGAT4A", "CD55", "AIM1", "PER1",
                "FOSL2", "EGFR1", "TSPYL2", "YPEL5", "CSRNP1", "REL", "SKIL", "PIK3R1",
                "FOXP1", "RGCC", "PFKFB3", "MYADM", "ZFP36L2", "USP36", "TC2N", "FAM177A1",
                "BTG2", "TSC22D2", "FAM65B", "STAT4", "RGP5", "NEU1", "IRFD1", "PDE4B",
                "NR4A1"
        )
        cd8_fehr <- c("CD8A", "EOMES", "PRF1", "IFNG", "CD274")
        ifng <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG")
        myeloid_inf <- c("CXCL1", "CXCL2", "CXCL3", "CXCL8", "IL6", "PTGS2")
        
        signatures <- list(
                b_cell = b_cell, cd8_rose = cd8_rosen,
                cd8_prat = cd8_prat, cd8_fehr = cd8_fehr, ifng = ifng,
                myeloid_inf = myeloid_inf
        )
        
        write_rds(signatures, paste0(signatures_file, "signatures.Rds"))
        paste0(signatures_file, "signatures.Rds")
}

run_gsva <- function(data_path, signatures_path) {
        data <- read_rds(data_path)
        signatures <- read_rds(signatures_path)
        
        # Foresight: KIAA0125 --> FAM30A
        signatures$b_cell[(signatures$b_cell == "KIAA0125")] <- "FAM30A"
        
        genes <- unlist(signatures) |>
                as_tibble()
        no_match <- as.character(genes[!(genes$value %in% rownames(data)), ])
        
        if (no_match == "character(0)") {
                print("All signature genes have match in dataset.")
        } else {
                warning(paste(no_match, "is in a signature but has no match in the dataset\n"))
        }
        
        scores <- gsva(assay(data, 2),
                       signatures,
                       mx.diff = T,
                       kcdf = "Gaussian"
        )
        
        write_rds(scores, paste0(str_remove(data_path, "norm-counts.Rds"), "scores.Rds"))
        paste0(str_remove(data_path, "norm-counts.Rds"), "scores.Rds")
}

merge_gsva <- function(data_path, scores) {
        data <- read_rds(data_path)
        
        data_coldata <- data |>
                colData() |>
                as_tibble(rownames = "sample")
        
        gsva <- read_rds(scores) |>
                t() |>
                as_tibble(rownames = "sample")
        
        joined <- inner_join(data_coldata, gsva, by = "sample")
        
        colData(data) <- DataFrame(joined)
        
        write_rds(data, paste0(str_remove(data_path, "norm-counts.Rds"), "dds-w-scores.Rds"))
        paste0(str_remove(data_path, "norm-counts.Rds"), "dds-w-scores.Rds")
}

theme_tufte <- function(font_size = 30) {
        theme(
                panel.grid = element_blank(),
                panel.background = element_rect(fill = "#FFFFF8", color = "#CCCCCC"),
                plot.background = element_rect(fill = "#FFFFF8"),
                strip.background = element_rect(fill = "#BBBBB0"),
                legend.background = element_rect(fill = "#FFFFF8"),
                legend.position = "top",
                legend.key = element_blank(),
                text = element_text(family = "GillSans", size = font_size)
        )
}

get_gill <- function() {
        showtext_auto()
        
        if (Sys.info()[["sysname"]] == "Windows") {
                font_add("GillSans", "GIL_____.TTF")
        } else if (Sys.info()[["sysname"]] == "Darwin") {
                font_add("GillSans", "GillSans.ttc")
        } else if (Sys.info()[["sysname"]] == "Linux") {
                font_add("GillSans", "Gill Sans.otf")
        }
}


test_plot <- function(tcga_dds = list()) {
        get_gill()
        sigs <- tcga_dds |> 
                map(\(x) {
                        x <- x |> 
                                read_rds()
                        x <- cbind(x$b_cell, x$cd8_rose) |>
                                as_tibble() |> 
                                setNames(c("b", "cd8t"))
                }) |> 
                bind_rows(.id = "proj")
        ggplot(sigs, aes(x = b, y = cd8t, color = proj)) + 
                scale_color_viridis_d(option = "plasma", end = 0.8) +
                geom_point(size = 3, alpha = 0.5, shape = 16) + 
                facet_wrap(~proj) +
                theme_tufte() + 
                theme(legend.position = "none")
        ggsave("./test-plot.png", width = 10, height = 10)
                
}
