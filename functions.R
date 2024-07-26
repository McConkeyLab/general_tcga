# Check for dependencies
check_if_gdc_dttc_exists <- function() {
  exists <- suppressWarnings(
    system2(
      "command", c("-v", "gdc-client"),
      stdout = FALSE,
      stderr = FALSE
    )
  )
  # Returns 0 if it exists

  if (exists != 0) {
    cli::cli_abort(
      c(
        "Could not find gdc-client",
        i = "gdc-client doesn't appear to be installed on this machine",
        i = "It might not be in your PATH",
        i = "You can download it here:",
        i = "https://bio-formats.readthedocs.io/en/stable/users/comlinetools/index.html" #nolint
      )
    )
  }
}


# File Preparation and Download -------------------------------------------

make_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
  path
}

make_rm <- function(tcga_common_dir) {

  base_path <- "https://gdac.broadinstitute.org/runs/stddata__latest/samples_report"

  download.file(
    glue("{base_path}/redactions_2016_01_28__00_00_08.tsv"),
    file.path(tcga_common_dir, "redacted.tsv"),
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

  redacted <- read_tsv(file.path(tcga_common_dir, "redacted.tsv"), show_col_types = FALSE)
  blacklisted <- read_tsv(file.path(tcga_common_dir, "blacklisted.tsv"), show_col_types = FALSE)
  ffpe <- read_tsv(file.path(tcga_common_dir, "ffpe.tsv"), show_col_types = FALSE)
  replicated <- read_tsv(file.path(tcga_common_dir, "replicated.txt"), show_col_types = FALSE)

  all <- tibble(id = c(redacted$Barcode,
                       blacklisted$`TCGA ID`,
                       ffpe$Barcode,
                       replicated$`Removed Samples`)) |>
    separate_rows(sep = ",")

  write_tsv(all, file.path(tcga_common_dir, "rm.tsv"))
  paste0(tcga_common_dir, "rm.tsv")
}

get_clin <- function(tcga_data_dir, tcga_proj) {

  file_to_download <-
    glue("https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/{toupper(tcga_proj)}/20160128/gdac.broadinstitute.org_{toupper(tcga_proj)}.Merge_Clinical.Level_1.2016012800.0.0.tar.gz")

  dest_file_loc <- paste0(tcga_data_dir, tcga_proj, "_clin.tar.gz")

  download.file(file_to_download, dest_file_loc)
  untar(dest_file_loc, exdir = paste0("./01_data/", tcga_proj, "/"))

  clin_dir <- str_extract(file_to_download, "[^/]*(?=.tar.gz)")

  clin_file <- list.files(paste0(tcga_data_dir, clin_dir), pattern = "clin.merged.txt")

  paste0(tcga_data_dir, clin_dir, "/", clin_file)
}

# Wrangle and Tidy -------------------------------------------------------------

get_latest <- function(data) {
  rowwise(data) |>
    mutate(across(everything(), as.numeric)) |>
    mutate(collapsed = if_else(!is.infinite(max(c_across(), na.rm = TRUE)),
                               max(c_across(), na.rm = TRUE),
                               NA_real_
    )) |>
    dplyr::select(collapsed)
}

tidy_clin <- function(clin_path, project) {

  clinical <- read_tsv(clin_path, col_names = FALSE, show_col_types = FALSE) |>
    t()

  colnames(clinical) <- clinical[1, ]

  clinical <- clinical |>
    as_tibble() |>
    dplyr::slice(-1)

  # Create Survival Time Column
  death <- clinical |>
    dplyr::select(contains("days_to_death")) |>
    get_latest() |>
    setNames("death_days")

  follow_up <- clinical |>
    dplyr::select(contains("days_to_last_followup")) |>
    get_latest() |>
    set_names("followUp_days")

  all_clin <-
    tibble(death, follow_up) |>
    cbind(clinical) |>
    mutate(
      follow_up_time = if_else(is.na(death_days), followUp_days, death_days),
      death = if_else(patient.vital_status == "dead" | !is.na(death_days), 1, 0)
    ) |>
    dplyr::select(-contains("days_to_death"),
                  -contains("days_to_last_followup"),
                  -patient.vital_status,
                  -death_days,
                  -followUp_days) |>
    dplyr::filter((follow_up_time >= 0) & !is.na(follow_up_time)) |>
    dplyr::filter((follow_up_time > 0) | (death == 0))
    # For reasoning as to the above filtering rule, see
    # https://www.graphpad.com/support/faq/events-deaths-at-time-zero-in-survival-analysis/

  clin_filtered <- all_clin |>
    recipe() |>
    step_zv(everything()) |>
    prep() |>
    bake(new_data = NULL)

  clin_filtered <- clin_filtered[, which(colSums(is.na(clin_filtered)) < nrow(clin_filtered) * 0.9)]

  clin_filtered$project <- project

  clin_filtered
}

make_man <- function(rm_cases_file, tcga_project) {

  tcga_project_gdc <- paste0("TCGA-", toupper(tcga_project))

  rm_cases <- read_tsv(rm_cases_file, show_col_types = FALSE)

  manifest <- files() |>
    GenomicDataCommons::filter(cases.project.project_id == tcga_project_gdc) |>
    GenomicDataCommons::filter(type == "gene_expression") |>
    GenomicDataCommons::filter(analysis.workflow_type == "STAR - Counts") |>
    GenomicDataCommons::filter(access == "open") |>
    manifest() |>
    # Creates a 'path' feature that DESeq uses for its sampleTable argument
    mutate(path = paste(id, file_name, sep = "/"))

  barcodes <- manifest$id |>
    UUIDtoBarcode(from_type = "file_id") |>
    dplyr::rename(submitter_id = associated_entities.entity_submitter_id)

  uuids <- manifest$id |>
    UUIDtoUUID(to_type = "case_id")

  joined_ids <- full_join(barcodes, uuids, by = "file_id")

  manifest <- manifest |>
    full_join(joined_ids, by = c("id" = "file_id")) |>
    dplyr::rename(submitter_id = submitter_id.y) |>
    mutate(short_id = str_sub(submitter_id, 1, 12)) |>
    dplyr::filter(str_detect(submitter_id, "^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-0")) |>
    anti_join(rm_cases, by = c("submitter_id" = "id")) |>
    arrange(submitter_id) |>
    distinct(short_id, .keep_all = TRUE)
}

download_tcga_data <- function(manifest, dir) {
  manifest <- manifest |>
    dplyr::select(id, path, submitter_id, short_id, cases.case_id)
  write_tsv(manifest, file = paste0(tempdir(), "/temp.txt"))
  system2(
    "gdc-client",
    c("download",
      "-m", paste0(tempdir(), "/temp.txt"),
      "-d", dir
    )
  )
  manifest
}

man_to_dds <- function(clin, manifest) {

  # Some samples may be in the manifest but not in the clinical data.

  # This is because clinical data was filtered to remove patients with
  # nonsensical or missing followup time (among other filters).

  # Additionally, some samples may be in the clinical data but not the manifest.

  # This is because sometimes samples have clinical data but did not have RNA
  # seq performed on them (one of the filtering criteria for the manifest
  # generation)

  tumor_data <- manifest |>
    inner_join(clin, by = c("cases.case_id" = "patient.bcr_patient_uuid"))

  read_star_file <- function(star_file_entry) {
    file_path <- star_file_entry$path
    sample_name <- star_file_entry$id
    fs::path("01_data", "00_gdcdata", file_path) |>
      read_tsv(skip = 1, col_types = "cc_i_____") |>
      dplyr::slice(-(1:4)) |>
      dplyr::rename(!!sample_name := unstranded)
  }

  read_and_join_star_file <- function(prev_star_file, new_entry) {
    full_join(prev_star_file, read_star_file(new_entry), by = c("gene_id", "gene_name"))
  }

  counts <- purrr::reduce(split(tumor_data, 1:nrow(tumor_data))[-1], 
                          read_and_join_star_file, 
                          .init = read_star_file(tumor_data[1,]))
  rd <- DataFrame(counts[1:2])
  rd$gene_id_w_ver <- rd$gene_id

  # Strip Ensembl Version Number
  rd$ensembl_gene_id <- str_replace(rd$gene_id, "\\.[0-9]+$", "")

  cd <- DataFrame(sampleName = tumor_data$id, fileName = tumor_data$path)
  se <- SummarizedExperiment(rowData = rd, colData = cd, assays = list(counts = data.matrix(counts[-(1:2)])))
  dds <- DESeqDataSet(se, design = ~1)


  samples <- tumor_data[match(colnames(dds), tumor_data$id), ]
  colData(dds) <- cbind(colData(dds), samples)
  colnames(dds) <- colData(dds)$submitter_id
  rownames(dds) <- rowData(dds)$ensembl_gene_id

  dds
}

get_hgnc <- function(tcga_dds = list()) {
  genes <- tcga_dds |>
    map(\(x){
      x |>
        rownames() |>
        as_tibble()
    }) |>
    purrr::reduce(bind_rows) |>
    distinct()

  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ids <- getBM(
    attributes = c("hgnc_symbol", "entrezgene_id", "ensembl_gene_id"),
    filters = "ensembl_gene_id",
    values = genes$value,
    mart = mart,
    useCache = FALSE
  )
}

join_ids <- function(dds, ids) {
  ids <- ids[match(rownames(dds), ids$ensembl_gene_id), ]
  rowData(dds) <- cbind(rowData(dds), ids)
  rownames(dds) <- rowData(dds)$hgnc_symbol
  rownames(dds) <- make.names(rownames(dds), unique = TRUE)
  nas <- which(is.na(rowData(dds)$ensembl_gene_id))
  if (!identical(nas, integer(0))) {
    dds <- dds[-nas, ]
  }
  dds
}

normalize <- function(dds) {
  norm_counts <- dds |>
    estimateSizeFactors() |>
    vst() |>
    assay()
  assay(dds, 2) <- norm_counts
  assayNames(dds)[[2]] <- "vst"
  dds
}
