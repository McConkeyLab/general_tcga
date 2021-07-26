# File Preparation and Download ------------------------------------------------

make_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
  path
}

make_rm <- function(tcga_common_dir) {
  
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

get_clin<- function(tcga_data_dir, tcga_proj) {
  
  file_to_download <- 
    glue("http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/{toupper(tcga_proj)}/20160128/gdac.broadinstitute.org_{toupper(tcga_proj)}.Merge_Clinical.Level_1.2016012800.0.0.tar.gz")
  
  dest_file_loc <- paste0(tcga_data_dir, tcga_proj, "_clin.tar.gz")
  
  download.file(file_to_download, dest_file_loc)
  untar(dest_file_loc, exdir = paste0("./01_data/", tcga_proj, "/"))
  
  clin_dir <- str_extract(file_to_download, "[^/]*(?=.tar.gz)")
  
  clin_file <- list.files(paste0(tcga_data_dir, clin_dir), pattern = "clin.merged.txt")
  
  paste0(tcga_data_dir, clin_dir, "/", clin_file)
}

# Wrangle and Tidy -------------------------------------------------------------

tidy_clin <- function(tcga_clin_path) {
  
  clinical <-
    read_tsv(tcga_clin_path,
             col_names = F,
             show_col_types = FALSE) |>
    t() 
  
  colnames(clinical) <- clinical[1, ]
  
  clinical <- clinical |>
    as_tibble() |> 
    dplyr::slice(-1) |> 
    dplyr::select(contains("days_to_death"),
                  contains("days_to_last_followup"),
                  patient.vital_status,
                  patient.bcr_patient_uuid,
                  sex = patient.gender,
                  age = patient.age_at_initial_pathologic_diagnosis,
                  race = patient.race_list.race,
                  days_to_collection = matches("^patient.samples.sample.days_to_collection$"),
                  path_stage = matches("^patient.stage_event.pathologic_stage$"),
                  tnm_t = matches("^patient.stage_event.tnm_categories.pathologic_categories.pathologic_t$"),
                  tnm_n = matches("^patient.stage_event.tnm_categories.pathologic_categories.pathologic_n$"),
                  tnm_m = matches("^patient.stage_event.tnm_categories.pathologic_categories.pathologic_m$"),
                  anatomic_subdivision = matches("^patient.anatomic_neoplasm_subdivision$"))
  
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
  
  death <- clinical |>
    dplyr::select(contains("days_to_death")) |>
    get_latest() |> 
    setNames("death_days")
  
  fu <- clinical |>
    dplyr::select(contains("days_to_last_followup")) |>
    get_latest() |> 
    set_names("followUp_days")
  
  all_clin <-
    tibble(death, fu) |>
    cbind(clinical) |> 
    mutate(
      new_death = if_else(is.na(death_days), followUp_days, death_days),
      death_event = if_else(patient.vital_status == "dead" | !is.na(death_days), 1, 0)
    ) |>
    dplyr::select(-contains("days_to_death"),
                  -contains("days_to_last_followup"),
                  -patient.vital_status) |> 
    dplyr::filter((new_death >= 0) & !is.na(new_death)) |>
    dplyr::filter((new_death > 0) | (death_event == 0)) |> 
    # For reasoning as to the above filtering rule, see
    # https://www.graphpad.com/support/faq/events-deaths-at-time-zero-in-survival-analysis/
    dplyr::filter(!is.na(sex)) |> 
    mutate(sex = factor(sex, levels = c("male", "female"), labels = c("M", "F")),
           age = as.numeric(age),
           days_to_collection = as.numeric(days_to_collection))
}

make_man <- function(rm_cases_file, tcga_project) {
  
  tcga_project_gdc <- paste0("TCGA-", toupper(tcga_project))
  
  rm_cases <- read_tsv(rm_cases_file, show_col_types = FALSE)
  
  manifest <- files() |>
    GenomicDataCommons::filter(cases.project.project_id == tcga_project_gdc) |>
    GenomicDataCommons::filter(type == "gene_expression") |>
    GenomicDataCommons::filter(analysis.workflow_type == "HTSeq - Counts") |>
    manifest() |>
    # Creates a 'path' feature that DESeq uses for its sampleTable argument
    mutate(path = paste(id, filename, sep = "/"))
  
  barcodes <- manifest$id |>
    UUIDtoBarcode(from_type = "file_id") |>
    dplyr::rename(submitter_id = associated_entities.entity_submitter_id)
  
  uuids <- manifest$id |> 
    UUIDtoUUID(to_type = "case_id")
  
  joined_ids <- full_join(barcodes, uuids, by = "file_id")
  
  manifest <- manifest |>
    full_join(joined_ids, by = c("id" = "file_id")) |>
    mutate(short_id = str_sub(submitter_id, 1, 12)) |> 
    dplyr::filter(str_detect(submitter_id, "^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-0")) |> 
    anti_join(rm_cases, by = c("submitter_id" = "id")) |> 
    arrange(submitter_id) |> 
    distinct(short_id, .keep_all = TRUE)
}

man_to_dds <- function(clin, manifest) {
  manifest <- manifest |> 
    dplyr::select(id, path, submitter_id, short_id, cases.case_id)
  gdc_set_cache("./01_data/gdcdata", create_without_asking = T)
  lapply(manifest$id, gdcdata)
  
  
  # Some samples may be in the manifest but not in the clinical data.
  
  # This is because clinical data was filtered to remove patients with
  # nonsensical or missing followup time (among other filters).
  
  # Additionally, some samples may be in the clinical data but not the manifest.
  
  # This is because sometimes samples have clinical data but did not have RNA
  # seq performed on them (one of the filtering criteria for the manifest
  # generation)
  
  tumor_data <- manifest |>
    inner_join(clin, by = c("cases.case_id" = "patient.bcr_patient_uuid"))
  
  sample_table <- data.frame(sampleName = tumor_data$id,
                             fileName = tumor_data$path)
  
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table,
                                    directory = gdc_cache(),
                                    design = ~1)
  
  samples <- tumor_data[match(colnames(dds), tumor_data$id), ]
  colData(dds) <- cbind(colData(dds), samples)
  colnames(dds) <- colData(dds)$submitter_id
  
  # Strip Ensembl Version Number
  rownames(dds) <- str_replace(rownames(dds), ".[0-9]+$", "")
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

norm <- function(dds, ids) {
  
  ids <- ids[match(rownames(dds), ids$ensembl_gene_id), ]
  rowData(dds) <- cbind(rowData(dds), ids)
  rownames(dds) <- rowData(dds)$hgnc_symbol
  rownames(dds) <- make.names(rownames(dds), unique = TRUE)
  nas <- which(is.na(rowData(dds)$ensembl_gene_id))
  dds <- dds[-nas, ]
  
  norm_counts <- dds |>
    estimateSizeFactors() |>
    vst() |>
    assay()
  
  assay(dds, 2) <- norm_counts
  assayNames(dds)[[2]] <- "vst"
  dds
}

tidy_signatures <- function(signatures_file) {
  
  # Pan B-Cell Signature: https://jitc.bmj.com/content/5/1/18.long
  
  # CD8+ T Effector Signature (cd8_rosen): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5480242/#SD1
  # Supplementary figure 6
  
  # CD8+ T (cd8_prat): https://cancerres.aacrjournals.org/content/canres/77/13/3540/F1.large.jpg?width=800&height=600&carousel=1
  
  # CD8+ T (cd8_sade): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6641984/bin/NIHMS1510803-supplement-11.xlsx
  
  # CD8+ T (cd8_fehr): https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(16)00587-0/fulltext
  
  # IFNg Signature: https://www.jci.org/articles/view/91190/table/2
  
  # Expanded Immune Signature: https://www.jci.org/articles/view/91190
  # Table 2
  
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
  
  exp_immune <- c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB")
  
  signatures <- list(
    b_cell = b_cell, cd8_rose = cd8_rosen,
    cd8_prat = cd8_prat, cd8_fehr = cd8_fehr, ifng = ifng,
    myeloid_inf = myeloid_inf, exp_immune = exp_immune
  )
  
}

run_gsva <- function(dds, signatures) {
  
  # Foresight:
  signatures$b_cell[(signatures$b_cell == "KIAA0125")] <- "FAM30A"
  signatures$ifng[(signatures$ifng == "HLA-DRA")] <- "HLA.DRA"
  signatures$exp_immune[(signatures$exp_immune == "HLA-DRA")] <- "HLA.DRA"
  signatures$exp_immune[(signatures$exp_immune == "HLA-E")] <- "HLA.E"
  
  
  genes <- unlist(signatures) |>
    as_tibble()
  no_match <- as.character(genes[!(genes$value %in% rownames(dds)), ])
  
  if (no_match == "character(0)") {
    print("All signature genes have match in dataset.")
  } else {
    warning(paste(no_match, "is in a signature but has no match in the dataset\n"))
  }
  
  scores <- gsva(
    assay(dds, 2),
    signatures,
    mx.diff = T,
    kcdf = "Gaussian"
  )
}

merge_gsva <- function(data, scores) {
  
  col_data <- data |>
    colData() |>
    as_tibble(rownames = "sample")
  
  scores <- scores |> 
    t() |>
    as_tibble(rownames = "sample")
  
  joined <- 
    inner_join(col_data, scores, by = "sample") |> 
    DataFrame()
  
  colData(data) <- joined
  
  data
}

bin_gsva <- function(dds) {
  
  col_data <- dds |>  
    colData() |> 
    as_tibble() |> 
    mutate(b_bin = if_else(b_cell > 0, "Hi", "Lo"),
           cd8_bin = if_else(cd8_rose > 0, "Hi", "Lo")) |> 
    unite(b8t, b_bin, cd8_bin, remove = FALSE) |> 
    DataFrame()
  
  colData(dds) <- col_data
  
  dds
  
}

# Plotting Helpers -------------------------------------------------------------

theme_tufte <- function(font_size = 30) {
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "#FFFFF8", color = "#CCCCCC"),
    plot.background = element_rect(fill = "#FFFFF8"),
    strip.background = element_rect(fill = "#BBBBB0", size = 0),
    strip.text = element_text(size = font_size),
    legend.background = element_rect(fill = "#FFFFF8"),
    legend.position = "top",
    legend.key = element_blank(),
    axis.line = element_line(size = 0.1, color = "#BBBBB0"),
    axis.ticks = element_line(size = 0.3, color = "#BBBBB0"),
    text = element_text(family = "Gill Sans", size = font_size)
  )
}

# Taken from user allan-cameron on StackOverflow
ggsurvplot_facet2 <- function(pval.size = 5, ...)
{
  newcall <- bquote(
    p <- p + geom_text(data = pvals.df, aes(x = pval.x, y = pval.y, 
                                            label = pval.txt), size = .(pval.size), hjust = 0)
  )
  
  body(ggsurvplot_facet)[[20]][[3]][[8]] <- newcall
  ggsurvplot_facet(...)
}


# Make Plots -------------------------------------------------------------------

## Individual Project Plots -----------------------------------------------------

make_clin_table <- function(dds, project) {
  
  proj_fig_dir <- tar_read_raw(paste0("fig_dir_", project))
  file_name <- paste0(proj_fig_dir, "clin-table.png")
  
  dds |> 
    colData() |>
    as_tibble() |> 
    dplyr::select(-(sample:cases.case_id)) |> 
    tbl_summary(by = sex) |> 
    add_overall() |> 
    as_gt() |> 
    gtsave(file_name)
  file_name
}

dense_ind <- function(data_path, x, file_name, color = NULL, facet = NULL) {
  get_gill()
  enq_fac <- enquo(facet)
  df <- data_path |> 
    read_rds() |> 
    ggplot(aes(x = {{ x }}, color = {{ color }})) + 
    geom_density() + 
    facet_grid(rows = vars(!!enq_fac)) + 
    coord_cartesian(xlim = c(-1, 1)) + 
    theme_tufte()
  
  ggsave(paste0(str_remove(data_path, "tidy-clin-dat.Rds"), file_name), 
         df, 
         width = 3, height = 3)
  paste0(str_remove(data_path, "tidy-clin-dat.Rds"), file_name)
}

surv_ind <- function(data, strata, file_name, project, 
                     legend_labs = NULL, 
                     color = NULL, 
                     linetype = NULL, 
                     facet = NULL, 
                     legend_title = NULL,
                     width = 4,
                     height = 4) {
  
  proj_fig_dir <- tar_read_raw(paste0("fig_dir_", project))
  file_name <- paste0(proj_fig_dir, file_name)
  
  data <- data |> 
    colData() |> 
    as_tibble()
  
  if (is.null(legend_labs)){
    legend_labs <- 
      setdiff(strata, facet)
    if (length(legend_labs) == 1) {
      legend_labs <- levels(data[[legend_labs]]) 
    } else {
      df <- data |> 
        dplyr::select(all_of(strata)) |>
        mutate(across(everything(), fct_drop))
      
      make_labs <- function(reduced, will_reduce) {
        if (is_tibble(reduced)) {
          reduced <- factor(unlist(reduced))
        }
        temp <- expand_grid(a = levels(reduced), b = levels(will_reduce))
        c <- unite(temp, b, a, b, sep = "/")
      }
      
      legend_labs <- purrr::reduce(df, make_labs) |>
        unlist() |>
        unname()
    }
  }
  
  if (is.null(legend_title) & is.null(color)){
    legend_title <- strata |>
      str_replace("_", " ") |> 
      str_to_title() |> 
      paste(collapse = "/")
  }
  
  ivs <- paste(strata, collapse = " + ")
  form <- as.formula(paste("Surv(new_death, death_event) ~", ivs))
  
  fit <- surv_fit(form, data = data)
  data <- as.data.frame(data)
  
  args <- list(fit = fit, data = data, pval = TRUE, pval.coord = c(0, 0.1), 
               legend.labs = legend_labs, legend.title = legend_title,
               break.time.by = 365.25, xscale = "d_m", size = 0.5, 
               pval.size = 6, font.x = 15, font.y = 15)
  
  if (is.null(facet) || data[[facet]] |> unique() |> length() == 1) {
    ggsurv <- do.call(ggsurvplot, c(args))
    ggsurv <- ggsurv$plot
  } else {
    ggsurv <- do.call(ggsurvplot_facet2, c(args, facet.by = facet, short.panel.labs = TRUE))
    width <- 9
  }
  
  ggsurv <- ggsurv +
    theme_tufte(10) + 
    theme(legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)) +
    scale_color_viridis_d(option = "plasma", end = 0.8) 
  
  agg_png(file_name, width = width, height = height, units = "in", res = 288)
  print(ggsurv)
  dev.off()
  file_name
}

density_pan_immune_cell <- function(tcga_dds = list()) {
  get_gill() 
  b_dens <- tcga_dds |> 
    map(\(x) {
      x <- x |> 
        read_rds()
      x <- cbind(exp_immune = x$exp_immune, sex = x$sex) |>
        as_tibble() |> 
        mutate(exp_immune = as.numeric(exp_immune))
    }) |> 
    bind_rows(.id = "proj")
  print(b_dens)
  b_dens <- ggplot(b_dens, aes(x = exp_immune)) +
    geom_density() + 
    facet_wrap(~proj) + 
    theme_tufte(20) + 
    theme(legend.position = "none")
  ggsave("denseplot_exp-immune.png", b_dens, width = 10, height = 10)
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
