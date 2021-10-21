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

get_clin <- function(tcga_data_dir, tcga_proj) {
  
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

get_latest <- function(data) {
  rowwise(data) |>
    mutate(across(everything(), as.numeric)) |>
    mutate(collapsed = if_else(!is.infinite(max(c_across(), na.rm = T)),
                               max(c_across(), na.rm = T),
                               NA_real_
    )) |>
    dplyr::select(collapsed)
}

tidy_clin <- function(clin_path) {
  
  clinical <- read_tsv(clin_path, col_names = F, show_col_types = FALSE) |>
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
                  path_stage = matches("^patient.stage_event.pathologic_stage$"),
                  clin_stage = matches("^patient.stage_event.clinical_stage$"),
                  pathologic_tnm_t = matches("^patient.stage_event.tnm_categories.pathologic_categories.pathologic_t$"),
                  pathologic_tnm_n = matches("^patient.stage_event.tnm_categories.pathologic_categories.pathologic_n$"),
                  pathologic_tnm_m = matches("^patient.stage_event.tnm_categories.pathologic_categories.pathologic_m$"),
                  clinical_tnm_t = matches("^patient.stage_event.tnm_categories.clinical_categories.clinical_t$"),
                  clinical_tnm_n = matches("^patient.stage_event.tnm_categories.clinical_categories.clinical_n$"),
                  clinical_tnm_m = matches("^patient.stage_event.tnm_categories.clinical_categories.clinical_m$"),
                  grade = matches("^patient.neoplasm_histologic_grade$"),
                  pack_years = matches("^patient.number_pack_years_smoked$"),
                  hpv_status = matches("^patient.hpv_test_results.hpv_test_result.hpv_status$"),
                  h_pylori_status = matches("^patient.h_pylori_infection$")) |>
    mutate(across(contains("pack_years"), as.numeric)) |> 
    mutate(across(contains("hpv"), ~ fct_relevel(., "negative"))) |> 
    mutate(across(contains("tnm"), ~ str_remove(., "(?<=[:digit:])[:alpha:]$"))) |> 
    mutate(across(matches("path_stage|clin_stage"), ~str_remove(., "[^0ivs]$"))) |> 
    mutate(across(contains("tnm_n"), ~ str_replace(., "[1234]$", "\\+"))) |> 
    mutate(across(contains("tnm_n"), ~ fct_relevel(., "n+")))
  
  # Create Survival Time Column
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
      follow_up_time = if_else(is.na(death_days), followUp_days, death_days),
      death = if_else(patient.vital_status == "dead" | !is.na(death_days), 1, 0)
    ) |>
    dplyr::select(-contains("days_to_death"),
                  -contains("days_to_last_followup"),
                  -patient.vital_status,
                  -death_days,
                  -followUp_days) |> 
    dplyr::filter((follow_up_time >= 0) & !is.na(follow_up_time)) |>
    dplyr::filter((follow_up_time > 0) | (death == 0)) |> 
    # For reasoning as to the above filtering rule, see
    # https://www.graphpad.com/support/faq/events-deaths-at-time-zero-in-survival-analysis/
    dplyr::filter(!is.na(sex)) |> 
    mutate(sex = factor(sex, levels = c("male", "female"), labels = c("M", "F")),
           age = as.numeric(age)) |> 
    mutate(across(matches("grade|race"), fct_rev))
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

download_tcga_data <- function(manifest) {
  manifest <- manifest |> 
    dplyr::select(id, path, submitter_id, short_id, cases.case_id)
  gdc_set_cache("./01_data/00_gdcdata", create_without_asking = T)
  lapply(manifest$id, gdcdata)
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
  
  sample_table <- data.frame(sampleName = tumor_data$id,
                             fileName = tumor_data$path)
  
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table,
                                    directory = "./01_data/00_gdcdata",
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

normalize <- function(dds, ids) {
  
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

write_norm <- function(dds, project) {
  write_rds(dds, paste0("./01_data/", project, "/norm.Rds"))
}

tidy_signatures <- function() {
  
  # Pan B-Cell Signature: https://jitc.bmj.com/content/5/1/18.long
  
  # CD8+ T Effector Signature: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5480242/#SD1
  # Supplementary figure 6
  
  # Expanded Immune Signature: https://www.jci.org/articles/view/91190
  # Table 2
  
  b_cell <- c("BLK", "CD19", "FCRL2", "MS4A1", "KIAA0125", "TNFRSF17", "TCL1A",
              "SPIB", "PNOC")
  cd8 <- c("CD8A", "GZMA", "GZMB", "IFNG", "CXCL9", "CXCL10", "PRF1", "TBX21")
  exp_immune <- c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", 
                  "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", 
                  "TAGAP", "CXCL10", "STAT1", "GZMB")
  
  signatures <- list(b_cell = b_cell, cd8 = cd8, exp_immune = exp_immune)
  
}

run_gsva <- function(dds, signatures) {
  
  # Foresight:
  signatures$b_cell[(signatures$b_cell == "KIAA0125")] <- "FAM30A"
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
           b_bin = factor(b_bin, levels = c("Lo", "Hi")),
           cd8_bin = if_else(cd8 > 0, "Hi", "Lo"),
           cd8_bin = factor(cd8_bin, levels = c("Lo", "Hi")),
           imm_bin = if_else(exp_immune > 0, "Hi", "Lo"),
           imm_bin = factor(imm_bin, levels = c("Lo", "Hi"))) |>
    DataFrame()
  
  colData(dds) <- col_data
  
  dds
  
}

# Survival ---------------------------------------------------------------------

tidy_for_survival <- function(dds, project) {
  
  data <- dds |> 
    colData() |> 
    as_tibble() |> 
    dplyr::select(-c(sample:cases.case_id)) # ID columns, not used for survival
  
  # Remove columns that are all NA
  na_sums <- apply(data, 2, \(x) is.na(x) |> sum())
  just_nas <- which(na_sums == nrow(data))
  if(length(just_nas) == 0) {
    no_nas <- data
  } else {
    no_nas <- data[,-just_nas]
  }
  
  # Remove TNM 'x'
  no_x <- no_nas |>
    mutate(across(contains("tnm"), ~ str_replace(.x, ".*x$", NA_character_))) |> 
    mutate(across(contains("tnm_n"), ~ if_else(.x != "n0", "n+", .x))) |> 
    mutate(across(contains("tnm_n"), ~ fct_relevel(.x, "n0"))) |> 
    mutate(across(matches("^grade$"), ~ str_replace(.x, ".*x$", NA_character_))) |> 
    mutate(across(matches("^grade$"), ~ fct_relevel(.x, "low")))
  
  data <- no_x |> 
    dplyr::select(-matches("days_to_collection|anatomic")) |> 
    dplyr::select(-b_cell, -cd8, -exp_immune) # Just using discretized 
  
  if (project %in% c("coad", "hnsc", "lihc", "luad")) {
    data <- data |> 
      mutate(race = if_else(race %in% c("asian", "american indian or alaska native"), 
                            "asian, american indian or alaska native", 
                            as.character(race))) #if not done, causes issues with Anova in univariate
  }
  
  # Project specific cases
  if (project == "blca") {
    data <- data |> 
      dplyr::filter(path_stage != "stage i") |> 
      dplyr::filter(!(pathologic_tnm_t %in% c("t0", "t1"))) |> 
      mutate(grade = fct_rev(grade))
  }
  
  if (project == "coad") {
    data <- data |> 
      mutate(pathologic_tnm_t = if_else(pathologic_tnm_t %in% c("t1", "t2"), "t1/2", pathologic_tnm_t)) |> 
      dplyr::filter(pathologic_tnm_t != "tis") |> 
      mutate(pathologic_tnm_t = fct_relevel(pathologic_tnm_t, "t1/2"))
  }
  
  if (project == "hnsc") {
    data <- data |> 
      mutate(path_stage = if_else(path_stage %in% c("stage 0", "stage i", "stage ii"), "stage 0/i/ii", path_stage),
             clin_stage = if_else(clin_stage %in% c("stage 0", "stage i", "stage ii"), "stage 0/i/ii", clin_stage),
             pathologic_tnm_t = if_else(pathologic_tnm_t %in% c("t0", "t1", "t2"), "t0/1/2", pathologic_tnm_t), 
             clinical_tnm_t = if_else(clinical_tnm_t %in% c("t1", "t2"), "t1/2", clinical_tnm_t),
             grade = if_else(grade %in% c("g3", "g4"), "g3/4", as.character(grade)),
             hpv_status = if_else(hpv_status %in% c("negative", "positive"), as.character(hpv_status), NA_character_),
             hpv_status = fct_relevel(hpv_status, "negative"))
  }
  
  if (project == "kirc") {
    data <- data |> 
      mutate(pathologic_tnm_t = if_else(pathologic_tnm_t %in% c("t3", "t4"), "t3/4", pathologic_tnm_t),
             grade = if_else(grade %in% c("g1", "g2"), "g1/2", as.character(grade))) |> 
      dplyr::select(-clinical_tnm_m)
  }
  
  if (project == "lihc") {
    data <- data |> 
      mutate(path_stage = if_else(path_stage %in% c("stage iii", "stage iv"), "stage iii/iv", path_stage),
             path_stage = fct_relevel(path_stage, "stage iii/iv", after = Inf),
             pathologic_tnm_t = if_else(pathologic_tnm_t %in% c("t3", "t4"), "t3/4", pathologic_tnm_t),
             pathologic_tnm_t = fct_relevel(pathologic_tnm_t, "t3/4", after = Inf)) #if not done, causes issues with Anova in univariate
  } 
  
  if (project == "lusc") {
    data <- data |> 
      mutate(path_stage = if_else(path_stage %in% c("stage iii", "stage iv"), "stage iii/iv", path_stage))
  }

  if (project == "skcm") {
    data <- data |> 
      dplyr::filter(!(path_stage %in% c("stage 0", "i/ii nos"))) |> # rm'ing stage 0  removes the one AA individual
      dplyr::filter(!(pathologic_tnm_t %in% c("t0", "tis"))) |> 
      mutate(path_stage = fct_drop(path_stage))
  }
  
  if (project == "stad") {
    data <- data |> 
      mutate(race = if_else(race %in% c("asian", "native hawaiian or other pacific islander"), 
                            "asian, native hawaiian or other pacific islander", as.character(race))) |> #if not done, causes issues with Anova in univariate
      dplyr::select(-hpv_status) # Too few to do any surv on
  }
  
  data <- data |> 
    mutate(race = race |> fct_drop() |> fct_relevel("white"))
  
  data
  
}

run_univariate <- function(data, stratum, sex_arg = "all") {
  
  if (sex_arg != "all") {
    data <- dplyr::filter(data, sex == sex_arg)
  }
  
  #Checks for factors with one level
  if (data[[stratum]] |> as.factor() |> droplevels() |> levels() |> length() <= 1) {
    return()
  }
  
  cox_model <- function(df) {
    as.formula(paste("Surv(follow_up_time, death) ~ ", stratum)) |> 
      coxph(data = df)
  }
  

  
  tibble(data = list(data)) |> 
    mutate(fit_cox = map(data, cox_model),
           fit_tidy_cox = map(fit_cox, tidy, conf.int = TRUE, exponentiate = TRUE),
           fit_tidy_anova = map(fit_cox, car::Anova, type = "III"),
           fit_tidy_anova = map(fit_tidy_anova, tidy)) |> 
    dplyr::select(-data)
}

run_all_uni_combos <- function(data, project) {
  strata <- colnames(data)
  strata <- strata[!strata %in% c("follow_up_time", "death")]
  
  eg <- expand_grid(sex_arg = c("all", "M", "F"), stratum = strata)

  tibble(data = list(data), eg) |> 
    mutate(project = project,
           res = pmap(list(data, stratum, sex_arg), run_univariate)) |> 
    dplyr::select(-data) |> 
    unnest(res) |> 
    unnest(cols = fit_tidy_cox, names_sep = "_") |> 
    unnest(cols = fit_tidy_anova, names_sep = "_") |> 
    mutate(base_level = map(fit_cox, function(x) {
      x <- x$xlevels |> unlist() |> unname()
      x <- x[1]
      x <- if_else(is.null(x), NA_character_, x)
    })) |> 
    unnest(cols = base_level)
}

get_multivariable_names <- function(univariate, project) {
  univariate <- univariate |> 
    dplyr::filter(sex_arg == "all",
                  fit_tidy_anova_p.value <= 0.05)
  multi_names <- union(univariate$stratum, c("age", "path_stage"))
  # Since we're going with path_stage, remove the pathologic_tnm
  multi_names <- multi_names[!str_detect(multi_names, "_tnm")]
  # These factors will be added in individually later and must be removed so
  # they can be added in a controlled fashion
  multi_names <- multi_names[!(multi_names %in% c("sex", "cd8", "b_cell", "exp_immune", 
                                                  "cd8_bin", "b_bin", "imm_bin"))]
  
  if (project == "lusc") {
    multi_names <- c(multi_names, "pack_years")
  }
  
  multi_names
}

run_multivariable <- function(data, multivariable_names, sex_arg = "all", signature = "none") {
  
  if (sex_arg != "all") {
    data <- dplyr::filter(data, sex == sex_arg)
  }
  
  if (signature != "none") {
    multivariable_names <- c(multivariable_names, signature)
  }

  naming_function <- function(var, lvl, ordinal = FALSE, sep = "|") {
    paste0(var, "|", lvl)
  }
  
  rec <- recipe(data) |>   
    step_select(all_of(multivariable_names), follow_up_time, death) |> 
    step_dummy(all_nominal(), naming = naming_function) 
  
  prepared_data <- rec |> 
    prep() |> 
    bake(new_data = NULL)
  
  fit_mod <- function(df) {
    coxph(Surv(follow_up_time, death) ~ ., data = df)
  }
  
  tibble(data = list(prepared_data)) |> 
    mutate(fit = map(data, fit_mod),
           tidy_coxph = map(fit, tidy, conf.int = TRUE, exponentiate = TRUE),
           tidy_anova = map(fit, \(x) x |> car::Anova(type = "III") |> tidy())) |> 
    dplyr::select(-data)
}

run_all_multi_combos <- function(data, names, project) {
  
  eg <- expand_grid(sex_arg = c("all", "M", "F"), signature = c("none", "cd8_bin", "b_bin", "imm_bin"))
  
  tibble(data = list(data), names = list(names), eg) |> 
    rowwise() |> 
    mutate(res = list(run_multivariable(data, names, sex_arg, signature))) |> 
    dplyr::select(-data) |> 
    mutate(project = project) |> 
    unnest(res) 
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
    text = element_text(family = "Gill Sans MT", size = font_size)
  )
}

# Taken from user allan-cameron on StackOverflow
ggsurvplot_facet2 <- function(pval.size = 5, ...) {
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
  
  clin <- dds |> 
    colData() |>
    as_tibble() |> 
    dplyr::select(-(sample:cases.case_id))
  
  how_many_nas <- apply(clin, 2, \(x) x |> is.na() |> sum())
  all_nas <- how_many_nas == nrow(clin)
  
  clin <- clin[,!all_nas]
  
  clin <- dplyr::select(clin, -cd8, -b_cell, -exp_immune)
  
  colnames(clin) <- colnames(clin) |> 
    str_replace_all("_", " ") |> 
    str_to_title() |> 
    str_replace_all("Cd8", "CD8+") |> 
    str_replace_all("Bin", "Signature") |> 
    str_replace_all("Tnm", "TNM")
  
  clin |> 
    tbl_summary(by = Sex) |> 
    add_overall(last = TRUE) |> 
    as_gt() |> 
    gtsave(file_name)
  file_name
}

make_univariate_plot <- function(univariate_data, project) {
  proj_fig_dir <- tar_read_raw(paste0("fig_dir_", project))
  file_name <- paste0(proj_fig_dir, "univariate-plot.png")

  b <- univariate_data |> 
    rowwise() |>
    mutate(level = if_else(fit_tidy_cox_term != stratum, 
                           str_remove(fit_tidy_cox_term, paste0("^", stratum)),
                           fit_tidy_cox_term),
           subtitle = sex_arg |> 
             str_replace("^M$", "Males") |> 
             str_replace("^F$", "Females") |> 
             str_replace("^all$", "All"),
           est_label = paste(round(fit_tidy_cox_conf.low, 2), round(fit_tidy_cox_estimate, 2), round(fit_tidy_cox_conf.high, 2))) |> 
    relocate(level, .before = fit_tidy_cox_term) |> 
    mutate(plot_conf_low = if_else(fit_tidy_cox_conf.low < 0.03, 0.03, fit_tidy_cox_conf.low),
           plot_conf_high = if_else(fit_tidy_cox_conf.high > 25, 25, fit_tidy_cox_conf.high))

  
  # all
  all <- b |>
    dplyr::filter(sex_arg == "all") |> 
    dplyr::filter(fit_tidy_anova_term != "NULL") |> 
    group_by(stratum) |> 
    mutate(mean_est = mean(fit_tidy_cox_estimate)) |> 
    arrange(fit_tidy_cox_estimate) |> 
    ungroup() |> 
    mutate(
      anova_stars = case_when(fit_tidy_anova_p.value < 0.001 ~ "***",
                              fit_tidy_anova_p.value < 0.01 ~ "**",
                              fit_tidy_anova_p.value < 0.05 ~ "*",
                              TRUE ~ ""),
      indiv_stars = case_when(fit_tidy_cox_p.value < 0.001 ~ "***",
                              fit_tidy_cox_p.value < 0.01 ~ "**",
                              fit_tidy_cox_p.value < 0.05 ~ "*",
                              TRUE ~ ""),
      y = if_else(is.na(base_level),
                  level,
                  paste0(level, " vs ", base_level)),
      y = y |> 
        str_replace_all("_", " ") |> 
        str_to_title() |> 
        str_replace_all("(\\bI[iv]*\\b)", str_to_upper) |> 
        str_replace("(\\bVs\\b)", str_to_lower),
      y = paste(indiv_stars, y),
      stratum = paste(stratum, anova_stars),
      stratum = stratum |> 
        str_replace_all("_", " ") |> 
        str_to_title() |> 
        str_replace_all("Cd8", "CD8+") |> 
        str_replace_all("\\bB\\b", "B-cell") |> 
        str_replace_all("\\bBin\\b", "Signature") |> 
        str_replace_all("Tnm", "TNM") |> 
        fct_inorder())
  height <- (all$y |> length())/3
  all <- all |>  
    ggplot(aes(x = fit_tidy_cox_estimate, y = y)) + 
    geom_vline(xintercept = 1, color = "#FF0000", size = 0.2) + 
    geom_linerange(aes(xmin = plot_conf_low, xmax = plot_conf_high)) +
    geom_point() + 
    facet_grid(rows = "stratum", shrink = T, scales = "free", space = "free") +
    bladdr::theme_tufte(10) +
    scale_x_log10(breaks = c(0.05, 0.1, 0.25, 0.5, 2, 4, 10, 20), labels = c(0.05, 0.1, 0.25, 0.5, 2, 4, 10, 20)) +
    theme(panel.grid.major.x = element_line(color = "#CCCCCC"),
          strip.placement = "inside",
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10, face = "bold"),
          axis.text.y = element_text(size = 9),
          axis.title = element_blank(),
          axis.ticks = element_blank()) + 
    coord_cartesian(xlim = c(0.05, 20))
  

  agg_png(file_name, width = 9, height = height,  units = "in", res = 288)
  print(all)
  dev.off()
  
  file_name
  
}

make_multivariable_plot <- function(multi_data, project) {
  proj_fig_dir <- tar_read_raw(paste0("fig_dir_", project))
  file_name <- paste0(proj_fig_dir, "multivariable-plot.png")
  
  # Temp
  multi_data <- tar_read(multivariable_hnsc) |> 
    dplyr::select(-fit) |> 
    unnest(tidy_coxph) |> 
    dplyr::mutate(term = str_remove_all(term, "`")) |> 
    separate(term, into = c("stratum", "level"), sep = "\\|")
  
  b <- multi_data |> 
    unnest(cols = names) |> 
    rowwise() |>
    mutate(level = if_else(tidy_coxph_term != names, 
                           str_remove(tidy_coxph_term, paste0("^", names)),
                           tidy_coxph_term),
           subtitle = sex_arg |> 
             str_replace("^M$", "Males") |> 
             str_replace("^F$", "Females") |> 
             str_replace("^all$", "All"),
           est_label = paste(round(tidy_coxph_conf.low, 2), round(tidy_coxph_estimate, 2), round(tidy_coxph_conf.high, 2))) |> 
    relocate(level, .before = tidy_coxph_term) |> 
    mutate(plot_conf_low = if_else(tidy_coxph_conf.low < 0.03, 0.03, tidy_coxph_conf.low),
           plot_conf_high = if_else(tidy_coxph_conf.high > 25, 25, tidy_coxph_conf.high))
  
  
  # all
  all <- b |>
    dplyr::filter(sex_arg == "all") |> 
    dplyr::filter(tidy_anova_term != "NULL") |> 
    group_by(names) |> 
    mutate(mean_est = mean(tidy_coxph_estimate)) |> 
    arrange(tidy_coxph_estimate) |> 
    ungroup() |> 
    mutate(
      anova_stars = case_when(tidy_anova_p.value < 0.001 ~ "***",
                              tidy_anova_p.value < 0.01 ~ "**",
                              tidy_anova_p.value < 0.05 ~ "*",
                              TRUE ~ ""),
      indiv_stars = case_when(tidy_coxph_p.value < 0.001 ~ "***",
                              tidy_coxph_p.value < 0.01 ~ "**",
                              tidy_coxph_p.value < 0.05 ~ "*",
                              TRUE ~ ""),
      y = if_else(is.na(base_level),
                  level,
                  paste0(level, " vs ", base_level)),
      y = y |> 
        str_replace_all("_", " ") |> 
        str_to_title() |> 
        str_replace_all("(\\bI[iv]*\\b)", str_to_upper) |> 
        str_replace("(\\bVs\\b)", str_to_lower),
      y = paste(indiv_stars, y),
      stratum = paste(names, anova_stars),
      stratum = stratum |> 
        str_replace_all("_", " ") |> 
        str_to_title() |> 
        str_replace_all("Cd8", "CD8+") |> 
        str_replace_all("\\bB\\b", "B-cell") |> 
        str_replace_all("\\bBin\\b", "Signature") |> 
        str_replace_all("Tnm", "TNM") |> 
        fct_inorder())
  height <- (all$y |> length())/3
  all |>  
    ggplot(aes(x = tidy_coxph_estimate, y = y)) + 
    geom_vline(xintercept = 1, color = "#FF0000", size = 0.2) + 
    geom_linerange(aes(xmin = plot_conf_low, xmax = plot_conf_high)) +
    geom_point() + 
    facet_grid(rows = "stratum", shrink = T, scales = "free", space = "free") +
    bladdr::theme_tufte(10) +
    scale_x_log10(breaks = c(0.05, 0.1, 0.25, 0.5, 2, 4, 10, 20), labels = c(0.05, 0.1, 0.25, 0.5, 2, 4, 10, 20)) +
    theme(panel.grid.major.x = element_line(color = "#CCCCCC"),
          strip.placement = "inside",
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10, face = "bold"),
          axis.text.y = element_text(size = 9),
          axis.title = element_blank(),
          axis.ticks = element_blank()) + 
    coord_cartesian(xlim = c(0.05, 20))
  
  
  agg_png(file_name, width = 9, height = height,  units = "in", res = 288)
  print(all)
  dev.off()
  
  file_name
  
}

dense_ind <- function(data, x, file_name, project, color = NULL, 
                      facet = NULL, width = 3, height = 3) {
  
  proj_fig_dir <- tar_read_raw(paste0("fig_dir_", project))
  file_name <- paste0(proj_fig_dir, file_name)
  
  enq_fac <- enquo(facet)
  
  this_plot <- data |> 
    colData() |> 
    as_tibble() |> 
    ggplot(aes(x = .data[[x]], color = .data[[color]])) + 
    geom_density() + 
    facet_grid(rows = vars(!!enq_fac)) + 
    coord_cartesian(xlim = c(-1, 1)) + 
    theme_tufte(10)
  
  agg_png(file_name, width = width, height = height,  units = "in", res = 288)
  print(this_plot)
  dev.off()
  
  file_name
}

dens_ind_all <- function(data, xs, project, color) {

  expand_grid(data = list(data), xs) |> 
    mutate(file_name = paste0("density_", xs, ".png"),
           project = project,
           color = color,
           file_path = pmap(list(data = data, x = xs, file_name = file_name, project = project, color = color), dense_ind)) |> 
    pull(file_path) |> 
    unlist()
}

surv_ind <- function(data, strata, file_name, project, 
                     legend_labs = NA, 
                     color = NA, 
                     linetype = NA, 
                     facet = NA, 
                     legend_title = NA,
                     width = 4,
                     height = 4) {
  
  proj_fig_dir <- tar_read_raw(paste0("fig_dir_", project))
  file_name <- paste0(proj_fig_dir, file_name)
  
  data <- data |> 
    colData() |> 
    as_tibble()
  
  if (is.na(legend_labs)){
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
  
  if (is.na(legend_title) & is.na(color)){
    legend_title <- strata |>
      str_replace("_", " ") |> 
      str_to_title() |> 
      paste(collapse = "/")
  }
  
  ivs <- paste(strata, collapse = " + ")
  form <- as.formula(paste("Surv(follow_up_time, death) ~", ivs))
  
  fit <- surv_fit(form, data = data)
  data <- as.data.frame(data)
  
  args <- list(fit = fit, data = data, pval = TRUE, pval.coord = c(0, 0.1), 
               legend.labs = legend_labs, legend.title = legend_title,
               break.time.by = 365.25, xscale = "d_m", size = 0.5, 
               pval.size = 6, font.x = 15, font.y = 15, font.tickslab = 7)
  if (is.na(facet) || data[[facet]] |> unique() |> length() == 1) {
    ggsurv <- do.call(ggsurvplot, c(args))
    ggsurv <- ggsurv$plot
  } else {
    ggsurv <- do.call(ggsurvplot_facet2, c(args, facet.by = facet, short.panel.labs = TRUE))
    width <- 8
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

surv_ind_all <- function(data, strata, project, facet) {
  expand_grid(data = list(data), strata, facet) |> 
    mutate(
      
      file_name = if_else(is.na(facet),
                          paste("survival", strata, sep = "_") |> paste0(".png"),
                          paste("survival", strata, facet, sep = "_") |> paste0(".png")),
      project = project,
      file_path = pmap(list(data = data, 
                            strata = strata, 
                            file_name = file_name,
                            facet = facet, 
                            project = project), surv_ind)) |> 
    pull(file_path) |> 
    unlist()
}

## Aggregate Project Plots -----------------------------------------------------

make_hr_plot <- function(data) {
  hr_plot <- data |> 
    dplyr::filter(str_detect(term, "cd8_bin|b_bin|imm_bin")) |> 
    mutate(project = toupper(project),
           sex_arg = factor(sex_arg, levels = c("F", "M", "all")),
           term = case_when(term == "b_binHi" ~ "B-cell Signature",
                            term == "cd8_binHi" ~ "CD8+ T-cell Signature",
                            term == "imm_binHi" ~ "Pan-Immune Signature")) |> 
    ggplot(aes(x = estimate, y = sex_arg, color = sex_arg)) +
    geom_vline(xintercept = 1, alpha = 0.5) + 
    scale_color_manual(values = c("#F8B7CD", "#0671B7", "black")) + 
    geom_linerange(aes(xmin = conf.low, xmax = conf.high)) + 
    geom_point() + 
    facet_grid(project~term) + 
    theme_tufte(10) + 
    labs(x = "Hazard Ratio") + 
    theme(legend.position = "none",
          axis.title.y = element_blank())
  agg_png("./00_common/hr-plot.png", width = 6, height = 8, units = "in", res = 288)
  print(hr_plot)
  dev.off()
  "./00_common/hr-plot.png"
}

make_hr_plot_multi <- function(data) {
  hr_plot <- data |> 
    dplyr::filter(str_detect(tidy_coxph_term, "cd8_bin|b_bin|imm_bin")) |> 
    mutate(project = toupper(project),
           sex = factor(sex_arg, levels = c("F", "M", "all")),
           term = case_when(tidy_coxph_term == "b_binHi" ~ "B-cell Signature",
                            tidy_coxph_term == "cd8_binHi" ~ "CD8+ T-cell Signature",
                            tidy_coxph_term == "imm_binHi" ~ "Pan-Immune Signature")) |> 
    ggplot(aes(x = tidy_coxph_estimate, y = sex, color = sex)) +
    geom_vline(xintercept = 1, alpha = 0.5) + 
    scale_color_manual(values = c("#F8B7CD", "#0671B7", "black")) + 
    geom_linerange(aes(xmin = tidy_coxph_conf.low, xmax = tidy_coxph_conf.high)) + 
    geom_point() + 
    facet_grid(project~term) + 
    theme_tufte(10) + 
    labs(x = "Hazard Ratio") + 
    theme(legend.position = "none",
          axis.title.y = element_blank())
  
  agg_png("./00_common/hr-plot_multi.png", width = 6, height = 8, units = "in", res = 288)
  print(hr_plot)
  dev.off()
  "./00_common/hr-plot_multi.png"
}