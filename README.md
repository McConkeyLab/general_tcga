# general_tcga

This repository contains a `targets` workflow that automatically:
- Downloads TCGA RNAseq counts and clinical data
- Removes blacklisted, FFPE, or otherwise excluded samples (at the manifest level, before downloading them)
- Tidies clinical data
- Obtains HGNC gene symbols from Entrez IDs
- Coerces counts and clinical data into a `SummarizedExperiment`
- Normalizes counts using `DESeq2::vst` (supplied as the second assay slot of the resultant `SummarizedExperiment`)
- Performs `GSVA` on a set of gene signatures and appends those results to the `colData` of the `SummarizedExperiment` (This analysis may be removed later)
- Writes resultant `SummarizedExperiment` files as an `.Rds` file.

# Requirements

- R >=4.1.0 (relies on native pipe)
- The following packages:
  - `tidyverse`, `broom`, `glue`, `rvest`, `GenomicDataCommons`, `TCGAutils`, `biomaRt`, `SummarizedExperiment`, `DESeq2`, `GSVA`, `tidymodels`, `targets`, `tarchetypes`
  - You can download these quickly and easily using the following:
  ```
  install.packages("pak")
  pak::pkg_install("tidyverse", "broom", "glue", "rvest", "GenomicDataCommons", "TCGAutils", "biomaRt", "SummarizedExperiment", "DESeq2", "GSVA", "tidymodels", "targets", "tarchetypes")
  ```

# Directions
1. Create a new version control project in RStudio (git, not subversion)
2. Supply the repository URL
3. Ensure you have met the requirements above
4. Run `targets::tar_make()`
5. Wait (this process will take quite a while as it has to download all the tcga files as well as perform data manipulations. Expect over an hour.)
6. After this process has run, you will have `dds.rds` files for each TCGA project in 01_data/{project_name}. Additionally, you can load any object into your global environment by running `targets::tar_read(target_name)`. For instance, running `a <-  targets::tar_read(dds_w_bin_scores)` will assign the final `SummarizedExperiment` objects to the variable `a`.
