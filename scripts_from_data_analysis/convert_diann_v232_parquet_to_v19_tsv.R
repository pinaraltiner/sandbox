## ============================================================================
## convert_diann_v232_parquet_to_v19_tsv.R
##
## Purpose:  Read a DIA-NN v2.3.2 report in Parquet format and write a TSV file
##           that is directly consumable by the existing R post-processing
##           pipeline (exp1_DIANN_data_analysis.R, exp2_DIANN_data_analysis.R,
##           parser_func.R for Exp3).
##
## Usage:    source("convert_diann_v232_parquet_to_v19_tsv.R")
##           convert_diann_parquet(
##             parquet_path = "path/to/report.parquet",
##             output_tsv   = "path/to/report.tsv"
##           )
##
## What it does:
##   1. Reads the Parquet file using the 'arrow' package.
##   2. Renames 'Precursor.Normalised' -> 'PG.Normalised' (the only column
##      the R pipeline references that was renamed between v1.9 and v2.3.2).
##   3. Selects the columns that the downstream R scripts actually use,
##      plus a few extras that are useful for QC / debugging.
##   4. Writes a tab-delimited TSV file identical in format to v1.9 output.
##
## Requirements:
##   install.packages("arrow")   # for read_parquet()
##
## Compatibility notes:
##   - The downstream R scripts read this TSV with readr::read_tsv() and
##     reference: Run, Stripped.Sequence, Modified.Sequence, PG.Normalised,
##     PTM.Site.Confidence, Protein.Names.
##   - v2.3.2 stores Run as the filename stem (e.g. "OXPAL230121_13"),
##     same as v1.9 — no path stripping needed.
##   - Modified.Sequence uses the same UniMod notation as v1.9.
##   - The script preserves all original v2.3.2 columns by default.
##     Set keep_all_columns = FALSE to output only the minimal set.
##
## Author:   David / Claude — 2026-03-21
## Context:  Pinar Altiner phosphoproteomics benchmarking paper,
##           DIA-NN v2.3.2 re-processing (Tasks T01–T09)
## ============================================================================

library(arrow)

convert_diann_parquet <- function(parquet_path,
                                  output_tsv = NULL,
                                  keep_all_columns = TRUE) {

  ## --- Validate input -------------------------------------------------------
  if (!file.exists(parquet_path)) {
    stop("Parquet file not found: ", parquet_path)
  }

  ## --- Default output path: same directory, .tsv extension ------------------
  if (is.null(output_tsv)) {
    output_tsv <- sub("\\.parquet$", ".tsv", parquet_path)
    if (output_tsv == parquet_path) {
      output_tsv <- paste0(parquet_path, ".tsv")
    }
  }

  ## --- Read parquet ---------------------------------------------------------
  message("Reading: ", parquet_path)
  df <- arrow::read_parquet(parquet_path)
  message("  ", nrow(df), " rows, ", ncol(df), " columns")

  ## --- Column rename: Precursor.Normalised -> PG.Normalised -----------------
  ## This is the ONLY breaking rename between v1.9 and v2.3.2 that affects
  ## the downstream R pipeline. All three DIA-NN scripts (exp1, exp2, exp3)
  ## do: rename("Intensity" = "PG.Normalised")
  if ("Precursor.Normalised" %in% colnames(df)) {
    colnames(df)[colnames(df) == "Precursor.Normalised"] <- "PG.Normalised"
    message("  Renamed: Precursor.Normalised -> PG.Normalised")
  } else if ("PG.Normalised" %in% colnames(df)) {
    message("  PG.Normalised already present (v1.9 format?)")
  } else {
    warning("Neither Precursor.Normalised nor PG.Normalised found in the parquet file. ",
            "The downstream R scripts will fail at rename(\"Intensity\" = \"PG.Normalised\").")
  }

  ## --- Verify all required columns are present ------------------------------
  required_cols <- c("Run", "Stripped.Sequence", "Modified.Sequence",
                     "PG.Normalised", "PTM.Site.Confidence", "Protein.Names")

  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    warning("Missing columns required by R pipeline: ",
            paste(missing, collapse = ", "))
  } else {
    message("  All 6 required columns verified: ",
            paste(required_cols, collapse = ", "))
  }

  ## --- Optionally trim to minimal column set --------------------------------
  if (!keep_all_columns) {
    ## Minimal set: what the R scripts actually touch, plus useful extras
    minimal_cols <- c(
      ## Used directly by exp1/exp2/exp3 R scripts:
      "Run",
      "Stripped.Sequence",
      "Modified.Sequence",
      "PG.Normalised",
      "PTM.Site.Confidence",
      "Protein.Names",
      ## Useful for QC / extended analysis:
      "Precursor.Id",
      "Precursor.Charge",
      "Precursor.Quantity",
      "Precursor.Mz",
      "Protein.Ids",
      "Protein.Group",
      "Genes",
      "Q.Value",
      "Global.Q.Value",
      "Peptidoform.Q.Value",
      "Global.Peptidoform.Q.Value",
      "Lib.PTM.Site.Confidence",
      "Site.Occupancy.Probabilities",
      "Protein.Sites",
      "PEP",
      "RT",
      "IM"
    )
    keep <- intersect(minimal_cols, colnames(df))
    df <- df[, keep]
    message("  Trimmed to ", ncol(df), " columns (minimal mode)")
  }

  ## --- Write TSV ------------------------------------------------------------
  message("Writing: ", output_tsv)
  write.table(df,
              file      = output_tsv,
              sep       = "\t",
              row.names = FALSE,
              quote     = FALSE,
              na        = "")  ## DIA-NN v1.9 TSV uses empty string for NA

  message("Done. ", nrow(df), " rows written to ", output_tsv)

  invisible(output_tsv)
}


## ============================================================================
## Batch conversion helper
## ============================================================================

#' Convert all DIA-NN v2.3.2 parquet reports in a directory
#'
#' @param dir_path       Directory containing .parquet files
#' @param pattern        Regex to match parquet files (default: report parquets)
#' @param keep_all_columns  Keep all columns (TRUE) or minimal set (FALSE)
#' @param recursive      Search subdirectories
#'
#' @return Character vector of output TSV paths (invisible)

batch_convert_diann_parquet <- function(dir_path,
                                        pattern = "^report.*\\.parquet$",
                                        keep_all_columns = TRUE,
                                        recursive = TRUE) {

  parquet_files <- list.files(dir_path,
                              pattern    = pattern,
                              full.names = TRUE,
                              recursive  = recursive)

  ## Exclude speclib and site_report parquets — only convert main reports
  parquet_files <- parquet_files[!grepl("report-lib|site_report|predicted\\.speclib",
                                        parquet_files)]

  if (length(parquet_files) == 0) {
    message("No matching parquet files found in: ", dir_path)
    return(invisible(character(0)))
  }

  message("Found ", length(parquet_files), " parquet report(s) to convert:\n",
          paste(" ", parquet_files, collapse = "\n"))

  output_paths <- character(length(parquet_files))
  for (i in seq_along(parquet_files)) {
    message("\n--- [", i, "/", length(parquet_files), "] ---")
    output_paths[i] <- convert_diann_parquet(parquet_files[i],
                                              keep_all_columns = keep_all_columns)
  }

  message("\n=== Batch conversion complete: ", length(output_paths), " files ===")
  invisible(output_paths)
}
