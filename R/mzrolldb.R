#' Table of samples from file
#'
#' @param mzrolldb_file_path: file path to mzrolldb file
#'
#' @return samples table from mzrolldb file
#'
#' @export
PDB_sample_list <- function(mzrolldb_file_path) {
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = mzrolldb_file_path)

  samples <- dplyr::tbl(con, "samples") %>%
    dplyr::collect() %>%
    dplyr::arrange(sampleId)

  DBI::dbDisconnect(conn = con)

  return(samples)
}

#' Table of peaks from file
#'
#' @param mzrolldb_file_path: file path to mzrolldb file
#'
#' @return peakgroups table from mzrolldb file
#'
#' @export
PDB_peaks <- function(mzrolldb_file_path) {
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = mzrolldb_file_path)

  cols_to_cast <- c(
    # Identifiers and Integer Positions
    "peakId", "groupId", "sampleId", "pos", "minpos", "maxpos",
    "scan", "minscan", "maxscan", "width", "label", "fromBlankSample",

    # Retention Time and m/z (Real/Numeric)
    "rt", "rtmin", "rtmax", "mzmin", "mzmax", "peakMz", "medianMz", "baseMz",

    # Intensities and Areas
    "peakArea", "peakAreaCorrected", "peakAreaTop", "peakAreaFractional",
    "peakRank", "peakIntensity", "peakBaseLineLevel", "smoothedIntensity",
    "smoothedPeakArea", "smoothedPeakAreaCorrected", "smoothedPeakAreaTop",

    # Quality and Fit Stats
    "quality", "gaussFitSigma", "gaussFitR2", "noNoiseObs", "noNoiseFraction",
    "symmetry", "signalBaselineRatio", "groupOverlap", "groupOverlapFrac",
    "localMaxFlag", "smoothedSignalBaselineRatio",

    # FWHM Specific Columns
    "minPosFWHM", "maxPosFWHM", "minScanFWHM", "maxScanFWHM",
    "rtminFWHM", "rtmaxFWHM", "peakAreaFWHM", "smoothedPeakAreaFWHM"
  )

  peakgroups <- dplyr::tbl(con, "peaks") %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(cols_to_cast), as.numeric)) %>%
    dplyr::collect() %>%
    dplyr::arrange(peakId)

  DBI::dbDisconnect(conn = con)

  return(peakgroups)
}

#' Table of peakgroups from file
#'
#' @param mzrolldb_file_path: file path to mzrolldb file
#'
#' @return peakgroups table from mzrolldb file
#'
#' @export
PDB_peakgroups <- function(mzrolldb_file_path) {
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = mzrolldb_file_path)

  peakgroups <- dplyr::tbl(con, "peakgroups") %>%
    dplyr::collect() %>%
    dplyr::arrange(groupId)

  DBI::dbDisconnect(conn = con)

  return(peakgroups)
}
