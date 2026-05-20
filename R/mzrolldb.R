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
    "peakId",
    "groupId",
    "sampleId",
    "pos",
    "minpos",
    "maxpos",
    "scan",
    "minscan",
    "maxscan",
    "width",
    "label",
    "fromBlankSample",

    # Retention Time and m/z (Real/Numeric)
    "rt",
    "rtmin",
    "rtmax",
    "mzmin",
    "mzmax",
    "peakMz",
    "medianMz",
    "baseMz",

    # Intensities and Areas
    "peakArea",
    "peakAreaCorrected",
    "peakAreaTop",
    "peakAreaFractional",
    "peakRank",
    "peakIntensity",
    "peakBaseLineLevel",
    "smoothedIntensity",
    "smoothedPeakArea",
    "smoothedPeakAreaCorrected",
    "smoothedPeakAreaTop",

    # Quality and Fit Stats
    "quality",
    "gaussFitSigma",
    "gaussFitR2",
    "noNoiseObs",
    "noNoiseFraction",
    "symmetry",
    "signalBaselineRatio",
    "groupOverlap",
    "groupOverlapFrac",
    "localMaxFlag",
    "smoothedSignalBaselineRatio",

    # FWHM Specific Columns
    "minPosFWHM",
    "maxPosFWHM",
    "minScanFWHM",
    "maxScanFWHM",
    "rtminFWHM",
    "rtmaxFWHM",
    "peakAreaFWHM",
    "smoothedPeakAreaFWHM"
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

#' Label mzrollDB file
#'
#' @description
#' Based on a \code{top_hits} results tibble, update an mzrollDB file
#' by labeling peak groups appropriately.
#'
#' Labeling follows the following rules:
#' label 'c' only - incorporation (only for parents), no sig isotopes
#' label 'a' incorporation + score diff > 1 - used for individual isotopes
#' label 'ca' sig isotopes and parents
#'
#' @param mzrolldb_file_path file path to mzrolldb file
#' @param top_hits tibble consisting of at least columns \code{groupId}, \code{isotope}, and \code{groupRank} columns,
#' where \code{groupRank} refers to an isotope-specific score
#' @param sig_isotopic_incorporation_scores tibble consisting of at least columns \code{groupId}, \code{groupRank}, where \code{groupRank}
#' refers only to a score associated with evidence of isotopic incorporation.
#'
#' @export
label_isotopes_by_top_hits <- function(
  mzrolldb_file_path,
  top_hits,
  sig_isotopic_incorporation_scores
) {
  all_groups <- clamshell::PDB_peakgroups(mzrolldb_file_path)

  sig_matches <- top_hits %>%
    dplyr::select(groupId, isotope, groupRank) %>%
    dplyr::mutate(label_add = "a") %>%
    dplyr::distinct() %>%
    dplyr::mutate(groupId = as.numeric(groupId))

  sig_parents <- sig_isotopic_incorporation_scores %>%
    dplyr::select(groupId, groupRank) %>%
    dplyr::group_by(groupId) %>%
    dplyr::mutate(groupRank = max(groupRank)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      label_add = ifelse(groupId %in% sig_matches$groupId, "ca", "c")
    )

  # [1] Significant isotopes
  sig_isotopes <- all_groups %>%
    dplyr::select(-groupRank) %>%
    dplyr::inner_join(
      sig_matches,
      by = c("parentGroupId" = "groupId", "tagString" = "isotope")
    ) %>%
    dplyr::mutate(label = paste0(label, label_add)) %>%
    dplyr::select(dplyr::all_of(colnames(all_groups)))

  # [2] corresponding parents
  sig_parents <- all_groups %>%
    dplyr::select(-groupRank) %>%
    dplyr::inner_join(sig_parents, by = c("groupId")) %>%
    dplyr::mutate(label = paste0(label, label_add)) %>%
    dplyr::select(dplyr::all_of(colnames(all_groups)))

  sig_parents_and_isotopes <- rbind(sig_parents, sig_isotopes)

  # [3] isotopes of corresponding parents not covered in significant isotopes
  all_other_groups <- all_groups %>%
    dplyr::filter(!groupId %in% sig_parents_and_isotopes$groupId) %>%
    dplyr::mutate(groupRank = 0)

  updated_peakgroups <- rbind(sig_parents_and_isotopes, all_other_groups) %>%
    dplyr::mutate(searchTableName = "clamDB") %>%
    dplyr::arrange(groupId)

  conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = mzrolldb_file_path)
  DBI::dbWriteTable(conn, "peakgroups", updated_peakgroups, overwrite = TRUE)
  DBI::dbExecute(conn, "VACUUM;")
  DBI::dbDisconnect(conn = conn)

  return(invisible(0))
}
