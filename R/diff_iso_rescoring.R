#' Re-score Differential Abundance Isotopic Peak Groups.
#'
#' @description
#' Compare differences in the isotopic envelopes between two sample sets - a
#' representative unlabeled set, and a representative labeled set.
#' Override the 'groupRank' column in an mzrollDB file with the new score value.
#'
#' @param mzrolldb_file MzrollDB file path, where \code{groupRank} column of
#' \code{peakgroups} will be updated after peak groups diff iso re-scoring.
#' @param mzML_dir local directory containing mzML files.
#' @param rescoring_function function for re-scoring peak groups based on
#'    differential isotopic abundance.
#' @param unlabeled_samples sample names corresponding to samples in
#'    representative unlabeled sample set.
#' @param labeled_samples sample names corresponding to smaples in
#'    representative labeled sample set.
#' @param diff_iso_params named list containing parameters that are used in the
#'    construction of individual peak group isotopic matrices. Primarily used
#'    by \code{mzkitcpp::ISO_isotope_matrices()}.
#' @param verbose if \code{TRUE}, print additional messages to the console.
#'
#' @return a tibble comparing the original diff iso scores to the re-scored
#' values. The output tibble has columns \code{groupId}, \code{compoundId},
#' \code{groupRank} (the original score), \code{label}, \code{rescoredRank}, and
#' \code{diff}, which is \code{groupRank - rescoredRank}. This way, the original
#' score and the updated score can be compared directly.
#'
#' @export
diff_iso_rescore <- function(
  mzrolldb_file,
  mzML_dir,
  rescoring_function,
  unlabeled_samples,
  labeled_samples,
  diff_iso_params,
  verbose = TRUE
) {
  # [1] Import saved sample, peakgroup, and peaks data
  samples <- PDB_sample_list(mzrolldb_file)
  groups <- PDB_peakgroups(mzrolldb_file)
  peaks <- PDB_peaks(mzrolldb_file)

  # [2] Generate iso matrices
  iso_matrices <- mzkitcpp::ISO_isotope_matrices(
    mzML_dir,
    samples,
    peaks,
    groups,
    unlabeled_samples,
    labeled_samples,
    diff_iso_params,
    FALSE
  )

  if (verbose) {
    cat(paste0(
      "mzkitcpp::ISO_isotope_matrices() returned ",
      nrow(iso_matrices),
      " isotopic measurements.\n"
    ))
  }

  # [3] Reshape outputs to list of isotopic matrices
  iso_matrices_reshaped <- to_iso_matrices(iso_matrices)

  if (verbose) {
    cat(paste0(
      "to_iso_matrices() returned ",
      length(iso_matrices_reshaped),
      " diff iso comparable isotopic matrices.\n"
    ))
  }

  # [4] diff iso re-scoring
  scores <- purrr::map(
    iso_matrices_reshaped,
    rescoring_function,
    unlabeled_samples,
    labeled_samples
  )
  scores_tibble <- tibble::tibble(
    groupId = as.integer(names(scores)),
    groupRank = unlist(scores)
  ) %>%
    dplyr::filter(groupRank > 0)

  if (verbose) {
    cat(paste0("Found ", nrow(scores_tibble), " diff iso scores >0.\n"))
  }

  # [5] update groups and re-save mzrollDB file
  groups_updated <- groups %>%
    dplyr::select(-groupRank) %>%
    dplyr::left_join(scores_tibble, by = c("groupId")) %>%
    dplyr::mutate(groupRank = ifelse(is.na(groupRank), 0, groupRank)) %>%
    dplyr::arrange(groupRank)

  conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = mzrolldb_file)

  DBI::dbWriteTable(conn, "peakgroups", groups_updated, overwrite = TRUE)

  DBI::dbDisconnect(conn)

  if (verbose) {
    cat(paste0(
      "Successfuly completed rescoring and saved updated peak groups to mzrollDB file!\n"
    ))
  }

  # [6] compute comparison of the original to new scoring approach.
  original_groups <- groups %>%
    dplyr::filter(parentGroupId == 0) %>%
    dplyr::select(groupId, compoundId, groupRank, label)

  updated_groups <- groups_updated %>%
    dplyr::filter(parentGroupId == 0) %>%
    dplyr::select(groupId, compoundId, groupRank) %>%
    dplyr::rename(rescoredRank = groupRank)

  group_comparison <- original_groups %>%
    dplyr::inner_join(updated_groups, by = c("groupId", "compoundId")) %>%
    dplyr::mutate(diff = groupRank - rescoredRank) %>%
    dplyr::arrange(-groupRank)

  # [7] return comparison
  return(group_comparison)
}

#' Re-score and label Differential Abundance Isotopic Peak Groups.
#'
#' @description
#' convenience function to copy an mzrollDB, re-score the values (via \code{diff_iso_rescore()}),
#' identify compounds of interest, label them, and save the rescored and labeled
#' peak groups into a new mzrollDB named from the initial mzrollDB with the
#' \code{-rescored.mzrollDB} suffix.
#'
#' @param original_mzrolldb_file MzrollDB file path, where \code{groupRank} column of
#' \code{peakgroups} will be updated after peak groups diff iso re-scoring.
#' @param mzML_dir local directory containing mzML or mzXML files.
#' @param rescoring_function function for re-scoring peak groups based on
#'    differential isotopic abundance.
#' @param unlabeled_samples_pattern string pattern to identify samples in
#'    representative unlabeled sample set.
#' @param labeled_samples_pattern string pattern to identify samples in
#'    representative labeled sample set.
#' @param rescore_suffix name to be appended to the \code{original_mzrolldb_file}
#'    (before the \code{.mzrollDB} file extension) in the renamed output file.
#' @param is_correct_natural_abundance When determining isotope matrices, optionally
#'    correct for natural abundance.
#' @param ms2_score_threshold minimum value for MS2 score to label peak group.
#'   This should be adjusted based on the scoring type, but for a cosine score
#'   or similar default value of \code{0.7} is fine.
#' @param rank_thresh minimum value for diff iso scoring to label peak group.
#'   Default value is \code{1.30103}, which corresponds to \code{-log10(0.05)}.
#' @param verbose if \code{TRUE}, print additional messages to the console.
#'
#' @return \code{invisible(0)}, as side effect, creates a new \code{mzrollDB} file
#' with the rescored values, named the same as the \code{original_mzrolldb_file}
#' except \code{rescore_suffix} is added to the name.
#'
#' @export
diff_iso_rescore_and_label <- function(
  original_mzrolldb_file,
  mzML_dir,
  rescoring_function,
  unlabeled_samples_pattern,
  labeled_samples_pattern,
  rescore_suffix = "-rescored",
  is_correct_natural_abundance = FALSE,
  ms2_score_threshold = 0.7,
  rank_thresh = 1.30103, # = -log10(.05)
  verbose = TRUE
) {
  rescored_mzrolldb_file <- file.path(
    dirname(original_mzrolldb_file),
    paste0(
      gsub(".mzrollDB", "", basename(original_mzrolldb_file)),
      rescore_suffix,
      ".mzrollDB"
    )
  )

  # start by copying original file into rescored file path
  cmd <- glue::glue("cp {original_mzrolldb_file} {rescored_mzrolldb_file}")
  if (verbose) {
    cat(paste0(
      "Rescoring will be executed in copied mzrollDB file:\n`",
      cmd,
      "`\n"
    ))
  }
  system(cmd)

  diff_iso_params <- list()
  diff_iso_params[["diffIsoScoringFractionOfSampleTotal"]] <- TRUE
  diff_iso_params[[
    "diffIsoScoringCorrectNatAbundance"
  ]] <- is_correct_natural_abundance

  all_samples_files <- list.files(mzML_dir, pattern = "*.mzX?ML")

  unlabeled_samples <- all_samples_files[grepl(
    unlabeled_samples_pattern,
    all_samples_files
  )]
  labeled_samples <- all_samples_files[grepl(
    labeled_samples_pattern,
    all_samples_files
  )]

  group_comparison <- metisotopes::diff_iso_rescore(
    mzrolldb_file = rescored_mzrolldb_file,
    mzML_dir = mzML_dir,
    rescoring_function = metisotopes::diff_iso_m_plus_zero_fraction_WelchTTest,
    unlabeled_samples = unlabeled_samples,
    labeled_samples = labeled_samples,
    diff_iso_params = diff_iso_params,
    verbose = verbose
  )

  rescored_groups <- PDB_peakgroups(rescored_mzrolldb_file)

  sig_hits_above_thresh <- group_comparison %>%
    dplyr::inner_join(rescored_groups, by = c("groupId")) %>%
    dplyr::filter(
      ms2Score >= ms2_score_threshold & rescoredRank >= rank_thresh
    ) %>%
    dplyr::select(groupId)

  labeled_groups <- rescored_groups %>%
    dplyr::mutate(
      label = ifelse(groupId %in% sig_hits_above_thresh$groupId, "c", label)
    )

  conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = rescored_mzrolldb_file)

  DBI::dbWriteTable(conn, "peakgroups", labeled_groups, overwrite = TRUE)

  DBI::dbDisconnect(conn)

  if (verbose) {
    cat(paste0(
      "Successfully identified and saved ",
      length(sig_hits_above_thresh$groupId),
      " compounds of interest in rescored file.\n"
    ))
  }

  return(invisible(0))
}

#' Re-score and label for isotopic incorporation and covariates of interest
#'
#' @description
#' First, score peak groups based on the likelihood of isotopic incorporation.
#' This may work, for example, by comparing the fractional abundance of the [M+0] isotope
#' among unlabeled samples and labeled samples.
#' Then, among those peak groups that exhibit sufficient labeling incorporation,
#' Compare the differential abundance observed between different sample subsets
#' for each individual isotopic peak, for only the labeled samples.
#' The manner by which labeled samples may be divided into groups is specified
#' via the input parameter \code{experimental_design}.
#' In the initial 'isotopic incorporation' stage, isotope matrices are calculated
#' by using fractional abundance.
#' In the subsequent labeled stage, isotope matrices are calculated using the full
#' abundance.
#' A single argument, \code{is_correct_natural_abundance}, may be used to determine
#' if isotopes should correct natural abundance or not.
#'
#' @param mzrolldb_file mzrollDB file path, where \code{groupRank} column of
#' \code{peakgroups} will be updated after peak groups diff iso re-scoring.
#' This file should have been generated with isotopes.
#' @param mzML_dir local directory containing mzML or mzXML files.
#' @param unlabeled_samples_pattern string pattern to identify samples in
#'    representative unlabeled sample set.
#' @param labeled_samples_pattern string pattern to identify samples in
#'    representative labeled sample set.
#' @param experimental_design tibble describing various ways to subset the dataset.
#'    Each subset will be evaluated for significance (only among labeled samples)
#'    via application of the \code{condition_rescoring_function}.
#' @param incorporation_rescoring_function function for re-scoring peak groups based on
#'    differential isotopic abundance.
#' @param condition_rescoring_function function for identifying individual isotopic
#'   peak groups that are significance relative to one condition versus another.
#'   These will be compared only among the labeled samples for peak groups
#'   that exhibit isotopic incorporation.
#' @param is_correct_natural_abundance When determining isotope matrices, optionally
#'    correct for natural abundance.
#' @param incorporation_score_threshold minimum value for diff iso scoring to label peak group.
#'   Default value is \code{1.30103}, which corresponds to \code{-log10(0.05)}.
#' @param condition_score_threshold minimum value for a given peak group to be
#'   considered rescored.
#'   Default value is \code{1.30103}, which corresponds to \code{-log10(0.05)}.
#' @param rescore_suffix name to be appended to the \code{original_mzrolldb_file}
#'    (before the \code{.mzrollDB} file extension) in the renamed output file.
#' @param verbose if \code{TRUE}, print additional messages to the console.
#'
#' @return
#' (1) A long table listing out all peak groups where incorporation
#'    is thought to have been observed.
#'    This is the same set of peak groups that is used as a filter
#'    for labeled-only sample significance testing.
#' (2) A long table listing out all peak groups (including isotopic peak groups)
#'    where a significant difference by one condition is observed.
#'
#' . Note also that a side-effect is the generation of a new mzrollDB file,
#' . named according to the \code{rescore_suffix} parameter, which is generated
#' . in the same directory as the input \code{original_mzrolldb_file}.
#'
#' @export
diff_iso_conditions_rescore_and_label <- function(
  mzrolldb_file,
  mzML_dir,
  unlabeled_samples_pattern,
  labeled_samples_pattern,
  experimental_design,
  incorporation_rescoring_function = metisotopes::diff_iso_m_plus_zero_fraction_WelchTTest,
  condition_rescoring_function = metisotopes::diff_iso_all_isotopes_WelchTTest,
  is_correct_natural_abundance = FALSE,
  rescore_suffix = "-rescored",
  incorporation_score_threshold = 1.30103, # = -log10(.05)
  condition_score_threshold = 1.30103, # = -log10(.05)
  verbose = TRUE
) {
  # [1] Isotopic Incorporation
  samples <- PDB_sample_list(mzrolldb_file)
  groups <- PDB_peakgroups(mzrolldb_file)
  peaks <- PDB_peaks(mzrolldb_file)

  parent_to_group <- groups %>%
    dplyr::select(groupId, parentGroupId, tagString) %>%
    dplyr::filter(parentGroupId != 0) %>%
    dplyr::rename(isotope = tagString, childGroupId = groupId)

  parent_groups <- groups %>%
    dplyr::filter(parentGroupId == 0) %>%
    dplyr::select(groupId, compoundName, adductName)

  group_summaries <- peaks %>%
    dplyr::inner_join(parent_groups, by = c("groupId")) %>%
    dplyr::group_by(groupId) %>%
    dplyr::mutate(groupMz = mean(peakMz), groupRt = mean(rt)) %>%
    dplyr::ungroup() %>%
    dplyr::select(groupId, groupMz, groupRt, compoundName, adductName) %>%
    dplyr::distinct()

  all_samples_files <- list.files(mzML_dir, pattern = "*.mzML")
  unlabeled_samples <- all_samples_files[grepl(
    unlabeled_samples_pattern,
    all_samples_files
  )]
  labeled_samples <- all_samples_files[grepl(
    labeled_samples_pattern,
    all_samples_files
  )]

  incorporation_diff_iso_params <- list()
  incorporation_diff_iso_params[["diffIsoScoringFractionOfSampleTotal"]] <- TRUE
  incorporation_diff_iso_params[[
    "diffIsoScoringCorrectNatAbundance"
  ]] <- is_correct_natural_abundance

  iso_matrices <- mzkitcpp::ISO_isotope_matrices(
    mzML_dir,
    samples,
    peaks,
    groups,
    unlabeled_samples,
    labeled_samples,
    incorporation_diff_iso_params,
    FALSE
  )

  iso_matrices_reshaped <- to_iso_matrices(iso_matrices)

  incorporation_scores <- purrr::map(
    iso_matrices_reshaped,
    incorporation_rescoring_function,
    unlabeled_samples,
    labeled_samples
  )

  incorporation_scores_tibble <- tibble::tibble(
    groupId = as.integer(names(incorporation_scores)),
    groupRank = unlist(incorporation_scores)
  ) %>%
    dplyr::arrange(desc(groupRank)) %>%
    dplyr::mutate(
      incorporation_label = ifelse(
        groupRank >= incorporation_score_threshold,
        "c",
        ""
      )
    )

  incorporation_scores_filtered <- incorporation_scores_tibble %>%
    dplyr::filter(groupRank >= incorporation_score_threshold)

  incorporation_scores_w_group <- group_summaries %>%
    dplyr::inner_join(incorporation_scores_filtered, by = c("groupId")) %>%
    dplyr::arrange(desc(groupRank))

  # [2] Condition Rescoring/Evaluation

  sig_incorporation_groups <- groups %>%
    dplyr::filter(
      groupId %in%
        incorporation_scores_w_group$groupId |
        parentGroupId %in% incorporation_scores_w_group$groupId
    )

  sig_incorporation_peaks <- peaks %>%
    dplyr::filter(groupId %in% sig_incorporation_groups$groupId)

  sig_diff_iso_params <- list()
  sig_diff_iso_params[["diffIsoScoringFractionOfSampleTotal"]] <- FALSE
  sig_diff_iso_params[[
    "diffIsoScoringCorrectNatAbundance"
  ]] <- is_correct_natural_abundance

  sig_iso_matrices <- mzkitcpp::ISO_isotope_matrices(
    mzML_dir,
    samples,
    sig_incorporation_peaks,
    sig_incorporation_groups,
    unlabeled_samples,
    labeled_samples,
    sig_diff_iso_params,
    FALSE
  )

  sig_incorporation_iso_matrices_labeled <- sig_iso_matrices %>%
    dplyr::filter(sample %in% labeled_samples)

  sig_incorporation_iso_matrices_reshaped <- to_iso_matrices(
    sig_incorporation_iso_matrices_labeled
  )

  conditions_comparisons <- vector(
    mode = "list",
    length = nrow(experimental_design)
  )
  names(conditions_comparisons) <- experimental_design$name

  # Iterate through each binary factor covariate, and assess the dataset
  # using Cartesian product approach.
  # [6] Define new sets, based on every possible covariate
  for (i in 1:nrow(experimental_design)) {
    # grab conditions set for this round of t-tests
    condition_1 <- experimental_design[[i, "condition_1"]]
    condition_2 <- experimental_design[[i, "condition_2"]]
    label <- experimental_design[[i, "label"]]

    condition_name <- paste0(condition_1, "_vs_", condition_2)

    condition_1_samples <- labeled_samples[grepl(condition_1, labeled_samples)]
    condition_2_samples <- labeled_samples[grepl(condition_2, labeled_samples)]

    condition_i_scores <- purrr::map(
      sig_incorporation_iso_matrices_reshaped,
      condition_rescoring_function,
      condition_1_samples,
      condition_2_samples
    )

    condition_i_scores_w_names <- dplyr::bind_rows(
      condition_i_scores,
      .id = "groupId"
    ) %>%
      dplyr::filter(score >= condition_score_threshold) %>%
      dplyr::mutate(comparison = condition_name, condition_label = label)

    conditions_comparisons[[i]] <- condition_i_scores_w_names
  }

  # Flatten conditions comparisons to single table
  conditions_comparisons_flattened <- purrr::reduce(
    conditions_comparisons,
    rbind
  ) %>%
    dplyr::mutate(groupId = as.integer(groupId))

  # Add back group summaries
  conditions_w_group_summaries <- conditions_comparisons_flattened %>%
    dplyr::inner_join(group_summaries, by = c("groupId")) %>%
    dplyr::arrange(desc(score))

  # [3] Create Modified mzrollDB (with labels/rescored values)

  child_group_label_updates <- conditions_w_group_summaries %>%
    dplyr::inner_join(
      parent_to_group,
      by = c("groupId" = "parentGroupId", "isotope")
    ) %>%
    dplyr::select(groupId, childGroupId, condition_label)

  child_group_reshaped <- tibble::tibble(
    groupId = c(
      child_group_label_updates$groupId,
      child_group_label_updates$childGroupId
    ),
    condition_label = rep(child_group_label_updates$condition_label, 2)
  ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(groupId, condition_label) %>%
    dplyr::group_by(groupId) %>%
    dplyr::mutate(condition_label = paste0(condition_label, collapse = "")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

  combined_peakgroup_updates <- incorporation_scores_tibble %>%
    dplyr::select(groupId, groupRank, incorporation_label) %>%
    dplyr::full_join(child_group_reshaped, by = c("groupId")) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        !is.na(incorporation_label) & !is.na(condition_label) ~ paste0(
          incorporation_label,
          condition_label
        ),
        !is.na(incorporation_label) &
          is.na(condition_label) ~ incorporation_label,
        is.na(incorporation_label) & !is.na(condition_label) ~ condition_label,
        is.na(incorporation_label) & is.na(condition_label) ~ "",
        TRUE ~ ""
      )
    ) %>%
    dplyr::mutate(groupRank = ifelse(!is.na(groupRank), groupRank, 0)) %>%
    dplyr::select(groupId, groupRank, label)

  updated_peakgroups <- groups %>%
    dplyr::select(-groupRank, -label) %>%
    dplyr::left_join(combined_peakgroup_updates, by = c("groupId"))

  rescored_mzrolldb_file <- file.path(
    dirname(mzrolldb_file),
    paste0(
      gsub(".mzrollDB", "", basename(mzrolldb_file)),
      rescore_suffix,
      ".mzrollDB"
    )
  )

  cmd <- glue::glue("cp {mzrolldb_file} {rescored_mzrolldb_file}")
  system(cmd)
  if (verbose) {
    cat(paste0(
      "Successfully created rescored mzrolldb file: ",
      rescored_mzrolldb_file
    ))
  }

  conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = rescored_mzrolldb_file)
  DBI::dbWriteTable(conn, "peakgroups", updated_peakgroups, overwrite = TRUE)
  DBI::dbDisconnect(conn)

  # [4] Return output tables
  return(
    list(
      incorporation_results = incorporation_scores_w_group,
      condition_results = conditions_w_group_summaries
    )
  )
}
