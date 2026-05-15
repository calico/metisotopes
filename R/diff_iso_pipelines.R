#' Pipeline Diff Iso Search
#'
#' @description Full Diff Iso pipeline - provide some configuration parameters,
#' search parameters, and iso parameters.
#' Note that this pipeline does not involve conditional variables - only focused
#' On determining which samples show evidence of isotopic incorporation
#' (between labeled and unlabeled samples)
#' This function will create two files in the designated \code{output_directory}
#' : (1) a basic mzrollDB results file, and (2) a rescored mzrollDB file.
#' If there is any file named \code{peakdetector.mzrollDB} in the \code{output_directory},
#' it will be deleted.
#'
#' @param peakdetector_executable absolute path to compiled peakdetector executable.
#' @param peakdetector_methods_folder absolute path to peakdetector methods folder.
#' @param peakdetector_params list of peakdetector parameters.
#' @param sample_directory absolute path to directory containing all \code{mzML} files.
#' @param output_directory output directory for mzrollDB files.
#' @param unscored_file_name name of mzrollDB output file.
#' @param rescore_suffix added to file name in output directory for re-scored version
#'   of output mzrollDB file.
#' @param is_correct_natural_abundance When determining isotope matrices, optionally
#'    correct for natural abundance.
#' @param unlabeled_samples_pattern string pattern to identify samples in
#'    representative unlabeled sample set.
#' @param labeled_samples_pattern string pattern to identify samples in
#'    representative labeled sample set.
#' @param rank_thresh minimum value for diff iso scoring to label peak group.
#'   Default value is \code{1.30103}, which corresponds to \code{-log10(0.05)}.
#' @param verbose if \code{TRUE}, print additional messages to the console.
#'
#' @export
pipeline_diff_iso_search <- function(
  peakdetector_executable,
  peakdetector_methods_folder,
  peakdetector_params,
  sample_directory,
  output_directory,
  unscored_file_name,
  rescore_suffix,
  is_correct_natural_abundance,
  unlabeled_samples_pattern,
  labeled_samples_pattern,
  rank_thresh = 1.30103, # = -log10(.05)
  verbose = TRUE
) {
  # [1] Prepare peakdetector command line
  cmd <- peakdetector_command_line(
    peakdetector_executable,
    peakdetector_methods_folder,
    sample_directory,
    output_directory,
    peakdetector_params
  )

  # [2] Remove any stale plain 'peakdetector.mzrollDB' files (these are overwritten anyway)
  default_mzrolldb_file <- file.path(output_directory, "peakdetector.mzrollDB")
  if (file.exists(default_mzrolldb_file)) {
    system(glue::glue("rm {default_mzrolldb_file}"))
  }

  # [3] Create unscored mzrollDB file
  tictoc::tic("Peakgroup detection, isotope extraction, and compound identification")
  system(cmd, ignore.stdout = !verbose, ignore.stderr = TRUE)
  tictoc::toc()

  # [4] rename output file to desired name
  unscored_mzrolldb_file <- file.path(output_directory, paste0(unscored_file_name, ".mzrollDB"))
  system(glue::glue("mv {default_mzrolldb_file} {unscored_mzrolldb_file}"))

  ## Start Re-scoring Part

  # [5] Manipulate pipeline to appear more like GUI pipeline
  peakgroups <- PDB_peakgroups(unscored_mzrolldb_file)
  peakgroups_updated <- peakgroups %>% dplyr::mutate(searchTableName = "clamDB")

  samples <- PDB_sample_list(unscored_mzrolldb_file)
  samples_updated <- diff_iso_color_samples(samples, unlabeled_samples_pattern, labeled_samples_pattern)

  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = unscored_mzrolldb_file)
  DBI::dbWriteTable(con, "peakgroups", peakgroups_updated, overwrite = TRUE)
  DBI::dbWriteTable(con, "samples", samples_updated, overwrite = TRUE)
  DBI::dbDisconnect(conn = con)

  # [6] Perform Diff Iso Scoring
  tictoc::tic("Diff Iso Scoring")
  diff_iso_rescore_and_label(
    original_mzrolldb_file = unscored_mzrolldb_file,
    mzML_dir = sample_directory,
    rescoring_function = metisotopes::diff_iso_m_plus_zero_fraction_WelchTTest,
    unlabeled_samples_pattern = unlabeled_samples_pattern,
    labeled_samples_pattern = labeled_samples_pattern,
    is_correct_natural_abundance = is_correct_natural_abundance,
    rescore_suffix = rescore_suffix,
    ms2_score_threshold = 0.0, # since we are searching unknowns
    rank_thresh = rank_thresh,
    verbose = verbose
  )
  tictoc::toc()

  # [7] Rename re-scored file
  rescored_file_name <- paste0(gsub(".mzrollDB", "", basename(unscored_mzrolldb_file)), rescore_suffix, ".mzrollDB")
  rescored_mzrolldb_file <- file.path(output_directory, rescored_file_name)

  # [8] Collect Results, return as table
  groups <- PDB_peakgroups(rescored_mzrolldb_file)
  peaks <- PDB_peaks(rescored_mzrolldb_file)

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

  groups_of_interest <- PDB_peakgroups(rescored_mzrolldb_file) %>%
    dplyr::filter(label == "c") %>%
    dplyr::select(groupId, ms2Score, groupRank) %>%
    dplyr::inner_join(group_summaries, by = c("groupId")) %>%
    dplyr::select(compoundName, adductName, groupMz, groupRt, ms2Score, groupRank)

  return(groups_of_interest)
}

#' Pipeline Diff Iso Search
#'
#' @description Full Conditional Diff Iso pipeline - provide some configuration parameters,
#' search parameters, and iso parameters.
#' This function will create two files in the designated \code{output_directory}
#' : (1) a basic mzrollDB results file, and (2) a rescored mzrollDB file.
#' If there is any file named \code{peakdetector.mzrollDB} in the \code{output_directory},
#' it will be deleted.
#' This pipeline will use the conditional re-score to investigate any features
#' that are predicted to have undergone isotopic incorporation, and apply this to the supplied
#' experimental design structure to find overrepresentations.
#'
#' @param peakdetector_executable absolute path to compiled peakdetector executable.
#' @param peakdetector_methods_folder absolute path to peakdetector methods folder.
#' @param peakdetector_params list of peakdetector parameters.
#' @param sample_directory absolute path to directory containing all \code{mzML} files.
#' @param output_directory output directory for mzrollDB files.
#' @param unscored_file_name name of mzrollDB output file.
#' @param rescore_suffix added to file name in output directory for re-scored version
#'   of output mzrollDB file.
#' @param is_correct_natural_abundance When determining isotope matrices, optionally
#'    correct for natural abundance.
#' @param unlabeled_samples_pattern string pattern to identify samples in
#'    representative unlabeled sample set.
#' @param labeled_samples_pattern string pattern to identify samples in
#'    representative labeled sample set.
#' @param experimental_design tibble describing various ways to subset the dataset.
#'    Each subset will be evaluated for significance (only among labeled samples)
#'    via application of the \code{condition_rescoring_function}.
#' @param incorporation_score_threshold minimum value for diff iso scoring to label peak group.
#'   Default value is \code{1.30103}, which corresponds to \code{-log10(0.05)}.
#' @param condition_score_threshold minimum value for a given peak group to be
#'   considered rescored.
#'   Default value is \code{1.30103}, which corresponds to \code{-log10(0.05)}.
#' @param verbose if \code{TRUE}, print additional messages to the console.
#'
#' @export
pipeline_diff_iso_conditions_search <- function(
  peakdetector_executable,
  peakdetector_methods_folder,
  peakdetector_params,
  sample_directory,
  output_directory,
  unscored_file_name,
  rescore_suffix,
  is_correct_natural_abundance,
  unlabeled_samples_pattern,
  labeled_samples_pattern,
  experimental_design,
  incorporation_score_threshold = 1.30103, # = -log10(.05)
  condition_score_threshold = 1.30103, # = -log10(.05)
  verbose = TRUE
) {
  # [1] Prepare peakdetector command line
  cmd <- peakdetector_command_line(
    peakdetector_executable,
    peakdetector_methods_folder,
    sample_directory,
    output_directory,
    peakdetector_params
  )

  # [2] Remove any stale plain 'peakdetector.mzrollDB' files (these are overwritten anyway)
  default_mzrolldb_file <- file.path(output_directory, "peakdetector.mzrollDB")
  if (file.exists(default_mzrolldb_file)) {
    system(glue::glue("rm {default_mzrolldb_file}"))
  }

  # [3] Create unscored mzrollDB file
  tictoc::tic("Peakgroup detection, isotope extraction, and compound identification")
  system(cmd, ignore.stdout = !verbose, ignore.stderr = TRUE)
  tictoc::toc()

  # [4] rename output file to desired name
  unscored_mzrolldb_file <- file.path(output_directory, paste0(unscored_file_name, ".mzrollDB"))
  system(glue::glue("mv {default_mzrolldb_file} {unscored_mzrolldb_file}"))

  ## Start re-scoring Part

  # [5] Manipulate pipeline to appear more like GUI pipeline
  peakgroups <- PDB_peakgroups(unscored_mzrolldb_file)
  peakgroups_updated <- peakgroups %>% dplyr::mutate(searchTableName = "clamDB")

  samples <- PDB_sample_list(unscored_mzrolldb_file)
  samples_updated <- diff_iso_color_samples(samples, unlabeled_samples_pattern, labeled_samples_pattern)

  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = unscored_mzrolldb_file)
  DBI::dbWriteTable(con, "peakgroups", peakgroups_updated, overwrite = TRUE)
  DBI::dbWriteTable(con, "samples", samples_updated, overwrite = TRUE)
  DBI::dbDisconnect(conn = con)

  conditions_rescore_results <- diff_iso_conditions_rescore_and_label(
    mzrolldb_file = unscored_mzrolldb_file,
    mzML_dir = sample_directory,
    unlabeled_samples_pattern = unlabeled_samples_pattern,
    labeled_samples_pattern = labeled_samples_pattern,
    experimental_design = experimental_design,
    is_correct_natural_abundance = is_correct_natural_abundance,
    rescore_suffix = rescore_suffix,
    incorporation_score_threshold = incorporation_score_threshold,
    condition_score_threshold = condition_score_threshold,
    verbose = verbose
  )

  return(conditions_rescore_results)
}

#' Pipeline Time-Emergent Differential Abundance
#'
#' @description Time-Emergent Differential Abundance pipeline - given files and parameters,
#' return an mzrollDB file contaning peaks, groups, and isotopic envelopes.
#' Envelopes where isotopic incorporation are indicated with blue circles.
#' Envelopes where time-emergent differential abundance are indicated with magenta circles.
#' This is all encoded in maven as 'c' and 'a' labeled peak groups.
#'
#' Please see constituent functions \code{compute_isotopic_incorporation()}, \code{compute_diff_scores()},
#' and \code{compute_time_emergent_diff_linear_model()} for more details.
#'
#' @param peakdetector_executable absolute path to compiled peakdetector executable.
#' @param peakdetector_methods_folder absolute path to peakdetector methods folder.
#' @param peakdetector_params list of peakdetector parameters.
#' @param sample_directory absolute path to directory containing all \code{mzML} files.
#' @param output_directory output directory for mzrollDB files.
#' @param output_file_name name of mzrollDB output file.
#' @param isotope_quant_measurement_type Specific quant type used for evaluation.
#' For example, \code{peakArea}, \code{peakAreaTop}, \code{smoothedPeakArea}, or \code{peakIntensity}.
#' @param t_early_control_samples sample names corresponding to \code{(Time=t_early, Treatment=control)}
#' @param t_late_control_samples sample names corresponding to \code{(Time=t_late, Treatment=control)}
#' @param t_early_treatment_samples sample names corresponding to \code{(Time=t_early, Treatment=treatment)}
#' @param t_late_treatment_samples sample names corresponding to \code{(Time=t_late, Treatment=treatment)}
#' @param sample_order vector of all sample names, in desired order. Will be used for various sample display metrics.
#' @param peakdetector_file default: \code{NULL}. If provided and not NULL, skip the peakdetector
#' generation step, and instead start from the provided peakdetector file.
#' This file will be copied and rescored.
#' @param incorporation_score_threshold \code{-log10(p-value)} score threshold from Welch t-test,
#' carried out with \code{t_early} vs \code{t_late} each in \code{treatment} and \code{condition} sample subsets.
#' If a peak group is significant from either subset, the peak group passes.
#' @param diff_score_threshold Used for flagging things according to \code{compute_diff_scores()},
#' which is approximately a ratio of p-values score.
#' @param p_val_interaction_threshold - used for \code{compute_time_emergent_diff_linear_model()}.
#' Applies to interaction term from linear model. See \code{diff_iso_emergent_significance()}.
#' @param  p_val_t_early_threshold - used for \code{compute_time_emergent_diff_linear_model()}.
#' Applies to \code{t_early} comparison of \code{control} vs \code{treatment}.
#' See \code{compute_time_emergent_diff_linear_model()} and \code{diff_iso_emergent_significance()}.
#' @param p_val_t_late_threshold - used for \code{compute_time_emergent_diff_linear_model()}.
#' Applies to \code{t_early} comparison of \code{control} vs \code{treatment}.
#' See \code{compute_time_emergent_diff_linear_model()} and \code{diff_iso_emergent_significance()}.
#' the \code{p_val_early_diff} should be ABOVE a threshold parameter (p-value not significant)
#' the \code{p_val_late_diff} should be BELOW a threshold parameter (p-value significant).
#' @param is_M0_normalize if \code{TRUE}, every isotope value is divided by the [M+0] value from its respective sample,
#' in \code{compute_time_emergent_diff_linear_model()}. Note that values are always Log-transformed
#' for that linear model, this additional step first divides isotope values by corresponding [M+0]
#' value in the envelope.
#' @param verbose if \code{TRUE}, print additional messages to the console.
#'
#' @export
pipeline_time_emergent_differential_abundance <- function(
  peakdetector_executable,
  peakdetector_methods_folder,
  sample_directory,
  output_directory,
  output_file_name,
  peakdetector_params,
  isotope_quant_measurement_type,
  t_early_control_samples,
  t_early_treatment_samples,
  t_late_control_samples,
  t_late_treatment_samples,
  sample_order,
  peakdetector_file = NULL,
  incorporation_score_threshold = 1.30103, # -log10(0.05),
  diff_score_threshold = 1.30103, # -log10(0.05),
  p_val_interaction_threshold = 0.10,
  p_val_t_early_threshold = 0.05,
  p_val_t_late_threshold = 0.05,
  is_filter_M0_normalize = TRUE,
  verbose = TRUE
) {
  # need full path for output file
  mzrolldb_file_path <- file.path(output_directory, output_file_name)

  # hook that allows for recomputation, or skip peakdetector step
  # (only execute rescoring)
  if (is.null(peakdetector_file)) {
    # [1] Prepare peakdetector command line
    cmd <- clamshell::peakdetector_command_line(
      peakdetector_executable,
      peakdetector_methods_folder,
      sample_directory,
      output_directory,
      peakdetector_params
    )

    # [3] Create unscored mzrollDB file
    tictoc::tic("Peakgroup detection, isotope extraction, and compound identification")
    system(cmd, ignore.stdout = !verbose, ignore.stderr = TRUE)
    tictoc::toc()

    # [4] rename output file to desired name
    default_mzrolldb_file <- file.path(output_directory, "peakdetector.mzrollDB")
    file.rename(default_mzrolldb_file, mzrolldb_file_path)
  } else {
    file.copy(peakdetector_file, mzrolldb_file_path)
  }

  # [5] isotopic incorporation
  isotopic_incorporation_scores <- compute_isotopic_incorporation(
    mzrolldb_file_path,
    isotope_quant_measurement_type,
    t_early_control_samples,
    t_early_treatment_samples,
    t_late_control_samples,
    t_late_treatment_samples,
    sample_order
  )

  sig_scores <- isotopic_incorporation_scores %>%
    dplyr::filter(is_isotopic_incorporation) %>%
    dplyr::select(-is_isotopic_incorporation)

  incorporation_group_ids <- as.character(unique(sig_scores$groupId))

  iso_matrices_conditional_df <- get_precomputed_iso_df(
    iso_mzrolldb_file = mzrolldb_file_path,
    isotope_quant_measurement_type = isotope_quant_measurement_type,
    is_fractional_abundance = FALSE,
    sample_order = sample_order
  )

  iso_matrices_conditional_df_list <- to_iso_matrices(iso_matrices_conditional_df)

  incorporation_subset <- iso_matrices_conditional_df_list[incorporation_group_ids]

  # [6] differential isotopic incorporation via diff scores
  diff_scores <- compute_diff_scores(
    incorporation_subset,
    t_early_control_samples,
    t_early_treatment_samples,
    t_late_control_samples,
    t_late_treatment_samples,
    sample_order,
    diff_score_threshold
  )

  # [7] Time-emergent differential abundance
  lm_scores <- compute_time_emergent_diff_linear_model(
    incorporation_subset,
    t_early_control_samples,
    t_early_treatment_samples,
    t_late_control_samples,
    t_late_treatment_samples,
    sample_order,
    p_val_interaction_threshold,
    p_val_t_early_threshold,
    p_val_t_late_threshold,
    is_filter_M0_normalize
  )

  # [8] update mzrolldb file, labeling peak groups and updating groupRank column
  # with new score values
  top_hits <- lm_scores %>%
    dplyr::filter(is_p_val_interaction & is_p_val_early_diff & is_p_val_late_diff) %>%
    dplyr::group_by(groupId, isotope) %>%
    dplyr::mutate(groupRank = max(emergentScore)) %>%
    dplyr::ungroup() %>%
    dplyr::select(groupId, isotope, groupRank) %>%
    dplyr::distinct()

  label_isotopes_by_top_hits(
    mzrolldb_file_path,
    top_hits,
    sig_scores
  )

  # [9] return scoring results as output
  scoring_results <- list(
    "isotopic_incorporation_scores" = isotopic_incorporation_scores,
    "diff_scores" = diff_scores,
    "lm_scores" = lm_scores
  )

  return(scoring_results)
}
