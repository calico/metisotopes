#' Add single command line parameter to args string
#'
#' @description
#' Add a single \code{(key, value)} parameter pair to \code{args} string.
#' Note that different kinds of keys require different handling in terms of their
#' formatting in the \code{args} string.
#'
#' @param args Existing formatting arguments ready for peakdetector executable.
#' @param key parameter key, e.g. name from \code{params} named list.
#' @param val parameter value, e.g. stored value from \code{params} named list.
#'
#' @returns updated \code{args} string, with data from single \code{(key, val)}
#' parameter formatted appropriately and appended to input \code{args}.
#'
#' @export
peakdetector_add_CL_argument <- function(args, key, val) {
  single_dash_keys_add_quotes <- c("-9")
  double_dash_keys_add_quotes <- c("--isotopeParameters")
  if (!is.null(key)) {
    if (key %in% single_dash_keys_add_quotes) {
      key_value_pair <- paste0(key, "'", val, "'")
    } else if (key %in% double_dash_keys_add_quotes) {
      key_value_pair <- paste0(key, " '", val, "'")
    } else if (stringr::str_starts(key, "--")) {
      key_value_pair <- paste0(key, " ", val)
    } else {
      key_value_pair <- paste0(key, val)
    }
    args <- c(args, key_value_pair)
  }
  return(args)
}

#' Add parameters to peakdetector command line
#'
#' @description
#' Given a list of parameters, where names are command line option keys,
#' Generate a formatted string consumable by peakdetector executable.
#' Note that double-dashed parameters must follow single-dashed parameters,
#' so the formatted string is constructed to meet this requirement.
#' The keys are not checked, so if the peakdetector executable does not recognize a
#' key, it may cause an error or be ignored.
#'
#' @param args Existing formatting arguments ready for peakdetector executable.
#' @param params Named list, containing command line options for peakdetector.
#'
#' @returns  updated \code{args} string, with data from \code{params} argument
#' formatted as a string and appended to input \code{args}.
#'
#' @export
peakdetector_add_params <- function(args, params) {
  param_names <- names(params)
  is_double_dash <- grepl("^--", param_names)

  # add regular arguments first
  for (i in 1:length(params)) {
    if (!is_double_dash[i]) {
      args <- peakdetector_add_CL_argument(args, param_names[i], params[[i]])
    }
  }

  # add double dash argument list
  for (i in 1:length(params)) {
    if (is_double_dash[i]) {
      args <- peakdetector_add_CL_argument(args, param_names[i], params[[i]])
    }
  }

  return(args)
}

#' Add samples to peakdetector command line
#'
#' @description
#' Given a folder containing one or more \code{.mzML} or \code{.mzXML} files,
#' generate a formatted string consumable by peakdetector executable.
#'
#' @param args Existing formatting arguments ready for peakdetector executable.
#' @param sample_directory Directory containing sample files.
#'
#' @returns updated \code{args} string, with formatted samples string appended to
#' input \code{args}.
#'
#' @export
peakdetector_add_samples <- function(args, sample_directory) {
  samples <- list.files(sample_directory, full.names = TRUE, pattern = "*.mzX?ML$")
  sample_str <- paste0(samples, collapse = " ")
  args <- c(args, sample_str)
  return(args)
}

#' Add RT alignment info to peakdetector command line
#'
#' @description
#' Given a folder containing one or more \code{.apts} or \code{.rt} files,
#' generate a formatted string consumable by peakdetector executable.
#'
#' @param args Existing formatting arguments ready for peakdetector executable.
#' @param sample_directory Directory containing RT alignment files
#'
#' @returns updated \code{args} string, with formatted samples string appended to
#' input \code{args}.
#'
#' @export
peakdetector_add_rt_file <- function(args, sample_directory) {
  rt_files <- list.files(sample_directory, full.names = TRUE, pattern = "*.rt$|*.apts$")
  rt_str <- paste0(rt_files, collapse = " ")
  args <- c(args, rt_str)
  return(args)
}

#' Create peakdetector command line
#'
#' @description
#' Creates a command line with appropriate arguments and values that can be properly
#' executed by \code{peakdetector} executable.
#' Arguments are re-ordered, with all double-dashed parameters following all
#' single-dashed parameters.
#'
#' @param peakdetector_executable absolute path to peakdetector executable.
#' @param peakdetector_methods_folder folder containing \code{.model} files, for
#' use with peak quality classification. By default, this should be
#'  \code{<peakdetector_executable>/Contents/Resources/methods/default.model}
#' @param sample_directory Directory containing set of \code{.mzML} or \code{.mzXML} files.
#' @param output_directory Folder where \code{peakdetector.mzrollDB} output file will be written.
#' @param params list of parameters consumed by peakdetector executable. See \code{peakdetector_default_parameters()}.
#' @param is_save_ms2_scans [FALSE] if \code{TRUE}, peak group ms2 scans are saved. Otherwise, they are not.
#' @param is_ms1_search [FALSE] if \code{TRUE}, do not include requirement to perform ms2-based peak detection.
#'   When this option is \code{TRUE}, this will use the same parameters as are used in a MAVEN non-MS2 search.
#'
#' @returns valid command line ready for \code{system(cmd)}
#'
#' @export
peakdetector_command_line <- function(
  peakdetector_executable,
  peakdetector_methods_folder,
  sample_directory,
  output_directory,
  params,
  is_save_ms2_scans = FALSE,
  is_ms1_search = FALSE
) {
  # Initialize to empty vector to avoid NA problems
  args <- character(0)

  # MS1 Search Options
  if (is_ms1_search) {
    params[["-e"]] <- "B"
  } else {
    args <- "-2"
  }

  # Option to save MS2 scans into mzrollDB file
  args <- c(args, ifelse(is_save_ms2_scans, "-l1", "-l0"))

  # output and method directory always added
  args <- c(
    args,
    paste0("-o", output_directory), # output will be written to this folder
    paste0("-m", peakdetector_methods_folder) # folder containing 'default.model' model file for peak quality scoring
  )

  # Add params arguments
  args <- peakdetector_add_params(args, params)

  # Optional RT alignment file(s)
  args <- peakdetector_add_rt_file(args, sample_directory)

  # Always add samples last
  args <- peakdetector_add_samples(args, sample_directory)

  # Create fully formatted string
  cmd <- paste0(paste0(peakdetector_executable, " "), paste0(args, collapse = " "))

  return(cmd)
}

#' Default Peakdetector Parameters
#'
#' @description
#' Return a named list (e.g, \code{params}), based on an input \code{spectral_library_file}.
#' This function also serves as documentation, containing the default value of parameters from
#' the MAVEN GUI.
#'
#' @param spectral_library_file Full absolute path of spectral library, to be searched by peakdetector.
#'
#' @returns named list \code{params} with parameters, ready for consumption of downstream functions,
#' e.g. \code{peakdetector_command_line()}.
#'
#' @export
peakdetector_default_parameters <- function(spectral_library_file) {
  params <- list()

  # Parameters always set by GUI
  params[["-f"]] <- "E" # Peak Grouping Algorithm Type - should always be 'E'
  params[["--mergedSmoothedMaxToBoundsIntensityPolicy"]] <- "MINIMUM" # GUI always sets this value
  params[["--mergedPeakRtBoundsMaxIntensityFraction"]] <- -1 # GUI Always sets this value

  # Parameters that can be set by GUI, but not by CL
  # [Peak Detection -> Peak Picking And Grouping -> EIC Smoother Type] must be 'Gaussian' (peakdetector always sets to value 'Gaussian')
  #    (UI: cmbSmootherType)
  # [Peak Scoring -> Matching Options -> Require compound's associated adduct and searched adduct match for compound matches]
  #    must be checked (GUI will enumerate extra compounds based on selected adducts, peakdetector uses adducts only from library)
  #    (UI: chkRequireAdductMatch) [matchingIsRequireAdductPrecursorMatch] [matchingIsRequireAdductPrecursorMatch]
  #
  # Configurable dialog parameters

  # Peak Detection Tab

  # Feature Detection
  params[["-p"]] <- 20 # Mass Slice m/z Merge Tolerance
  params[["-r"]] <- 100 # Mass Slice RT Merge Tolerance

  # Peak Picking And Grouping
  params[["--peakRtBoundsSlopeThreshold"]] <- 0.01 # Peak Boundary Slope Threshold
  params[["--mergedPeakRtBoundsSlopeThreshold"]] <- 0.01 # Peak Boundary Slope Threshold
  params[["--peakRtBoundsMaxIntensityFraction"]] <- 0 # Peak Boundary Intensity Frac Threshold
  params[["-y"]] <- 5 # EIC Smoothing Window
  params[["-g"]] <- 0.25 # Peak Group Max RT Difference
  params[["-u"]] <- 0.80 # Peak Group Merge Overlap
  params[["--mergedSmoothedMaxToBoundsMinRatio"]] <- 1 # Peak Group S/N Threshold

  # Peak Scoring Tab

  # Peak Scoring
  params[["-z"]] <- 1 # Min Signal/Noise Ratio
  params[["-i"]] <- 1e4 # Min Highest Peak Intensity
  params[["-b"]] <- 1 # Min Good Peak/Group
  params[["--minSignalBlankRatio"]] <- 2 # Min Signal/Blank Ratio
  params[["-w"]] <- 5 # Min Peak Width
  params[["-q"]] <- 0 # Min Peak Quality

  # Baseline Computation
  params[["-8"]] <- 80 # Drop top x% intensities from chromatogram
  params[["-7"]] <- 5 # Baseline Smoothing Window
  params[["--eicBaselineEstimationType"]] <- "EIC_NON_PEAK_MEDIAN_SMOOTHED_INTENSITY" # Baseline Computation Type

  # Fragmentation Matching
  params[["-0"]] <- "metaboliteSearch" # Scoring Algorithm
  params[["-1"]] <- spectral_library_file # Metabolite Library

  # ms2MinNumMatches: Peak Scoring -> MS/MS Matching Scoring Settings -> Spectral Matches
  # rtMatchTolerance: Peak Detection Tab -> Compound Database -> Compound RT Match tolerance
  # ms1PpmTolr: Peak Detection Tab -> Compound Database -> Compound m/z Match tolerance
  params[["-9"]] <- "ms2MinNumMatches=0;rtMatchTolerance=2;ms1PpmTolr=20;"

  return(params)
}

#' Default Peakdetector Isotope Parameters
#'
#' @description
#' Return a named list (e.g, \code{isotope_params_list}),containing default value of parameters from
#' the MAVEN GUI for an isotopes search. Note that the returned list does not actually specify any
#' isotopes, so this will not work without modification.
#'
#' @returns named list \code{isotope_params_list} with parameters, which should be modified
#' with the correct isotopes (and other parameters), and then used by \code{mzkitcpp::mzk_get_isotope_parameters()}.
#'
#' @export
peakdetector_default_isotope_parameters <- function() {
  isotope_params_list <- list()
  isotope_params_list[["ppm"]] <- 10
  isotope_params_list[["labeledIsotopeRetentionPolicy"]] <- "ONLY_CARBON_TWO_LABELS"
  isotope_params_list[["isotopicExtractionAlgorithm"]] <- "MEIC_FWHM_RT_BOUNDS_AREA"
  isotope_params_list[["isCombineOverlappingIsotopes"]] <- TRUE
  isotope_params_list[["isExtractNIsotopes"]] <- TRUE
  isotope_params_list[["maxIsotopesToExtract"]] <- 7
  isotope_params_list[["maxIsotopeScanDiff"]] <- 5
  isotope_params_list[["minIsotopicCorrelation"]] <- 0.6
  isotope_params_list[["natAbundanceThreshold"]] <- 0.01
  isotope_params_list[["eic_smoothingWindow"]] <- 5
  isotope_params_list[["isIgnoreNaturalAbundance"]] <- FALSE
  isotope_params_list[["isKeepEmptyIsotopes"]] <- FALSE

  # peakPickingAndGroupingParameters
  isotope_params_list[["peakIsReassignPosToUnsmoothedMax"]] <- FALSE
  isotope_params_list[["peakRtBoundsMaxIntensityFraction"]] <- 0.00
  isotope_params_list[["peakRtBoundsSlopeThreshold"]] <- 0.01
  isotope_params_list[["peakBaselineDropTopX"]] <- 80
  isotope_params_list[["mergedBaselineDropTopX"]] <- 80
  isotope_params_list[["mergedIsComputeBounds"]] <- TRUE
  isotope_params_list[["mergedPeakRtBoundsSlopeThreshold"]] <- 0.01
  isotope_params_list[["mergedSmoothedMaxToBoundsMinRatio"]] <- 1
  isotope_params_list[["mergedSmoothedMaxToBoundsIntensityPolicy"]] <- "MINIMUM"

  return(isotope_params_list)
}

#' Metabolomics Parameters
#'
#' @description
#' Convenience Function to generate metabolomics parameters
#'
#' @returns encoded mzkitchen parameters string
#'
#' @export
peakdetector_metabolite_search_params <- function(ms2MinNumMatches = 0, rtIsRequireRtMatch = 0, rtMatchTolerance = 1) {
  mzkitchen_search_params <- paste0(
    "ms2MinNumMatches=", ms2MinNumMatches, ";",
    "ms2PpmTolr=20;",
    "ms1PpmTolr=10;",
    "rtIsRequireRtMatch=", rtIsRequireRtMatch, ";",
    "rtMatchTolerance=", rtMatchTolerance, ";",
    "scanFilterMinFracIntensity=0.01;",
    "consensusMs2MzRemovedStr='202.077,203.084';",
    "consensusMs2MzRemovedTol=10;",
    "consensusMinFractionMs2Scans=0.33;",
    "grpMs2MaxScanRtTolFromApex=1;",
    "grpMs2PurityTopN=4;",
    "grpMs2PurityThresholdAfterTopN=0.95;",

    # silent parameters (no GUI option)
    "grpMs2PurityTopNCode='w';",
    "grpMs2LabelAvgPurityCode='z';"
  )
  return(mzkitchen_search_params)
}
