#!/usr/bin/env Rscript

###
identify_non_redundant_covariates <- function(data, cov_string) {
  all_covs <- strsplit(gsub("\\s*\\+\\s*$", "", cov_string), "\\s*\\+\\s*")[[1]]
  
  # check if two columns have the same level structure
  have_same_levels <- function(col1, col2) {
    levels1 <- unique(data[[col1]])
    levels2 <- unique(data[[col2]])
    return(length(levels1) == length(levels2) && 
           all(sapply(levels1, function(l) any(data[data[[col1]] == l, col2] == data[data[[col1]] == l, col2][1]))))
  }
  
  # Identify groups of columns with the same level structure
  redundant_groups <- list()
  for (i in 1:(length(all_covs) - 1)) {
    for (j in (i + 1):length(all_covs)) {
      if (have_same_levels(all_covs[i], all_covs[j])) {
        group <- c(all_covs[i], all_covs[j])
        redundant_groups[[length(redundant_groups) + 1]] <- group
      }
    }
  }
  
  merged_groups <- list()
  for (group in redundant_groups) {
    added <- FALSE
    for (i in seq_along(merged_groups)) {
      if (any(group %in% merged_groups[[i]])) {
        merged_groups[[i]] <- unique(c(merged_groups[[i]], group))
        added <- TRUE
        break
      }
    }
    if (!added) {
      merged_groups[[length(merged_groups) + 1]] <- group
    }
  }
  
  # Choose one representative from each group (the first one)
  to_keep <- sapply(merged_groups, function(group) group[1])
  to_keep <- c(to_keep, setdiff(all_covs, unlist(merged_groups)))
  
  # Remove any columns with only one unique value
  to_keep <- to_keep[sapply(to_keep, function(col) length(unique(data[[col]])) > 1)]
  
  # Return the non-redundant covariates as a string
  return(paste(to_keep, collapse = " + "))
}

### Define Function
inference_sceptre_m <- function(mudata, ...) {
  # convert MuData object to sceptre object
  sceptre_object <- sceptreIGVF::convert_mudata_to_sceptre_object(mudata)

  # extract set of discovery pairs to test
  pairs_to_test <- MultiAssayExperiment::metadata(mudata)$pairs_to_test |>
    as.data.frame()
  discovery_pairs <- pairs_to_test |>
    dplyr::rename(
      grna_target = intended_target_name,
      response_id = gene_id
    )

  df = sceptre_object@covariate_data_frame

  # assemble arguments to set_analysis_parameters()
  args_list <- list(...)
  formula_object <- args_list$formula_object
  cov_string <- args_list$cov_string

   # check nested covariates
  new_cov_string <- identify_non_redundant_covariates(df, cov_string)
  
  if("discovery_pairs" %in% names(args_list)){
    warning("The `discovery_pairs` argument is ignored. The `discovery_pairs` are set from the `pairs_to_test` metadata.")
  }
  args_list[["discovery_pairs"]] <- discovery_pairs
  if (!"formula_object" %in% names(args_list) || formula_object == "default") {
  args_list$formula_object <- stats::formula(~ log(response_n_nonzero) + log(response_n_umis))
} else {
  if (new_cov_string != "") {
    args_list$formula_object <- stats::formula(sprintf("~ %s + log(response_n_nonzero) + log(response_n_umis)", new_cov_string))
  }
}
  args_list$cov_string <- NULL
  args_list$sceptre_object <- sceptre_object

  # set analysis parameters
  sceptre_object <- do.call(sceptre::set_analysis_parameters, args_list)

  # extract gRNA assignment and turn off QC
  sceptre_object <- sceptre_object |>
    sceptre::assign_grnas(method = "thresholding", threshold = 1) |>
    sceptre::run_qc(n_nonzero_trt_thresh = 0L,
                    n_nonzero_cntrl_thresh = 0L,
                    p_mito_threshold = 1)

  # run discovery analysis
  sceptre_object <- sceptre_object |>
    sceptre::run_discovery_analysis()

  # get results
  discovery_results <- sceptre_object |>
    sceptre::get_result(analysis = "run_discovery_analysis") |>
    dplyr::select(response_id, grna_target, p_value, log_2_fold_change) |>
    dplyr::rename(gene_id = response_id,
                  intended_target_name = grna_target,
                  log2_fc = log_2_fold_change)

  # add results to MuData
  test_results <- pairs_to_test |>
    dplyr::left_join(discovery_results, by = c("intended_target_name", "gene_id"))
  MultiAssayExperiment::metadata(mudata)$test_results <- test_results

  # return MuData
  return(mudata)
}

### Run Command
args <- commandArgs(trailingOnly = TRUE)
args <- readLines(commandArgs(trailingOnly = TRUE)[1])

# obtain the command line arguments
mudata_fp <- args[1]
side <- args[2]
grna_integration_strategy <- args[3]
resampling_approximation <- args[4]
control_group <- args[5]
resampling_mechanism <- args[6]
formula_object <- args[7]
cov_string <- args[8]

# process formula object
# if (!identical(formula_object, "default")) {
#  formula_object <- stats::formula(formula_object)
# }

# read MuData
mudata_in <- MuData::readH5MU(mudata_fp)

# run sceptre inference
mudata_out <- inference_sceptre_m(
  mudata = mudata_in,
  side = side,
  grna_integration_strategy = grna_integration_strategy,
  resampling_approximation = resampling_approximation,
  control_group = control_group,
  resampling_mechanism = resampling_mechanism,
  formula_object = formula_object,
  cov_string = cov_string
)

# write MuData
MuData::writeH5MU(object = mudata_out, file = 'sceptre_inference_mudata.h5mu')

