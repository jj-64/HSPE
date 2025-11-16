#' Read LIS Microdata from .dta Files
#'
#' This function reads multiple LIS microdata files (Stata .dta format) from the
#' HSPE data directory, extracts a chosen variable and household weights, and
#' optionally expands the data using household weights.
#'
#' @param indices Integer vector of file indices to read.
#'   If `NA` (default), all files are read unless `country_filter` is used.
#'
#' @param variable Character string.
#'   The variable name to extract from each dataset (default = `"dhi"` for income).
#'
#' @param raw Logical.
#'   If `TRUE`, household weights are used to expand each record (replicate rows).
#'
#' @param country_filter Character or NULL.
#'   A text pattern or vector of patterns used to select only files whose
#'   filenames contain these substrings (e.g., `"at97"`).
#'   If multiple values: any matching file will be included.
#'
#' @return A named list where each element is a cleaned dataset for one file.
#' @export
#'
#' @examples
#' # Read all files:
#' # read_data()
#'
#' # Read only files containing "at97"
#' # read_data(country_filter = "at97")
#'
#' # Read 3 specific files with row replication
#' # read_data(indices = c(1, 3, 5), raw = TRUE)
#'
read_data <- function(indices = NA,
                          variable = "dhi",
                          raw = FALSE,
                          country_filter = NULL) {

  # ---- 1. List all .dta files ----
  files <- list.files(
    "C:/Users/User/Documents/HSPE/DATA/",
    pattern = "\\.dta$",
    full.names = TRUE
  )

  # Short file names (for labeling output)
  file_names <- basename(files)

  # ---- 2. Apply country filter (if provided) ----
  if (!is.null(country_filter)) {
    match_pattern <- paste(country_filter, collapse = "|")
    keep <- grepl(match_pattern, file_names, ignore.case = TRUE)

    files <- files[keep]
    file_names <- file_names[keep]

    if (length(files) == 0) {
      stop("No .dta files matched the provided country_filter.")
    }

    # If user supplied indices, apply to filtered files
    if (!is.na(indices[1])) {
      files <- files[indices]
      file_names <- file_names[indices]
    }

  } else {
    # No filter → use indices on full list
    if (is.na(indices[1])) indices <- seq_along(files)
    files <- files[indices]
    file_names <- file_names[indices]
  }

  # ---- 3. Prepare output list ----
  results <- vector("list", length(files))
  names(results) <- file_names

  # ---- 4. Process each file ----
  for (j in seq_along(files)) {

    file <- files[j]

    message("Reading file: ", file_names[j], " ...")

    data <- haven::read_dta(file)

    # Check variables exist
    if (!all(c("hwgt", variable) %in% names(data))) {
      warning("Skipping file ", file_names[j], ": required variables missing.")
      next
    }

    # Extract and clean
    microdata <- as.data.frame(data[, c("hwgt", variable)])
    microdata <- microdata[complete.cases(microdata), ]

    # ---- 5. Expand using household weights ----
    if (raw) {
      microdata$hwgt <- round(microdata$hwgt * 10, 0)
      microdata <- replicate_rows(microdata)
    }

    results[[j]] <- microdata
  }

  return(results)
}
