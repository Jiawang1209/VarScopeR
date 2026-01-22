#' Create a VarScope object
#'
#' Create a \code{VarScope} object from an ASV/feature table and corresponding sample metadata.
#' This function initializes the core container used throughout the VarScopeR workflow.
#'
#' @param counts Character, matrix, or data.frame. Path to a feature table file, or a matrix/data.frame
#'   with features as rows and samples as columns.
#' @param meta.data Character, or data.frame. Path to a metadata file, or a data.frame containing
#'   sample metadata. Must contain columns: \code{Sample} and \code{Group}.
#' @param project Character. Project name for the VarScope object (optional).
#' @param min.group Numeric. Minimum number of samples per group. Groups with fewer samples will
#'   be removed (default: 0, no filtering).
#' @param min.features Numeric. Minimum number of features (non-zero counts) per sample. Samples
#'   with fewer features will be removed (default: 0, no filtering).
#' @param ... Additional arguments passed to \code{read_counts_and_meta} when reading from files.
#'
#' @returns A \code{VarScope} object containing:
#'   \itemize{
#'     \item \code{assays$raw}: Raw feature count matrix
#'     \item \code{meta$sample}: Sample metadata
#'     \item \code{params}: Parameters used in object creation
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From file paths
#' vsobj <- CreateVarScopeObject(
#'   counts = "inst/extdata/otutab.txt",
#'   meta.data = "inst/extdata/group.csv",
#'   project = "ExampleProject"
#' )
#'
#' # From R objects
#' counts_matrix <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
#' rownames(counts_matrix) <- paste0("ASV_", 1:100)
#' colnames(counts_matrix) <- paste0("Sample_", 1:10)
#' 
#' meta_df <- data.frame(
#'   Sample = paste0("Sample_", 1:10),
#'   Group = rep(c("Control", "Treatment"), each = 5)
#' )
#' 
#' vsobj <- CreateVarScopeObject(
#'   counts = counts_matrix,
#'   meta.data = meta_df,
#'   project = "MyProject"
#' )
#' 
#' vsobj
#' }
CreateVarScopeObject <- function(counts,
                                 meta.data,
                                 project = NULL,
                                 min.group = 0,
                                 min.features = 0,
                                 ...){

  # Check if counts is a file path or data object
  if (is.character(counts) && length(counts) == 1 && file.exists(counts)) {
    # Read from file
    if (!is.character(meta.data) || length(meta.data) != 1 || !file.exists(meta.data)) {
      stop("If counts is a file path, meta.data must also be a valid file path.")
    }
    
    cat("Reading data from files...\n")
    tmp <- read_counts_and_meta(counts, meta.data, ...)
    counts <- tmp$counts
    group  <- tmp$meta
    
    # Check for mismatches
    if (length(tmp$check$missing_in_meta) > 0) {
      warning("Samples in counts table but not in metadata: ",
              paste(tmp$check$missing_in_meta, collapse = ", "))
    }
    if (length(tmp$check$missing_in_counts) > 0) {
      warning("Samples in metadata but not in counts table: ",
              paste(tmp$check$missing_in_counts, collapse = ", "))
    }
    
    # Filter to common samples
    common_samples <- intersect(colnames(counts), group$Sample)
    if (length(common_samples) == 0) {
      stop("No common samples found between counts table and metadata.")
    }
    counts <- counts[, common_samples, drop = FALSE]
    group <- group[group$Sample %in% common_samples, , drop = FALSE]
    
  } else {
    # Counts is a matrix or data.frame
    counts <- as.matrix(counts)
    storage.mode(counts) <- "double"
    
    # Validate counts matrix
    if (any(is.na(counts))) {
      warning("Counts matrix contains NA values. They will be converted to 0.")
      counts[is.na(counts)] <- 0
    }
    if (any(counts < 0)) {
      stop("Counts matrix contains negative values.")
    }
    
    # Validate metadata
    if (!is.data.frame(meta.data)) {
      stop("meta.data must be a data.frame when counts is a matrix/data.frame.")
    }
    if (!all(c("Sample", "Group") %in% colnames(meta.data))) {
      stop("meta.data must contain columns: 'Sample' and 'Group'.")
    }
    
    group <- meta.data
    
    # Check sample matching
    if (is.null(colnames(counts))) {
      stop("Counts matrix must have column names (sample IDs).")
    }
    if (is.null(rownames(counts))) {
      rownames(counts) <- paste0("Feature_", 1:nrow(counts))
      warning("Counts matrix has no row names. Auto-generated feature names.")
    }
    
    # Filter to common samples
    common_samples <- intersect(colnames(counts), group$Sample)
    if (length(common_samples) == 0) {
      stop("No common samples found between counts matrix and metadata.")
    }
    if (length(common_samples) < ncol(counts) || length(common_samples) < nrow(group)) {
      warning("Some samples are missing in either counts or metadata. ",
              "Using only common samples.")
    }
    counts <- counts[, common_samples, drop = FALSE]
    group <- group[group$Sample %in% common_samples, , drop = FALSE]
  }

  # Validate data dimensions
  if (nrow(counts) == 0) {
    stop("Counts matrix has no features (rows).")
  }
  if (ncol(counts) == 0) {
    stop("Counts matrix has no samples (columns).")
  }
  if (nrow(group) == 0) {
    stop("Metadata has no samples.")
  }

  # Filter samples by minimum features
  if (min.features > 0) {
    n_features_per_sample <- colSums(counts > 0)
    keep_samples <- n_features_per_sample >= min.features
    if (sum(keep_samples) < ncol(counts)) {
      removed <- sum(!keep_samples)
      cat("Removing", removed, "samples with <", min.features, "features.\n")
      counts <- counts[, keep_samples, drop = FALSE]
      group <- group[group$Sample %in% colnames(counts), , drop = FALSE]
    }
  }

  # Filter groups by minimum samples
  if (min.group > 0) {
    group_counts <- table(group$Group)
    keep_groups <- names(group_counts)[group_counts >= min.group]
    if (length(keep_groups) < length(group_counts)) {
      removed <- length(group_counts) - length(keep_groups)
      cat("Removing", removed, "groups with <", min.group, "samples.\n")
      group <- group[group$Group %in% keep_groups, , drop = FALSE]
      counts <- counts[, group$Sample, drop = FALSE]
    }
  }

  # Final validation
  if (ncol(counts) == 0) {
    stop("No samples remaining after filtering.")
  }
  if (nrow(counts) == 0) {
    stop("No features remaining after filtering.")
  }
  if (length(unique(group$Group)) < 2) {
    warning("Only one group remaining. Some analyses may require at least 2 groups.")
  }

  # Ensure sample order matches
  if (!identical(colnames(counts), group$Sample)) {
    group <- group[match(colnames(counts), group$Sample), , drop = FALSE]
  }

  # Create VarScope object
  object <- VarScopeObject()

  # Store data
  object$assays$raw <- counts
  object$meta$sample <- group

  # Store project name
  if (!is.null(project)) {
    object$misc$project <- project
  }

  # Store parameters and summary
  object$params$CreateVarScopeObject <- list(
    n_features = nrow(counts),
    n_samples = ncol(counts),
    n_groups = length(unique(group$Group)),
    groups = table(group$Group),
    min.group = min.group,
    min.features = min.features,
    project = project
  )

  cat("VarScope object created:\n")
  cat("  ", nrow(counts), "features across", ncol(counts), "samples\n")
  cat("  ", length(unique(group$Group)), "groups:", 
      paste(names(table(group$Group)), collapse = ", "), "\n")

  return(object)

}

VarScopeObject <- function(){
  structure(
    list(
      assays = list(), # matrix(raw/nrom/scale)
      meta = list(), # meta data
      features = list(), # features info
      reductions = list(), # PCA/UMAP/tSNE
      graphs = list(), # KNN / SNN
      clusters = list(), # Cluster / Modularity
      params = list(), # params in every step
      projcet = list(), # project info
      misc = list()
    ),
    class = "VarScope"
  )
}


