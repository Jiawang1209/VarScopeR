#' Run Principal Component Analysis (PCA) on ASV/feature data
#'
#' Perform PCA on ASV/feature abundance data. The analysis is performed on the
#' ASV/feature dimension, returning ASV coordinates in PC space (PC1-PC10 by default)
#' for ASV clustering. Each ASV is represented as a point in the PC space based on
#' its abundance pattern across samples.
#'
#' @param object A \code{VarScope} object containing normalized or scaled data.
#' @param assay Character. Name of the assay to use for PCA (default: "scale").
#' @param features Character vector or NULL. Features to use for PCA. If NULL,
#'   uses all features or variable features if available.
#' @param npcs Integer. Number of principal components to compute (default: 10).
#' @param reduction.name Character. Name to store the PCA reduction
#'   (default: "pca").
#' @param seed Integer. Random seed for reproducibility (default: 1115).
#'
#' @returns A \code{VarScope} object with PCA results stored in
#'   \code{object$reductions[[reduction.name]]}. The reduction contains:
#'   \itemize{
#'     \item \code{scores}: Matrix of ASV coordinates in PC space (ASVs × PCs)
#'       - This is the main output for ASV clustering
#'     \item \code{loadings}: Matrix of sample loadings (samples × PCs)
#'       - Contribution of each sample to each PC
#'     \item \code{sdev}: Standard deviations of principal components
#'     \item \code{var.explained}: Proportion of variance explained by each PC
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' vsobj <- CreateVarScopeObject(
#'   counts = "inst/extdata/otutab.txt",
#'   meta.data = "inst/extdata/group.csv"
#' )
#' vsobj <- NormalizeData(vsobj, normalization.method = "log2")
#' vsobj <- FindVariableFeatures(vsobj)
#' vsobj <- ScaleData(vsobj)
#' vsobj <- RunPCA(vsobj, npcs = 10)
#' 
#' # Access ASV coordinates in PC space (for ASV clustering)
#' asv_pc_coords <- vsobj$reductions$pca$scores  # ASVs × PCs
#' 
#' # Access sample loadings
#' sample_loadings <- vsobj$reductions$pca$loadings  # samples × PCs
#' }
RunPCA <- function(object,
                   assay = "scale",
                   features = NULL,
                   npcs = 10,
                   reduction.name = "pca",
                   seed = 1115) {

  set.seed(seed)

  # Check object
  stopifnot(inherits(object, "VarScope"))

  # Get data matrix (ASVs × samples)
  mat <- object$assays[[assay]]
  if (is.null(mat)) {
    stop("Assay '", assay, "' not found in object.")
  }

  # Select features
  if (is.null(features)) {
    # Use variable features if available, otherwise use all features
    if (!is.null(object$features$stats_filter) &&
        nrow(object$features$stats_filter) > 0) {
      features <- object$features$stats_filter$feature
      cat("Using", length(features), "variable features for PCA.\n")
    } else {
      features <- rownames(mat)
      cat("Using all", length(features), "features for PCA.\n")
    }
  }

  # Filter features
  features <- intersect(features, rownames(mat))
  if (length(features) == 0) {
    stop("No valid features found.")
  }
  if (length(features) < npcs) {
    warning("Number of features (", length(features),
            ") is less than npcs (", npcs, "). Setting npcs to ",
            length(features), ".")
    npcs <- length(features)
  }

  # Subset matrix to selected features
  mat_subset <- mat[features, , drop = FALSE]

  # Remove ASVs with zero variance across samples
  asv_vars <- apply(mat_subset, 1, var, na.rm = TRUE)
  valid_asvs <- asv_vars > 1e-10
  if (sum(valid_asvs) < npcs) {
    warning("Only ", sum(valid_asvs), " ASVs have non-zero variance.")
    npcs <- min(npcs, sum(valid_asvs))
  }
  mat_subset <- mat_subset[valid_asvs, , drop = FALSE]

  if (nrow(mat_subset) == 0) {
    stop("No ASVs with non-zero variance found.")
  }

  cat("Performing PCA on", nrow(mat_subset), "ASVs across", ncol(mat_subset), "samples...\n")

  # Perform PCA on ASV dimension
  # Input: ASVs × samples matrix
  # Each row is an ASV, each column is a sample
  # PCA will find principal components in the sample space
  # Result: ASV coordinates in PC space (for ASV clustering)
  
  # If data is already scaled, don't scale again
  if (assay == "scale" && !is.null(object$misc$ScaleData)) {
    # Data is already scaled per ASV (row-wise), so we don't need to scale again
    # But we still center (though it should already be centered)
    pca_result <- stats::prcomp(mat_subset, center = TRUE, scale. = FALSE)
  } else {
    # Data not scaled, center and scale
    # Scale each ASV (row) to have mean 0 and sd 1
    pca_result <- stats::prcomp(mat_subset, center = TRUE, scale. = TRUE)
  }

  # Extract components
  # Scores: ASVs × PCs (x matrix) - ASV coordinates in PC space (for clustering)
  asv_scores <- pca_result$x[, 1:min(npcs, ncol(pca_result$x)), drop = FALSE]
  rownames(asv_scores) <- rownames(mat_subset)
  colnames(asv_scores) <- paste0("PC", 1:ncol(asv_scores))

  # Loadings: samples × PCs (rotation matrix) - sample contributions to PCs
  sample_loadings <- pca_result$rotation[, 1:min(npcs, ncol(pca_result$rotation)), drop = FALSE]
  rownames(sample_loadings) <- colnames(mat_subset)
  colnames(sample_loadings) <- paste0("PC", 1:ncol(sample_loadings))

  # Standard deviations
  sdev <- pca_result$sdev[1:min(npcs, length(pca_result$sdev))]

  # Variance explained
  var_explained <- (sdev^2) / sum(pca_result$sdev^2)

  # Create reduction object
  reduction <- list(
    scores = asv_scores,        # ASVs × PCs - ASV coordinates for clustering
    loadings = sample_loadings,  # samples × PCs - sample contributions
    sdev = sdev,
    var.explained = var_explained,
    features.used = rownames(mat_subset)
  )

  # Store in object
  object$reductions[[reduction.name]] <- reduction

  # Store parameters
  object$params$RunPCA <- list(
    assay = assay,
    npcs = npcs,
    n_asvs_used = nrow(mat_subset),
    reduction.name = reduction.name
  )

  cat("PCA completed. Computed", ncol(asv_scores), "principal components.\n")
  cat("ASV coordinates stored in 'scores' (", nrow(asv_scores), " ASVs × ",
      ncol(asv_scores), " PCs).\n", sep = "")
  cat("Variance explained by PC1-PC", min(5, ncol(asv_scores)), ": ",
      paste(round(var_explained[1:min(5, length(var_explained))] * 100, 2),
            "%", collapse = ", "), "\n", sep = "")

  return(object)
}
