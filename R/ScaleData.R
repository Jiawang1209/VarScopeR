#' Scale feature abundance data in VarScope object
#'
#' Scale (z-score standardize) feature abundance data by centering and scaling
#' each feature to have zero mean and unit variance. This is typically performed
#' after normalization and before PCA to ensure all features are on the same scale.
#'
#' @param object A \code{VarScope} object containing normalized data.
#' @param assay Character. Name of the input assay to be scaled (default: "norm").
#' @param new.assay Character. Name of the scaled data (default: "scale").
#' @param features Character vector or NULL. Features to scale. If NULL, uses
#'   all features or variable features if available.
#' @param do.center Logical. Whether to center the data (subtract mean)
#'   (default: TRUE).
#' @param do.scale Logical. Whether to scale the data (divide by standard
#'   deviation) (default: TRUE).
#'
#' @returns A \code{VarScope} object containing the scaled data stored in
#'   \code{object$assays[[new.assay]]}. Scaling parameters (means and standard
#'   deviations) are stored in \code{object$misc$ScaleData}.
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
#' }
ScaleData <- function(object,
                      assay = "norm",
                      new.assay = "scale",
                      features = NULL,
                      do.center = TRUE,
                      do.scale = TRUE) {

  # Check object
  stopifnot(inherits(object, "VarScope"))

  # Get data matrix (features × samples)
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
      cat("Scaling", length(features), "variable features.\n")
    } else {
      features <- rownames(mat)
      cat("Scaling all", length(features), "features.\n")
    }
  }

  # Filter features
  features <- intersect(features, rownames(mat))
  if (length(features) == 0) {
    stop("No valid features found.")
  }

  # Subset matrix to selected features
  mat_subset <- mat[features, , drop = FALSE]

  # Calculate means and standard deviations for each feature (row)
  feature_means <- rowMeans(mat_subset, na.rm = TRUE)
  feature_sds <- apply(mat_subset, 1, sd, na.rm = TRUE)

  # Handle zero variance features
  zero_var <- feature_sds < 1e-10
  if (any(zero_var)) {
    warning("Found ", sum(zero_var), " features with zero variance. ",
            "These will not be scaled.")
    feature_sds[zero_var] <- 1  # Avoid division by zero
  }

  # Scale the data
  mat_scaled <- mat_subset

  if (do.center) {
    mat_scaled <- sweep(mat_scaled, 1, feature_means, "-")
  }

  if (do.scale) {
    mat_scaled <- sweep(mat_scaled, 1, feature_sds, "/")
  }

  # Store scaled data
  object$assays[[new.assay]] <- mat_scaled

  # Store scaling parameters
  object$misc$ScaleData <- list(
    assay = assay,
    new.assay = new.assay,
    features = features,
    means = feature_means,
    sds = feature_sds,
    do.center = do.center,
    do.scale = do.scale
  )

  # Store parameters
  object$params$ScaleData <- list(
    assay = assay,
    new.assay = new.assay,
    n_features = length(features),
    do.center = do.center,
    do.scale = do.scale
  )

  cat("Scaled data stored in assay '", new.assay, "'.\n", sep = "")

  return(object)
}
