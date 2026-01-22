#' Run Uniform Manifold Approximation and Projection (UMAP) on ASV/feature data
#'
#' Perform UMAP dimensionality reduction on ASV/feature data for ASV clustering.
#' UMAP can be computed from PCA results (recommended) or directly from scaled data.
#' The result provides a 2D or 3D embedding of ASVs for visualization and clustering.
#'
#' @param object A \code{VarScope} object containing PCA results or scaled data.
#' @param reduction Character. Name of the reduction to use as input. If "pca",
#'   uses PCA scores (ASV coordinates in PC space). If NULL, uses the assay
#'   specified by \code{assay} parameter (default: "pca").
#' @param assay Character. Name of the assay to use if \code{reduction} is NULL
#'   (default: "scale").
#' @param dims Integer vector. Which dimensions (PCs) to use if reduction is "pca".
#'   If NULL, uses all available PCs (default: NULL).
#' @param n.neighbors Integer. Size of local neighborhood for UMAP (default: 30).
#' @param min.dist Numeric. Minimum distance between points in the embedding
#'   (default: 0.3).
#' @param n.components Integer. Number of UMAP dimensions to compute (default: 2).
#' @param metric Character. Distance metric to use (default: "cosine").
#' @param reduction.name Character. Name to store the UMAP reduction
#'   (default: "umap").
#' @param seed Integer. Random seed for reproducibility (default: 1115).
#'
#' @returns A \code{VarScope} object with UMAP results stored in
#'   \code{object$reductions[[reduction.name]]}. The reduction contains:
#'   \itemize{
#'     \item \code{embedding}: Matrix of UMAP coordinates (ASVs × UMAP dimensions)
#'       - ASV coordinates in UMAP space for clustering
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
#' vsobj <- RunUMAP(vsobj, dims = 1:10)
#' # Access ASV UMAP coordinates: vsobj$reductions$umap$embedding
#' # This is ASVs × UMAP dimensions for ASV clustering
#' }
RunUMAP <- function(object,
                    reduction = "pca",
                    assay = "scale",
                    dims = NULL,
                    n.neighbors = 30,
                    min.dist = 0.3,
                    n.components = 2,
                    metric = "cosine",
                    reduction.name = "umap",
                    seed = 1115) {

  set.seed(seed)

  # Check object
  stopifnot(inherits(object, "VarScope"))

  # Check if uwot package is available
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for UMAP. Please install it using: install.packages('uwot')")
  }

  # Get input data
  if (!is.null(reduction) && reduction %in% names(object$reductions)) {
    # Use reduction (e.g., PCA scores)
    reduction_obj <- object$reductions[[reduction]]
    
    if (reduction == "pca" && "scores" %in% names(reduction_obj)) {
      # Use PCA scores
      input_data <- reduction_obj$scores
      
      # Select dimensions if specified
      if (!is.null(dims)) {
        dim_names <- paste0("PC", dims)
        available_dims <- intersect(dim_names, colnames(input_data))
        if (length(available_dims) == 0) {
          stop("Specified dimensions not found in PCA reduction.")
        }
        input_data <- input_data[, available_dims, drop = FALSE]
        cat("Using", length(available_dims), "PCs for UMAP:",
            paste(available_dims, collapse = ", "), "\n")
      } else {
        # Use all available PCs
        cat("Using all", ncol(input_data), "PCs for UMAP.\n")
      }
    } else {
      stop("Reduction '", reduction, "' does not contain scores or is not supported.")
    }
  } else {
    # Use assay data directly (ASVs × samples)
    mat <- object$assays[[assay]]
    if (is.null(mat)) {
      stop("Assay '", assay, "' not found in object.")
    }
    
    # Use ASV data directly (no transpose needed)
    # mat is already ASVs × samples
    input_data <- mat
    
    # Optionally select features
    if (!is.null(object$features$stats_filter) &&
        nrow(object$features$stats_filter) > 0) {
      features <- object$features$stats_filter$feature
      features <- intersect(features, rownames(input_data))
      if (length(features) > 0) {
        input_data <- input_data[features, , drop = FALSE]
        cat("Using", length(features), "variable features for UMAP.\n")
      }
    }
    
    cat("Computing UMAP directly from assay '", assay, "' (ASV dimension).\n", sep = "")
  }

  # Check input data
  if (nrow(input_data) < n.neighbors) {
    warning("Number of ASVs (", nrow(input_data),
            ") is less than n.neighbors (", n.neighbors, "). ",
            "Setting n.neighbors to ", nrow(input_data) - 1, ".")
    n.neighbors <- max(1, nrow(input_data) - 1)
  }

  # Run UMAP
  cat("Running UMAP with n.neighbors =", n.neighbors,
      ", min.dist =", min.dist, ", n.components =", n.components, "...\n")
  
  umap_result <- uwot::umap(
    X = input_data,
    n_neighbors = n.neighbors,
    min_dist = min.dist,
    n_components = n.components,
    metric = metric,
    verbose = FALSE
  )

  # Set row and column names
  rownames(umap_result) <- rownames(input_data)
  colnames(umap_result) <- paste0("UMAP_", 1:n.components)

  # Create reduction object
  reduction_obj <- list(
    embedding = umap_result  # ASVs × UMAP dimensions
  )

  # Store in object
  object$reductions[[reduction.name]] <- reduction_obj

  # Store parameters
  object$params$RunUMAP <- list(
    reduction = reduction,
    assay = if (is.null(reduction)) assay else NULL,
    dims = dims,
    n.neighbors = n.neighbors,
    min.dist = min.dist,
    n.components = n.components,
    metric = metric,
    reduction.name = reduction.name
  )

  cat("UMAP completed. Embedding stored in reduction '", reduction.name, "'.\n", sep = "")

  return(object)
}
