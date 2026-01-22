#' Find Nearest Neighbors and Build Graphs for ASV Clustering
#'
#' Construct K-nearest neighbor (KNN) and shared nearest neighbor (SNN) graphs
#' based on PCA results or other dimensional reductions. These graphs are used
#' for subsequent ASV clustering analysis. The graphs represent similarities
#' between ASVs based on their abundance patterns.
#'
#' @param object A \code{VarScope} object containing PCA or other reduction results.
#' @param reduction Character. Name of the reduction to use (default: "pca").
#' @param dims Integer vector. Which dimensions to use from the reduction
#'   (default: 1:10).
#' @param k.param Integer. Number of nearest neighbors to find (default: 20).
#' @param compute.SNN Logical. Whether to compute the shared nearest neighbor
#'   (SNN) graph (default: TRUE).
#' @param prune.SNN Numeric. Sets the cutoff for acceptable Jaccard index when
#'   computing the neighborhood overlap for the SNN construction (default: 1/15).
#' @param graph.name Character. Prefix for graph names (default: NULL, uses
#'   reduction name).
#' @param seed Integer. Random seed for reproducibility (default: 1115).
#'
#' @returns A \code{VarScope} object with graphs stored in \code{object$graphs}.
#'   The graphs represent ASV-ASV similarities and include:
#'   \itemize{
#'     \item \code{knn}: K-nearest neighbor graph (ASVs × ASVs)
#'     \item \code{snn}: Shared nearest neighbor graph (ASVs × ASVs, if compute.SNN = TRUE)
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
#' vsobj <- FindNeighbors(vsobj, dims = 1:10)
#' }
FindNeighbors <- function(object,
                         reduction = "pca",
                         dims = 1:10,
                         k.param = 20,
                         compute.SNN = TRUE,
                         prune.SNN = 1/15,
                         graph.name = NULL,
                         seed = 1115) {

  set.seed(seed)

  # Check object
  stopifnot(inherits(object, "VarScope"))

  # Check if reduction exists
  if (!reduction %in% names(object$reductions)) {
    stop("Reduction '", reduction, "' not found in object.")
  }

  reduction_obj <- object$reductions[[reduction]]

  # Get input data
  if ("scores" %in% names(reduction_obj)) {
    # Use PCA scores
    input_data <- reduction_obj$scores
    dim_names <- paste0("PC", dims)
    available_dims <- intersect(dim_names, colnames(input_data))
    if (length(available_dims) == 0) {
      stop("Specified dimensions not found in reduction.")
    }
    input_data <- input_data[, available_dims, drop = FALSE]
    cat("Using", length(available_dims), "dimensions from", reduction, "reduction.\n")
  } else if ("embedding" %in% names(reduction_obj)) {
    # Use UMAP or other embedding
    input_data <- reduction_obj$embedding
    cat("Using embedding from", reduction, "reduction.\n")
  } else {
    stop("Reduction '", reduction, "' does not contain scores or embedding.")
  }

  # Check k.param
  n_asvs <- nrow(input_data)
  if (k.param >= n_asvs) {
    warning("k.param (", k.param, ") is >= number of ASVs (", n_asvs, "). ",
            "Setting k.param to ", n_asvs - 1, ".")
    k.param <- max(1, n_asvs - 1)
  }

  # Compute distance matrix between ASVs
  cat("Computing distance matrix between ASVs...\n")
  dist_matrix <- as.matrix(dist(input_data, method = "euclidean"))

  # Find K nearest neighbors for each ASV
  cat("Finding", k.param, "nearest neighbors for each ASV...\n")
  knn_indices <- matrix(0, nrow = n_asvs, ncol = k.param)
  knn_distances <- matrix(0, nrow = n_asvs, ncol = k.param)

  for (i in 1:n_asvs) {
    # Get distances from ASV i to all others
    dists <- dist_matrix[i, ]
    dists[i] <- Inf  # Exclude self

    # Find k nearest neighbors
    knn_idx <- order(dists)[1:k.param]
    knn_indices[i, ] <- knn_idx
    knn_distances[i, ] <- dists[knn_idx]
  }

  # Build KNN graph (sparse adjacency matrix)
  # Create edge list
  edge_list <- data.frame(
    from = rep(1:n_asvs, each = k.param),
    to = as.vector(knn_indices),
    weight = as.vector(knn_distances)
  )

  # Create symmetric KNN graph
  # If ASV A is neighbor of ASV B, then ASV B is also neighbor of ASV A
  edge_list_symmetric <- rbind(
    edge_list,
    data.frame(
      from = edge_list$to,
      to = edge_list$from,
      weight = edge_list$weight
    )
  )

  # Remove duplicates and self-loops
  edge_list_symmetric <- edge_list_symmetric[
    edge_list_symmetric$from != edge_list_symmetric$to,
  ]
  edge_list_symmetric <- unique(edge_list_symmetric[, c("from", "to", "weight")])

  # Create adjacency matrix
  asv_names <- rownames(input_data)
  n_asvs <- length(asv_names)

  knn_graph <- matrix(0, nrow = n_asvs, ncol = n_asvs)
  rownames(knn_graph) <- asv_names
  colnames(knn_graph) <- asv_names

  for (i in 1:nrow(edge_list_symmetric)) {
    from_idx <- edge_list_symmetric$from[i]
    to_idx <- edge_list_symmetric$to[i]
    weight <- edge_list_symmetric$weight[i]
    knn_graph[from_idx, to_idx] <- weight
  }

  # Store KNN graph
  graph_name_knn <- if (is.null(graph.name)) {
    paste0(reduction, "_knn")
  } else {
    paste0(graph.name, "_knn")
  }
  object$graphs[[graph_name_knn]] <- knn_graph

  # Compute SNN graph if requested
  if (compute.SNN) {
    cat("Computing shared nearest neighbor (SNN) graph...\n")

    # SNN: For each pair of ASVs, count how many neighbors they share
    snn_graph <- matrix(0, nrow = n_asvs, ncol = n_asvs)
    rownames(snn_graph) <- asv_names
    colnames(snn_graph) <- asv_names

    for (i in 1:(n_asvs - 1)) {
      for (j in (i + 1):n_asvs) {
        # Get neighbors of ASV i and ASV j
        neighbors_i <- knn_indices[i, ]
        neighbors_j <- knn_indices[j, ]

        # Count shared neighbors
        shared <- sum(neighbors_i %in% neighbors_j)

        # Jaccard index: shared / (k + k - shared)
        jaccard <- shared / (2 * k.param - shared)

        # Prune based on threshold
        if (jaccard >= prune.SNN) {
          snn_graph[i, j] <- jaccard
          snn_graph[j, i] <- jaccard
        }
      }
    }

    # Store SNN graph
    graph_name_snn <- if (is.null(graph.name)) {
      paste0(reduction, "_snn")
    } else {
      paste0(graph.name, "_snn")
    }
    object$graphs[[graph_name_snn]] <- snn_graph
  }

  # Store parameters
  object$params$FindNeighbors <- list(
    reduction = reduction,
    dims = dims,
    k.param = k.param,
    compute.SNN = compute.SNN,
    prune.SNN = prune.SNN,
    graph.name = graph.name
  )

  cat("Graph construction completed.\n")
  if (compute.SNN) {
    cat("KNN graph stored as '", graph_name_knn, "'.\n", sep = "")
    cat("SNN graph stored as '", graph_name_snn, "'.\n", sep = "")
  } else {
    cat("KNN graph stored as '", graph_name_knn, "'.\n", sep = "")
  }

  return(object)
}
