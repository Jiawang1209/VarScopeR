#' Find ASV Clusters using Graph-based Community Detection
#'
#' Identify ASV clusters/modules using graph-based community detection algorithms
#' (Louvain or Leiden) on the shared nearest neighbor (SNN) graph. This function
#' groups ASVs into coherent modules based on their abundance pattern similarity.
#'
#' @param object A \code{VarScope} object containing a graph (typically SNN graph
#'   from \code{FindNeighbors}).
#' @param graph.name Character. Name of the graph to use for clustering
#'   (default: NULL, will search for SNN graph).
#' @param resolution Numeric. Resolution parameter for clustering. Higher values
#'   lead to more clusters (default: 0.5).
#' @param algorithm Character. Clustering algorithm to use: "louvain" or "leiden"
#'   (default: "louvain").
#' @param cluster.name Character. Name to store the cluster assignments
#'   (default: NULL, will be auto-generated).
#' @param seed Integer. Random seed for reproducibility (default: 1115).
#'
#' @returns A \code{VarScope} object with cluster assignments stored in
#'   \code{object$clusters}. Cluster assignments are stored as a named vector
#'   or factor with sample/feature names.
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
#' vsobj <- FindClusters(vsobj, resolution = 0.5)
#' }
FindClusters <- function(object,
                        graph.name = NULL,
                        resolution = 0.5,
                        algorithm = c("louvain", "leiden"),
                        cluster.name = NULL,
                        seed = 1115) {

  set.seed(seed)

  # Check object
  stopifnot(inherits(object, "VarScope"))

  # Check algorithm
  algorithm <- match.arg(algorithm)

  # Find graph to use
  if (is.null(graph.name)) {
    # Try to find SNN graph
    snn_graphs <- grep("snn", names(object$graphs), ignore.case = TRUE, value = TRUE)
    if (length(snn_graphs) > 0) {
      graph.name <- snn_graphs[1]
      cat("Using graph '", graph.name, "' for clustering.\n", sep = "")
    } else {
      # Try to find any graph
      if (length(object$graphs) == 0) {
        stop("No graphs found in object. Please run FindNeighbors first.")
      }
      graph.name <- names(object$graphs)[1]
      cat("Using graph '", graph.name, "' for clustering.\n", sep = "")
    }
  }

  if (!graph.name %in% names(object$graphs)) {
    stop("Graph '", graph.name, "' not found in object.")
  }

  graph <- object$graphs[[graph.name]]

  # Check if igraph is available for advanced algorithms
  if (algorithm == "leiden" && !requireNamespace("igraph", quietly = TRUE)) {
    warning("Package 'igraph' is required for Leiden algorithm. ",
            "Falling back to Louvain algorithm.")
    algorithm <- "louvain"
  }

  # Convert adjacency matrix to igraph object if needed
  if (algorithm == "louvain" || algorithm == "leiden") {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      # Simple clustering without igraph
      cat("Using simple graph-based clustering (igraph not available).\n")
      clusters <- simple_graph_clustering(graph, resolution = resolution)
    } else {
      # Use igraph for clustering
      cat("Using", algorithm, "algorithm for clustering...\n")

      # Convert to igraph
      # Create edge list from adjacency matrix
      graph[graph > 0] <- 1  # Binarize
      diag(graph) <- 0  # Remove self-loops

      # Create igraph object
      g <- igraph::graph_from_adjacency_matrix(
        graph,
        mode = "undirected",
        weighted = NULL,
        diag = FALSE
      )

      if (algorithm == "louvain") {
        # Louvain clustering
        cluster_result <- igraph::cluster_louvain(g, resolution = resolution)
        clusters <- as.numeric(igraph::membership(cluster_result))
      } else if (algorithm == "leiden") {
        # Leiden clustering
        if (requireNamespace("leiden", quietly = TRUE)) {
          # Use leiden package if available
          clusters <- leiden::leiden(g, resolution_parameter = resolution)
          clusters <- as.numeric(as.factor(clusters))
        } else {
          # Fall back to Louvain
          warning("Package 'leiden' not available. Using Louvain algorithm.")
          cluster_result <- igraph::cluster_louvain(g, resolution = resolution)
          clusters <- as.numeric(igraph::membership(cluster_result))
        }
      }

      # Name clusters
      names(clusters) <- rownames(graph)
    }
  } else {
    stop("Unknown algorithm: ", algorithm)
  }

  # Store cluster assignments
  if (is.null(cluster.name)) {
    cluster.name <- paste0("cluster_", graph.name, "_res", resolution)
  }

  object$clusters[[cluster.name]] <- clusters

  # Store parameters
  object$params$FindClusters <- list(
    graph.name = graph.name,
    resolution = resolution,
    algorithm = algorithm,
    cluster.name = cluster.name,
    n_clusters = length(unique(clusters))
  )

  cat("Clustering completed. Found", length(unique(clusters)),
      "clusters using", algorithm, "algorithm.\n")
  cat("Cluster assignments stored as '", cluster.name, "'.\n", sep = "")

  return(object)
}

# Helper function for simple graph clustering (when igraph is not available)
simple_graph_clustering <- function(graph, resolution = 0.5) {
  # Simple connected components clustering
  # This is a basic implementation
  n <- nrow(graph)
  clusters <- rep(0, n)
  names(clusters) <- rownames(graph)

  visited <- rep(FALSE, n)
  cluster_id <- 1

  # Binarize graph
  graph_bin <- (graph > 0) * 1
  diag(graph_bin) <- 0

  # Find connected components
  for (i in 1:n) {
    if (!visited[i]) {
      # BFS to find all connected nodes
      queue <- i
      visited[i] <- TRUE
      clusters[i] <- cluster_id

      while (length(queue) > 0) {
        current <- queue[1]
        queue <- queue[-1]

        neighbors <- which(graph_bin[current, ] > 0)
        for (neighbor in neighbors) {
          if (!visited[neighbor]) {
            visited[neighbor] <- TRUE
            clusters[neighbor] <- cluster_id
            queue <- c(queue, neighbor)
          }
        }
      }
      cluster_id <- cluster_id + 1
    }
  }

  return(clusters)
}
