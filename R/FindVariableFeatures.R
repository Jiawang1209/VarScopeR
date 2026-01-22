#' Identify variable features by variance decomposition in VarScope object
#'
#' Identify perturbation-driven variable features by decomposing feature-wise
#' variance into within-group and between-group components. Features showing
#' low within-group variability and high between-group variability are retained
#' as variable features.
#'
#' @param object A \code{VarScope} object containing normalized feature abundance data.
#' @param assay Character. Name of the assay to use for variable feature detection
#'   (default: "norm").
#' @param var_within_cutoff Numeric. Quantile threshold for within-group variance.
#'   Features with within-group variance below this quantile are retained
#'   (default: 0.4, i.e., bottom 40%).
#' @param var_between_cutoff Numeric. Lower threshold for between-group variance
#'   proportion (eta2_between). Features with eta2_between >= this value are
#'   retained (default: 0.5).
#' @param seed Integer. Random seed for reproducibility (default: 1115).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @returns A \code{VarScope} object with variable feature information stored in:
#'   \itemize{
#'     \item \code{object$features$stats}: Statistics for all features
#'     \item \code{object$features$stats_filter}: Filtered variable features
#'     \item \code{object$features$variable}: Character vector of variable feature names
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
#' vsobj <- FindVariableFeatures(
#'   vsobj,
#'   var_within_cutoff = 0.4,
#'   var_between_cutoff = 0.5
#' )
#' # Access variable features
#' variable_features <- vsobj$features$variable
#' # View statistics
#' head(vsobj$features$stats_filter)
#' }
FindVariableFeatures <- function(object,
                                 assay = "norm",
                                 var_within_cutoff = 0.4,
                                 var_between_cutoff = 0.5,
                                 seed = 1115,
                                 verbose = TRUE
                                 ){

  set.seed(seed)

  # Check object
  if (!inherits(object, "VarScope")) {
    stop("object must be a VarScope object.")
  }

  # Get data
  mat <- object$assays[[assay]]
  if (is.null(mat)) {
    stop("Assay '", assay, "' not found in object.")
  }

  # Get metadata
  group <- object$meta$sample
  if (is.null(group)) {
    stop("Sample metadata not found in object.")
  }

  # Validate dimensions
  if (ncol(mat) != nrow(group)) {
    stop("Number of columns in matrix (", ncol(mat), ") does not match ",
         "number of rows in metadata (", nrow(group), ").")
  }

  # Match samples
  if (is.null(colnames(mat))) {
    stop("Counts matrix must have column names (sample IDs).")
  }
  if (is.null(rownames(mat))) {
    stop("Counts matrix must have row names (feature IDs).")
  }

  idx <- match(colnames(mat), group$Sample)
  if (any(is.na(idx))) {
    missing <- colnames(mat)[is.na(idx)]
    stop("Samples in matrix not found in metadata: ",
         paste(missing, collapse = ", "))
  }

  # Check for groups
  grp <- factor(group$Group[idx])
  groups <- levels(grp)
  
  if (length(groups) < 2) {
    stop("Need at least 2 groups for variance decomposition. Found: ",
         length(groups), " group(s).")
  }

  if (verbose) {
    cat("Identifying variable features...\n")
    cat("  Using", nrow(mat), "features across", ncol(mat), "samples\n")
    cat("  Groups:", paste(groups, collapse = ", "), "\n")
  }

  # Handle missing values
  if (any(is.na(mat))) {
    warning("Matrix contains NA values. Replacing with 0.")
    mat[is.na(mat)] <- 0
  }

  # Calculate grand mean (per feature)
  grand_mean <- rowMeans(mat, na.rm = TRUE)

  # Calculate group means (feature × group)
  group_means <- do.call(
    cbind,
    lapply(groups, function(g) {
      group_data <- mat[, grp == g, drop = FALSE]
      rowMeans(group_data, na.rm = TRUE)
    })
  )
  colnames(group_means) <- groups

  # Calculate weights based on group sizes
  group_sizes <- table(grp)[groups]
  weights <- as.numeric(group_sizes) / sum(group_sizes)

  # Calculate between-group variance
  deviation2 <- sweep(group_means, 1, grand_mean, "-")^2
  var_between <- rowSums(deviation2 * rep(weights, each = nrow(mat)))

  # Calculate within-group variance
  within_by_group <- do.call(
    cbind,
    lapply(groups, function(g) {
      group_data <- mat[, grp == g, drop = FALSE]
      apply(group_data, 1, function(x) {
        if (sum(!is.na(x)) < 2) return(0)
        var(x, na.rm = TRUE)
      })
    })
  )
  colnames(within_by_group) <- groups

  var_within <- rowSums(within_by_group * rep(weights, each = nrow(mat)))

  # Calculate total variance and eta2
  var_total <- var_within + var_between
  eta2_between <- ifelse(var_total > 0, var_between / var_total, NA)

  # Create output data frame
  out <- data.frame(
    feature = rownames(mat),
    var_within = var_within,
    var_between = var_between,
    var_total = var_total,
    eta2_between = eta2_between,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Filter variable features
  var_within_threshold <- quantile(out$var_within, var_within_cutoff, na.rm = TRUE)
  
  out2 <- out %>%
    dplyr::filter(var_within <= var_within_threshold) %>%
    dplyr::filter(eta2_between >= var_between_cutoff) %>%
    dplyr::filter(!is.na(eta2_between)) %>%
    dplyr::arrange(desc(eta2_between))

  # Store results
  object$features$stats <- out
  object$features$stats_filter <- out2
  object$features$variable <- out2$feature

  # Store parameters
  object$params$FindVariableFeatures <- list(
    assay = assay,
    var_within_cutoff = var_within_cutoff,
    var_between_cutoff = var_between_cutoff,
    n_features_total = nrow(out),
    n_features_variable = nrow(out2),
    var_within_threshold = var_within_threshold
  )

  if (verbose) {
    cat("Variable feature identification completed:\n")
    cat("  Total features:", nrow(out), "\n")
    cat("  Variable features:", nrow(out2), "\n")
    cat("  Within-group variance threshold:", round(var_within_threshold, 4), "\n")
    cat("  Between-group variance threshold:", var_between_cutoff, "\n")
  }

  return(object)

}
