#' Normalize feature abundance data in VarScope object
#'
#' Normalize raw feature abundance data (e.g. OTU/ASV tables) stored in a
#' \code{VarScope} object. The normalized data will be saved as a new assay
#' within the object, while the original data remain unchanged.
#'
#' @param object A \code{VarScope} object containing raw feature abundance data.
#' @param assay Character. Name of the input assay to be normalized (default: "raw").
#' @param new.assay Character. Name of the normalized data assay (default: "norm").
#' @param normalization.method Character. Data transformation method. Options include:
#'   \itemize{
#'     \item \code{"none"}: No transformation (raw data)
#'     \item \code{"log2"}: Log2 transformation (log2(x + 1))
#'     \item \code{"log10"}: Log10 transformation (log10(x + 1))
#'     \item \code{"ln"}: Natural logarithm (log(x + 1))
#'     \item \code{"rrarefy"}: Random rarefaction to minimum sample depth
#'       (using \code{vegan::rrarefy})
#'     \item \code{"rrarefy_relative"}: Rarefaction followed by relative abundance
#'     \item \code{"clr"}: Centered log-ratio transformation
#'   }
#' @param seed Integer. Random seed for reproducibility (default: 1115).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @returns A \code{VarScope} object with normalized data stored in
#'   \code{object$assays[[new.assay]]}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' vsobj <- CreateVarScopeObject(
#'   counts = "inst/extdata/otutab.txt",
#'   meta.data = "inst/extdata/group.csv"
#' )
#' 
#' # Log2 transformation
#' vsobj <- NormalizeData(vsobj, normalization.method = "log2")
#' 
#' # CLR transformation
#' vsobj <- NormalizeData(vsobj, normalization.method = "clr", new.assay = "clr")
#' 
#' # Rarefaction
#' vsobj <- NormalizeData(vsobj, normalization.method = "rrarefy")
#' }
NormalizeData <- function(object,
                          assay = "raw",
                          new.assay = "norm",
                          normalization.method = c("none", "log2", "log10", "ln", 
                                                   "rrarefy", "rrarefy_relative", "clr"),
                          seed = 1115,
                          verbose = TRUE
                          ){

  set.seed(seed)

  # Check object
  if (!inherits(object, "VarScope")) {
    stop("object must be a VarScope object.")
  }

  # Argument check
  normalization.method <- match.arg(normalization.method)

  # Get matrix
  mat <- object$assays[[assay]]
  if (is.null(mat)) {
    stop("Assay '", assay, "' not found in object.")
  }

  # Validate matrix
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Matrix is empty.")
  }

  if (verbose) {
    cat("Normalizing data using method: '", normalization.method, "'...\n", sep = "")
    cat("  Input assay: '", assay, "'\n", sep = "")
    cat("  Output assay: '", new.assay, "'\n", sep = "")
  }

  # Handle missing values
  n_na_before <- sum(is.na(mat))
  if (n_na_before > 0) {
    warning("Matrix contains ", n_na_before, " NA values. They will be converted to 0.")
    mat[is.na(mat)] <- 0
  }

  # Check for negative values (should not occur in count data)
  if (any(mat < 0, na.rm = TRUE)) {
    warning("Matrix contains negative values. This may cause issues with some normalization methods.")
  }

  # Apply normalization
  mat_normalized <- switch(
    normalization.method,
    none = {
      if (verbose) cat("  No transformation applied.\n")
      mat
    },
    
    log2 = {
      if (verbose) cat("  Applying log2(x + 1) transformation...\n")
      log2(mat + 1)
    },
    
    log10 = {
      if (verbose) cat("  Applying log10(x + 1) transformation...\n")
      log10(mat + 1)
    },
    
    ln = {
      if (verbose) cat("  Applying natural log(x + 1) transformation...\n")
      log(mat + 1)
    },
    
    rrarefy = {
      if (!requireNamespace("vegan", quietly = TRUE)) {
        stop("Package 'vegan' is required for rarefaction. Please install it.")
      }
      if (verbose) cat("  Applying random rarefaction...\n")
      
      # Calculate minimum sample depth
      sample_depths <- colSums(mat)
      min_depth <- min(sample_depths)
      
      if (min_depth == 0) {
        stop("Some samples have zero total counts. Cannot perform rarefaction.")
      }
      
      if (verbose) {
        cat("  Rarefying to depth:", min_depth, "\n")
        cat("  Sample depths range:", min(sample_depths), "-", max(sample_depths), "\n")
      }
      
      # Rarefy (transpose because rrarefy works on rows)
      mat_rare <- t(vegan::rrarefy(t(mat), min_depth))
      mat_rare
    },
    
    rrarefy_relative = {
      if (!requireNamespace("vegan", quietly = TRUE)) {
        stop("Package 'vegan' is required for rarefaction. Please install it.")
      }
      if (verbose) cat("  Applying rarefaction + relative abundance...\n")
      
      sample_depths <- colSums(mat)
      min_depth <- min(sample_depths)
      
      if (min_depth == 0) {
        stop("Some samples have zero total counts. Cannot perform rarefaction.")
      }
      
      if (verbose) {
        cat("  Rarefying to depth:", min_depth, "\n")
      }
      
      # Rarefy
      mat_rare <- t(vegan::rrarefy(t(mat), min_depth))
      
      # Convert to relative abundance
      mat_rare <- sweep(mat_rare, 2, colSums(mat_rare), "/")
      mat_rare[is.na(mat_rare)] <- 0
      mat_rare
    },
    
    clr = {
      if (verbose) cat("  Applying centered log-ratio (CLR) transformation...\n")
      
      # CLR: log(x) - mean(log(x)) per sample (column)
      # Add small pseudocount to avoid log(0)
      pseudocount <- 1e-6
      mat_log <- log(mat + pseudocount)
      mat_clr <- sweep(mat_log, 2, colMeans(mat_log), "-")
      mat_clr
    }
  )

  # Validate output
  if (any(is.infinite(mat_normalized), na.rm = TRUE)) {
    warning("Normalized matrix contains infinite values.")
    mat_normalized[is.infinite(mat_normalized)] <- NA
  }

  # Store normalized data
  object$assays[[new.assay]] <- mat_normalized

  # Store parameters
  object$params$NormalizeData <- list(
    assay = assay,
    new.assay = new.assay,
    method = normalization.method,
    n_features = nrow(mat_normalized),
    n_samples = ncol(mat_normalized)
  )

  if (verbose) {
    cat("Normalization completed.\n")
    cat("  Output matrix dimensions:", nrow(mat_normalized), "×", ncol(mat_normalized), "\n")
    if (normalization.method != "none") {
      cat("  Value range:", round(min(mat_normalized, na.rm = TRUE), 4), "-",
          round(max(mat_normalized, na.rm = TRUE), 4), "\n")
    }
  }

  return(object)

}
