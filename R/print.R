#' Print VarScope Object
#'
#' @param x VarScope Object
#' @param ...
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
print.VarScope <- function(x, ...) {

  cat("An object of class VarScope\n")
  cat(strrep("-", 33), "\n")

  ## ---- basic size ----
  if (!is.null(x$assays$raw)) {
    n_feat <- nrow(x$assays$raw)
    n_samp <- ncol(x$assays$raw)
    cat(n_feat, "features across", n_samp, "samples\n")
  }

  ## ---- assays ----
  assays_present <- names(x$assays)
  if (length(assays_present) > 0) {
    cat("Assays:", paste(assays_present, collapse = ", "), "\n")
  }

  ## ---- variable features ----
  n_var <- if (!is.null(x$features$variable)) {
    length(x$features$variable)
  } else {
    0
  }
  cat("Variable features:", n_var, "\n")

  ## ---- reductions ----
  if (length(x$reductions) > 0) {
    cat("Reductions:", paste(names(x$reductions), collapse = ", "), "\n")
  }

  ## ---- clusters ----
  if (!is.null(x$clusters$module)) {
    cat("Feature modules:",
        length(unique(x$clusters$module)), "\n")
  }

  invisible(x)
}
