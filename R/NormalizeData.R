#' Normalize feature abundance data in VarScope object
#'
#' Normalize raw feature abundance data (e.g. OTU/ASV tables) stored in a
#' \code{VarScope} object. The normalized data will be saved as a new assay
#' within the object, while the original data remain unchanged.
#' @param object A \code{VarScope} object
#' containing raw feature abundance data .
#' @param assay Character.
#' Name of the input assay to be normalized
#' @param new.assay Character.
#' Name of the normalized data
#' @param normalization.method Character.
#' Data transformation methods applied before correlation analysis.
#' Options include:
#' "none" (raw data),
#' "scale" (z-score standardization),
#' "center" (mean centering only),
#' "log2" (log2 transfrom),
#' "log10" (log10 transfrom),
#' "ln" (natural transfrom ),
#' "rrarefy" (random rarefaction using \code{vegan::rrarefy}),
#' "rrarefy_relative" (rarefy then convert to relative abundance).
#' @param seed Integer (default = 1115).
#' Random seed for reproducibility.
#'
#' @returns A \code{VarScope} object containing the normalized data,
#' @export
#'
#' @examples NULL
NormalizeData <- function(object,
                          assay = "raw",
                          new.assay = "norm",
                          normalization.method = c("none", "log2", "log10", "ln", "rrarefy", "rrarefy_relative", "clr"),
                          seed = 1115
                          ){

  set.seed(seed)

  # check object
  stopifnot(inherits(object, "VarScope"))

  # argument check
  normalization.method <-  match.arg(normalization.method)

  # get matrix
  mat <- object$assays[[assay]]

  # data transfrom
  mat <- switch (
    normalization.method,
    none = mat,
    log2 = log2(mat + 1),
    log10 = log10(mat + 1),
    ln = log(mat + 1),
    rrarefy = t(vegan::rrarefy(t(mat), min(colSums(mat)))),
    rrarefy_relative = t(vegan::rrarefy(t(mat), min(colSums(mat)))) / colSums(t(vegan::rrarefy(t(mat), min(colSums(mat))))),
    clr = sweep(log(mat + 1e-6), 2, colMeans(log(mat + 1e-6)), "-")
  )

  # out
  object$assays[[new.assay]] <- mat

  # params
  object$params$NormalizeData <- list(
    assay = assay,
    new.assay = new.assay,
    method = normalization.method
  )

  # output
  return(object)

}
