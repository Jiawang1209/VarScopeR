#' Create a VarScope object
#'
#' Create a \code{VarScope} object from an ASV/feature table and corresponding sample metadata.
#' This function initializes the core container used throughout the VarScopeR workflow.
#' @param counts Character
#' The path to a feature table
#' @param meta.data Character
#' @param project Character
#' @param min.group Numeric
#' @param min.features Numeric
#'
#' @returns A \code{VarScope} object containing raw feature counts, sample metadata,
#' @export
#'
#' @examples
#' \dontrun{
#' vsobj <- CreateVarScopeObject(
#'   counts = "inst/extdata/otutab.txt",
#'   meta.data = "inst/extdata/group.csv",
#'   project = "ExampleProject"
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

  # 如果传入的是路径：自动读取 pass the path of file, read it
  if (is.character(counts) && length(counts) == 1 && file.exists(counts)) {
    tmp <- read_counts_and_meta(counts, meta.data, ...)
    counts <- tmp$counts
    group  <- tmp$meta
  } else {
    # pass the  counts is a matrix or data.frame，meta.data is data.frame
    counts <- as.matrix(counts)
    storage.mode(counts) <- "double"
    group <- meta.data
    if (!all(c("Sample","Group") %in% colnames(group))) {
      stop("meta.data must contain columns: Sample and Group.")
    }
  }

  # create a VarScope object
  object <- VarScopeObject()

  # assays
  object$assays$raw <- counts

  # meta data
  object$meta$sample <- group

  # params
  object$params$CreateVarScopeObject <- list(
    n_features = nrow(counts),
    n_samples = ncol(counts),
    groups = table(group$Group)
  )

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


