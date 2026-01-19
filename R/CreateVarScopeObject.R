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

CreateVarScopeObject <- function(counts,
                                 meta.data,
                                 project = NULL,
                                 min.group = 0,
                                 min.features = 0
                                 ){

  counts <- read.table(file = counts, header = T, sep = "\t", row.names = 1)
  group <- read.csv(file = meta.data, header = T) %>%
    purrr::set_names(c("Sample", "Group"))

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
