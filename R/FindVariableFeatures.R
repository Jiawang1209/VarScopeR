
FindVariableFeatures <- function(mat, group){

  # check input
  stopifnot(ncol(mat) == nrow(group))

  # check input
  stopifnot(identical(colnames(mat), group$Sample))
  stopifnot(dplyr::setequal(colnames(mat), group$Sample))


  # match colnames and Sample information
  idx <- match(colnames(mat), group$Sample)

  if (any(is.na(idx))) {
    missing <- colnames(mat)[is.na(idx)]
    stop("These samples in mat are missing in group$Sample: ",
         paste(missing, collapse = ", "))
  }

  grp <- factor(group$Group[idx])
  groups <- levels(grp)

  if (length(groups) < 2) stop("Need at least 2 groups.")


  # grand mean (per feature)
  grand_mean <- rowMeans(mat)

  # group means (feature x group)  -- force matrix
  group_means <- do.call(
    cbind,
    lapply(groups, function(g) rowMeans(mat[, grp == g, drop = FALSE]))
  )

  colnames(group_means) <- groups

  # weights
  group_sizes <- table(grp)[groups]
  weights <- as.numeric(group_sizes) / sum(group_sizes)

  # between-group variance
  deviation2 <- sweep(group_means, 1, grand_mean, "-")^2
  var_between <- rowSums(deviation2 * rep(weights, each = nrow(mat)))

  # within-group variance
  within_by_group <- do.call(
    cbind,
    lapply(groups, function(g) apply(mat[, grp == g, drop = FALSE], 1, var))
  )

  colnames(within_by_group) <- groups

  var_within <- rowSums(within_by_group * rep(weights, each = nrow(mat)))

  # total
  var_total <- var_within + var_between
  eta2_between <- ifelse(var_total > 0, var_between / var_total, NA)

  # output
  out <- data.frame(
    feature = rownames(mat),
    var_within = var_within,
    var_between = var_between,
    var_total = var_total,
    eta2_between = eta2_between,
    row.names = NULL,
    check.names = FALSE
  )


}
