# VarScopeR

VarScopeR is an R package for identifying perturbation-driven feature modules by decomposing within- and between-group variance. The package integrates feature-centric variance decomposition, dimensionality reduction, and graph-based community detection to reveal coherent response patterns across experimental conditions.

## Installation

```r
# Install from source (development version)
devtools::install_github("yourusername/VarScopeR")

# Or install from local source
devtools::install("path/to/VarScopeR")
```

## Dependencies

### Required Packages
- `dplyr`: Data manipulation
- `magrittr`: Pipe operator
- `readxl`: Reading Excel files
- `vegan`: Ecological analysis (rarefaction)

### Optional Packages (for specific features)
- `uwot`: UMAP dimensionality reduction
- `igraph`: Graph-based clustering

## Quick Start

```r
library(VarScopeR)

# 1. Create VarScope object
vsobj <- CreateVarScopeObject(
  counts = "path/to/feature_table.txt",
  meta.data = "path/to/group.csv"
)

# 2. Normalize data
vsobj <- NormalizeData(vsobj, normalization.method = "log2")

# 3. Find variable features
vsobj <- FindVariableFeatures(vsobj)

# 4. Scale data
vsobj <- ScaleData(vsobj)

# 5. Run PCA
vsobj <- RunPCA(vsobj, npcs = 10)

# 6. Run UMAP (optional)
vsobj <- RunUMAP(vsobj, dims = 1:10)

# 7. Find neighbors and clusters
vsobj <- FindNeighbors(vsobj, dims = 1:10)
vsobj <- FindClusters(vsobj, resolution = 0.5)
```

## Main Functions

### 1. CreateVarScopeObject

Create a VarScope object from feature table and metadata.

```r
vsobj <- CreateVarScopeObject(
  counts = "inst/extdata/otutab.txt",  # Path to feature table
  meta.data = "inst/extdata/group.csv" # Path to metadata
)
```

**Supported file formats:**
- CSV, TSV, TXT files
- Excel files (.xlsx, .xls)

**Metadata format:**
- Must contain columns: `Sample` and `Group`

### 2. NormalizeData

Normalize feature abundance data using various transformation methods.

```r
vsobj <- NormalizeData(
  object = vsobj,
  assay = "raw",
  new.assay = "norm",
  normalization.method = "log2"  # Options: "none", "log2", "log10", "ln", 
                                  # "rrarefy", "rrarefy_relative", "clr"
)
```

**Available methods:**
- `none`: Raw data
- `log2`, `log10`, `ln`: Logarithmic transformations
- `rrarefy`: Random rarefaction (using `vegan::rrarefy`)
- `rrarefy_relative`: Rarefaction followed by relative abundance conversion
- `clr`: Centered log-ratio transformation

### 3. FindVariableFeatures

Identify perturbation-driven variable features by decomposing variance into within-group and between-group components.

```r
vsobj <- FindVariableFeatures(
  object = vsobj,
  assay = "norm",
  var_within_cutoff = 0.4,   # Quantile threshold for within-group variance
  var_between_cutoff = 0.5   # Lower threshold for between-group variance
)
```

**Output:**
- `object$features$stats`: Statistics for all features
- `object$features$stats_filter`: Filtered variable features

### 4. ScaleData

Scale (z-score standardize) feature data for PCA analysis.

```r
vsobj <- ScaleData(
  object = vsobj,
  assay = "norm",
  new.assay = "scale",
  features = NULL  # NULL uses variable features if available
)
```

### 5. RunPCA

Perform Principal Component Analysis on ASV/feature data.

```r
vsobj <- RunPCA(
  object = vsobj,
  assay = "scale",
  npcs = 10,
  reduction.name = "pca"
)
```

**Access results:**
```r
# ASV loadings (ASVs × PCs)
loadings <- vsobj$reductions$pca$loadings

# Sample scores (samples × PCs)
scores <- vsobj$reductions$pca$scores

# Variance explained
var_explained <- vsobj$reductions$pca$var.explained
```

### 6. RunUMAP

Perform UMAP dimensionality reduction for visualization.

```r
vsobj <- RunUMAP(
  object = vsobj,
  reduction = "pca",      # Use PCA results
  dims = 1:10,            # Which PCs to use
  n.neighbors = 30,
  min.dist = 0.3,
  n.components = 2
)
```

**Note:** Requires `uwot` package: `install.packages("uwot")`

**Access results:**
```r
# UMAP coordinates (samples × UMAP dimensions)
umap_coords <- vsobj$reductions$umap$embedding
```

### 7. FindNeighbors

Construct K-nearest neighbor (KNN) and shared nearest neighbor (SNN) graphs.

```r
vsobj <- FindNeighbors(
  object = vsobj,
  reduction = "pca",
  dims = 1:10,
  k.param = 20,
  compute.SNN = TRUE
)
```

**Output:**
- `object$graphs$pca_knn`: KNN graph
- `object$graphs$pca_snn`: SNN graph

### 8. FindClusters

Identify clusters using graph-based community detection.

```r
vsobj <- FindClusters(
  object = vsobj,
  graph.name = NULL,      # Auto-detect SNN graph
  resolution = 0.5,      # Higher = more clusters
  algorithm = "louvain"   # Options: "louvain", "leiden"
)
```

**Note:** For Leiden algorithm, requires `igraph` package.

**Access results:**
```r
# Cluster assignments
clusters <- vsobj$clusters$cluster_pca_snn_res0.5
```

## Complete Workflow Example

```r
library(VarScopeR)

# Step 1: Create object
vsobj <- CreateVarScopeObject(
  counts = "inst/extdata/otutab.txt",
  meta.data = "inst/extdata/group.csv",
  project = "MyProject"
)

# Step 2: Normalize
vsobj <- NormalizeData(
  vsobj,
  normalization.method = "log2"
)

# Step 3: Find variable features
vsobj <- FindVariableFeatures(
  vsobj,
  var_within_cutoff = 0.4,
  var_between_cutoff = 0.5
)

# Step 4: Scale data
vsobj <- ScaleData(vsobj)

# Step 5: PCA
vsobj <- RunPCA(vsobj, npcs = 10)

# Access PCA results
pc_loadings <- vsobj$reductions$pca$loadings  # ASVs × PCs
pc_scores <- vsobj$reductions$pca$scores       # Samples × PCs

# Step 6: UMAP (optional)
vsobj <- RunUMAP(vsobj, dims = 1:10)
umap_coords <- vsobj$reductions$umap$embedding

# Step 7: Clustering
vsobj <- FindNeighbors(vsobj, dims = 1:10)
vsobj <- FindClusters(vsobj, resolution = 0.5)

# View results
print(vsobj)
```

## Data Structure

The VarScope object contains:

- **`assays`**: Data matrices (raw, norm, scale)
- **`meta`**: Sample metadata
- **`features`**: Feature statistics and variable features
- **`reductions`**: Dimensionality reduction results (PCA, UMAP)
- **`graphs`**: KNN and SNN graphs
- **`clusters`**: Cluster assignments
- **`params`**: Parameters used in each step
- **`misc`**: Miscellaneous data

## Key Features

1. **Variance Decomposition**: Identify features with high between-group and low within-group variance
2. **Dimensionality Reduction**: PCA and UMAP for visualization and analysis
3. **Graph-based Clustering**: Community detection using Louvain/Leiden algorithms
4. **Flexible Data Input**: Supports multiple file formats (CSV, TSV, Excel)
5. **Seurat-like Workflow**: Familiar interface for users of Seurat package

## Citation

If you use VarScopeR in your research, please cite:

```
Liu, Y. (2024). VarScopeR: Identifying perturbation-driven feature modules 
by variance decomposition. R package version 0.1.2.
```

## Author

**Yue Liu**
- Email: yueliu@iae.ac.cn

## License

This package is licensed under the MIT License. See `LICENSE` file for details.

## Version

Current version: 0.1.2

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Acknowledgments

This package is inspired by the Seurat package for single-cell RNA-seq analysis, adapted for microbiome and feature table analysis.
