# Step 1: Install and load the required packages
# install.packages("AUCell")
# install.packages("Seurat")
library(AUCell)
library(Seurat)

# Step 2: Load single-cell RNA-seq data using Seurat
# Replace 'sc_data' with your Seurat object containing the single-cell data
# Make sure your Seurat object contains the expression matrix in the 'assay' slot
sc_data <- Read10X(data.dir = "path_to_your_data")
sc_data <- CreateSeuratObject(counts = sc_data)

# Step 3: Prepare gene sets of interest
# Example gene sets (replace with your own gene sets)
gene_set1 <- c("GENE1", "GENE2", "GENE3", ...)
gene_set2 <- c("GENE4", "GENE5", "GENE6", ...)
# ...

# Step 4: Calculate AUCell scores for each gene set
calculate_aucell_scores_for_gene_set <- function(sc_data, gene_set) {
  # Check if gene set contains valid gene names (removes NAs or duplicates)
  gene_set <- gene_set[!is.na(gene_set) & !duplicated(gene_set)]
  
  # Calculate AUCell scores for the given gene set
  scores <- calculate_aucell_scores(assay(sc_data), gene_set)
  return(scores)
}

# Calculate AUCell scores for gene set 1
scores_gene_set1 <- calculate_aucell_scores_for_gene_set(sc_data, gene_set1)

# Calculate AUCell scores for gene set 2
scores_gene_set2 <- calculate_aucell_scores_for_gene_set(sc_data, gene_set2)

# Repeat for other gene sets if needed
# ...

# Step 5: Explore the results
# Plot the distribution of AUCell scores for gene set 1 with annotations
hist(scores_gene_set1, main = "AUCell Scores for Gene Set 1", xlab = "AUCell Scores", col = "lightblue")

# Identify cells with high enrichment for gene set 2 (e.g., top 10% of cells)
enriched_cells_gene_set2 <- which(scores_gene_set2 >= quantile(scores_gene_set2, 0.9))

# Annotate enriched cells in the Seurat object
sc_data$enriched_cells_gene_set2 <- FALSE
sc_data$enriched_cells_gene_set2[enriched_cells_gene_set2] <- TRUE

# Visualize enriched cells using UMAP or t-SNE plot from Seurat
# Replace 'UMAP' with your preferred dimensionality reduction method
DimPlot(sc_data, group.by = "enriched_cells_gene_set2", reduction = "UMAP")
