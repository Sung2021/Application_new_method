dir='~/Desktop/DF/DFCI_Barbie/DFCI_Barbie_NK_Marco/'
obj.srt = readRDS(paste0(dir,'rds/NK_MN3.MN4.23.07.27.rds'))

# Gene sets to subset out the expression matrix
# NK cell related genes
gene_vector <- c("NCR1", "NCR2", "NCR3", "FCGR3A", "IFNG", "TNF", "KLRC1", "LILRB1", "CD69", "IL2RA", "CCR7", "CXCR3", "PRF1", "GZMA", "GZMH", "GZMK", "GNLY", "FGFBP2",
                 "CD34", "IL2RB", "KLRD1",  "NCAM1", "KIT", "KLRB1")

# Expression data : use data slot
## expression data 
gene.exp <- obj.srt@assays$RNA@data[gene_vector,] %>% data.frame(check.names = F)
gene.exp[1:3,1:3]

# Function to normalize a vector to the range of -1 to 1
normalize_to_range <- function(x) {
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)
  normalized <- -1 + 2 * (x - min_x) / (max_x - min_x)
  return(normalized)
}

# Apply the normalization function to each row
normalized_gene.exp <- gene.exp
for (i in 1:nrow(gene.exp)) {
  normalized_gene.exp[i, ] <- normalize_to_range(gene.exp[i, ])
}

# Display the normalized data frame
normalized_gene.exp[,1:3]

# pheatmap

# row annotation
df.col = obj.srt@meta.data %>% select(orig.ident, cell_type)
df.col = df.col %>% arrange(cell_type, orig.ident)

pheatmap::pheatmap(normalized_gene.exp[,rownames(df.col)], 
                   show_colnames = F,
                   cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("#2874A6", "white", "#D35400"))(1000), 
                   annotation_col = df.col
                   )
