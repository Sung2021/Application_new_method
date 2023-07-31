## 2023 Takeda project
## human: patient cells
setwd('~/Desktop/DF/DFCI_Paweletz/2023_Takeda/')
# Load required libraries
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(tidymodels)

## import data
obj.srt = readRDS('rds/MT30269.30271.23.05.18.rds')

## target cell type
target_cell = 'macrophage'
### obtain Macrophage geneset from GSEA
gmt1 = clusterProfiler::read.gmt(gmtfile = '~/Downloads/GOBP_MACROPHAGE_CYTOKINE_PRODUCTION.v2023.1.Hs.gmt')
gmt1$gene %>% paste0(collapse = ',')
gmt2 = clusterProfiler::read.gmt(gmtfile = '~/Downloads/GOBP_MACROPHAGE_DIFFERENTIATION.v2023.1.Hs.gmt')
gmt2$gene %>% paste0(collapse = ',')
gmt3 = clusterProfiler::read.gmt(gmtfile = '~/Downloads/GOBP_MACROPHAGE_INFLAMMATORY_PROTEIN_1_ALPHA_PRODUCTION.v2023.1.Hs.gmt')
gmt3$gene %>% paste0(collapse = ',')
gmt4 = clusterProfiler::read.gmt(gmtfile = '~/Downloads/GOBP_MACROPHAGE_MIGRATION.v2023.1.Hs.gmt')
gmt4$gene %>% paste0(collapse = ',')
gmt5 = clusterProfiler::read.gmt(gmtfile = '~/Downloads/GOBP_MACROPHAGE_PROLIFERATION.v2023.1.Hs.gmt')
gmt5$gene %>% paste0(collapse = ',')
gmt=rbind(gmt1,gmt2,gmt3,gmt4,gmt5)

target_genes=gmt$gene[!(duplicated(gmt$gene))]
target_genes %>% paste0(collapse = ',')

my_genes <- c("APOE", "ATG7", "BCAT1", "CCL7", "CD163", "CD68", "CD84", "CHI3L1", "CHIT1",
              "CLEC5A", "COL8A2", "COLEC12", "CTSK", "CXCL5", "CYBB", "DNASE2B", "EMP1",
              "FDX1", "FN1", "GM2A", "GPC4", "KAL1", "MARCO", "ME1", "MS4A4A", "MSR1",
              "PCOLCE2", "PTGDS", "RAI14", "SCARB2", "SCG5", "SGMS1", "SULT1C2")

target_genes = union(target_genes, my_genes)
target_genes = target_genes[target_genes %in% rownames(obj.srt)]
# Generate data (cell id information to make split data)
data <- data.frame(cellid = colnames(obj.srt))

# Split the data into 20% training and 80% testing sets
data_split <- initial_split(data, prop = 0.2)

# Extract the training and testing sets
train_cells <- training(data_split)
train_cells[,target_genes] = obj.srt@assays$RNA@counts[target_genes,train_cells$cellid]
rownames(train_cells) = train_cells$cellid
train_data= train_cells[,-1]
columns_to_scale =colnames(train_data)

# Annotate target cells based on gene expression threshold
gene_vector <- c("MARCO", "APOE")
gene_vector = gene_vector[gene_vector %in% colnames(train_data)]
train_data <- train_data %>% 
  mutate(target_cell = ifelse((get(gene_vector[1])> 0 | get(gene_vector[2]) > 0 ),1,0))
train_data$target_cell %>% table()


# Log normalize the predictor variables in the train_data
train_data_log <- train_data %>%
  mutate(across(all_of(columns_to_scale), function(x) {
    if (all(is.numeric(x)) && all(x >= 0)) {
      log2(x + 1)
    } else {
      x  # Keep non-numeric or negative values as they are
    }
  }))

# Create a recipe for XGBoost regression
xgboost_recipe <- recipe(target_cell ~ ., data = train_data_log) %>%
  step_center(all_predictors())

# Create XGBoost regression model
xgboost_model <- boost_tree(mode = "regression", trees = 500, engine = "xgboost")

# Combine the recipe and model into a workflow
xgboost_workflow <- workflow() %>%
  add_recipe(xgboost_recipe) %>%
  add_model(xgboost_model)

# Fit the model to the train_data
xgboost_fit <- fit(xgboost_workflow, data = train_data_log)

# Make predictions on whole data
whole_cells = data.frame(t(as.data.frame(obj.srt@assays$RNA@counts[target_genes,])), check.names = F)
columns_to_scale =colnames(whole_cells)
whole_data = whole_cells

# Log normalize the predictor variables in the train_data
whole_data_log <- whole_data %>%
  mutate(across(all_of(columns_to_scale), function(x) {
    if (all(is.numeric(x)) && all(x >= 0)) {
      log2(x + 1)
    } else {
      x  # Keep non-numeric or negative values as they are
    }
  }))

# Predict using the fitted model
predictions <- predict(xgboost_fit, whole_data_log)

# Add predictions to the whole_data_log data frame
whole_data_log <- cbind(whole_data_log, predictions$.pred)

# Add prediction to Seurat object
obj.srt <- AddMetaData(obj.srt, metadata = whole_data_log$`predictions$.pred`, col.name = 'Macrophage')
obj.srt@meta.data[1:4,]
obj.srt@meta.data %>% ggplot(aes(Macrophage)) + geom_density()

# Visualize the results
FeaturePlot(obj.srt, features = 'Macrophage', cols = c('grey','red'))

# Extract variable importance scores
importance_df <- vip::vip(xgboost_fit, num_features = ncol(train_data_log) - 1)
importance_df$data %>% data.frame()

# View the variable importance scores
print(importance_df)
