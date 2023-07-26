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

# Create an R vector with the gene names
gene_vector1 <- c("APOE", "ATG7", "BCAT1", "CCL7", "CD163", "CD68", "CD84", "CHI3L1", "CHIT1",
                 "CLEC5A", "COL8A2", "COLEC12", "CTSK", "CXCL5", "CYBB", "DNASE2B", "EMP1",
                 "FDX1", "FN1", "GM2A", "GPC4", "KAL1", "MARCO", "ME1", "MS4A4A", "MSR1",
                 "PCOLCE2", "PTGDS", "RAI14", "SCARB2", "SCG5", "SGMS1", "SULT1C2")

gene_vector2 <- c("AIF1", "CCL1", "CCL14", "CCL23", "CCL26", "CD300LB", "CNR1", "CNR2",
                 "EIF1", "EIF4A1", "FPR1", "FPR2", "FRAT2", "GPR27", "GPR77", "RNASE2",
                 "MS4A2", "BASP1", "IGSF6", "HK3", "VNN1", "FES", "NPL", "FZD2", "FAM198B",
                 "HNMT", "SLC15A3", "CD4", "TXNDC3", "FRMD4A", "CRYBB1", "HRH1", "WNT5B")

gene_vector3 <- c("FUCA1", "MMP9", "LGMN", "HS3ST2", "TM4SF19", "CLEC5A", "GPNMB", "C11orf45", "CD68", "CYBB")

gene_vector4 <- c("CD163", "CD14", "CSF1R", "C1QC", "VSIG4", "C1QA", "FCER1G", "F13A1", "TYROBP", "MSR1",
                 "C1QB", "MS4A4A", "FPR1", "S100A9", "IGSF6", "LILRB4", "FPR3", "SIGLEC1", "LILRA1",
                 "LYZ", "HK3", "SLC11A1", "CSF3R", "CD300E", "PILRA", "FCGR3A", "AIF1", "SIGLEC9",
                 "FCGR1C", "OLR1", "TLR2", "LILRB2", "C5AR1", "FCGR1A", "MS4A6A", "C3AR1", "HCK",
                 "IL4I1", "LST1", "LILRA5", "CSTA", "IFI30", "CD68", "TBXAS1", "FCGR1B", "LILRA6",
                 "CXCL16", "NCF2", "RAB20", "MS4A7", "NLRP3", "LRRC25", "ADAP2", "SPP1", "CCR1",
                 "TNFSF13", "RASSF4", "SERPINA1", "MAFB", "IL18", "FGL2", "SIRPB1", "CLEC4A", "MNDA",
                 "FCGR2A", "CLEC7A", "SLAMF8", "SLC7A7", "ITGAX", "BCL2A1", "PLAUR", "SLCO2B1",
                 "PLBD1", "APOC1", "RNF144B", "SLC31A2", "PTAFR", "NINJ1", "ITGAM", "CPVL", "PLIN2",
                 "C1orf162", "FTL", "LIPA", "CD86", "GLUL", "FGR", "GK", "TYMP", "GPX1", "NPL", "ACSL1")

intersect(gene_vector1,gene_vector4)
macro_genes = union(gene_vector1, union(gene_vector2, union(gene_vector3, gene_vector4)))
macro_genes = macro_genes[macro_genes %in% rownames(obj.srt)]

# generate data
data <- data.frame(cellid=colnames(obj.srt))
data[1:3,]

# Split the data into 20% training and 80% testing sets
data_split <- initial_split(data, prop = 0.2)

# Extract the training and testing sets
train_cells <- training(data_split)
test_cells <- testing(data_split)

train_data = obj.srt@assays$RNA@counts[macro_genes,train_cells$cellid] %>% t() %>% data.frame(check.names = F)
columns_to_scale =colnames(train_data)

# add macrophage annotation manually
macro_short=c('CD163', 'CD68', 'MARCO', 'MRC1','APOE')
macro_short = macro_short[macro_short %in% colnames(train_data)]

# train_data = train_data %>% mutate(macrophage=ifelse((rowSums(train_data) >= quantile(rowSums(train_data), probs = 0.75)), 'macrophage', 'other'))
train_data = train_data %>% mutate(macrophage=ifelse((CD163+CD68+CSF1R+MS4A4A+MARCO+APOE >=10), 'macrophage', 'other'))
train_data$macrophage %>% table()

# Log normalize the predictor variables in the train_data
train_data_log <- train_data %>%
  mutate(across(all_of(columns_to_scale), function(x) {
    if (all(is.numeric(x)) && all(x >= 0)) {
      log2(x + 1)
    } else {
      x  # Keep non-numeric or negative values as they are
    }
  }))

## import libarary
library(tidymodels)

# Create a recipe for XGBoost classification
xgboost_recipe <- recipe(macrophage ~ ., data = train_data_log) %>%
  step_scale(all_predictors()) %>%
  step_center(all_predictors())

# Create XGBoost classification model
xgboost_model <- boost_tree(mode = "classification", trees = 100, engine = "xgboost")

# Combine the recipe and model into a workflow
xgboost_workflow <- workflow() %>%
  add_recipe(xgboost_recipe) %>%
  add_model(xgboost_model)

# Fit the model to the train_data
xgboost_fit <- fit(xgboost_workflow, data = train_data_log)

# Print the model summary
summary(xgboost_fit$fit)

# Make predictions on whole data
whole_data = obj.srt@assays$RNA@counts[macro_genes,] %>% t() %>% data.frame(check.names = F)
columns_to_scale =colnames(whole_data)

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

# Print the predictions
print(predictions)

# Add predictions to the data_log data frame
whole_data_log <- cbind(whole_data_log, predictions$.pred_class)
whole_data_log[1:3, ]
summary(whole_data_log$`predictions$.pred_class`)
tmp = whole_data_log$`predictions$.pred_class` %>% as.data.frame(row.names = rownames(whole_data_log))
tmp$. %>% table()

# Add prediction to seurat object
obj.srt=AddMetaData(obj.srt, metadata =tmp, col.name = 'Macrophage')

DimPlot(obj.srt, group.by = 'Macrophage', cols = c('red','grey'))

# Extract variable importance scores
importance_df <- vip::vip(xgboost_fit, num_features = ncol(train_data_log) - 1)

# View the variable importance scores
print(importance_df)

importance_df$data %>% data.frame()
