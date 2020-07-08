suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(library(data.tree))
suppressMessages(library(Seurat))
library(scPred)

reference <- readRDS('../immune_prediction/results/2020-06-01_5k_v3/5k_v3.RDS')
new <- readRDS('../immune_prediction/results/2020-06-01_5k_v3_nextgem/5k_v3_nextgem.RDS')

i <- reference$cell_type != "CD3+ CD14+ cell"
reference <- reference[,i]
DefaultAssay(reference) <- "RNA"

i <- new$cell_type != "CD3+ CD14+ cell"
new <- new[,i]

reference <- getFeatureSpace(reference, "cell_type")
reference <- trainModel(reference)
new <- scPredict(new, reference)

crossTab(new, "cell_type", "scpred_prediction")


res <- table(new$scpred_prediction == new$cell_type)

res / sum(res)


reference$dataset <- "5k_v3"
new$dataset <- "5k_v3_nextgem"

data <- merge(reference, new)
data <- data %>%  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

set.seed(66)
data <- harmony::RunHarmony(data, "dataset")
data <- RunUMAP(data, reduction = "harmony", dims = 1:30)
DimPlot(data, group.by = "dataset") + DimPlot(data, group.by = "cell_type") 



#### Seurat


data.list <- SplitObject(data, split.by = "dataset")
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = data.list, 
                                  reference = c(1, 2), 
                                  reduction = "rpca", 
                                  dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:30)
DimPlot(integrated, group.by = "dataset") + DimPlot(integrated, group.by = "cell_type") 







