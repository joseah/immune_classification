# Import libraries --------------------------------------------------------

library(scPred)
library(devtools)
library(Seurat)
library(tidyverse)
library(SeuratData)


# Read and process data ---------------------------------------------------

pbmcsca <- SeuratData::LoadData("pbmcsca")
all_data <- SplitObject(pbmcsca, "Method")

all_data <- all_data[c("10x Chromium (v2) A", "10x Chromium (v2) B")]

reference <- all_data$`10x Chromium (v2) A`
new <- all_data$`10x Chromium (v2) B`

process_data <- function(x, n = 30){
  x %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    RunUMAP(dims = 1:n)
}

reference <- process_data(reference)
new <- process_data(new)

reference <- getFeatureSpace(reference, "CellType")
reference <- trainModel(reference)

new <- scPredict(new, model)


DimPlot(new, group.by = "scpred_no_rejection")

