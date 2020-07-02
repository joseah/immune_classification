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


DimPlot(reference, group.by = "CellType")
DimPlot(new, group.by = "CellType")


# Integrate data ----------------------------------------------------------

library(harmony)

# Merging all data an creating joint PCA between reference and new dataset
data <- merge(reference, new)
data <- process_data(data)

# See unintegrated data
DimPlot(data, group.by = "Method")

# Run harmony with "free" (no reference-based) alignment
data <- RunHarmony(data, "Method")

# Run UMAP on integrated data
data <- RunUMAP(data, reduction = "harmony", 
                dims = 1:30, reduction.key = "umaph_", 
                reduction.name = "umap_harmony")

# See integrated data
DimPlot(data, group.by = "Method", reduction = "umap_harmony")


# Split data into original reference and new datasets 
# but now with integrated data in a new space

tt  <- SplitObject(data, "Method")
reference_2 <- tt$`10x Chromium (v2) A`
new_2 <- tt$`10x Chromium (v2) B`

# Visualize aligned training data (without new/test data)
DimPlot(reference_2, group.by = "CellType")

# Get feature space
reference_2 <- getFeatureSpace(reference_2, "CellType", reduction = "harmony")

# Train models
reference_2 <- trainModel(reference_2)


# In this case, new_2[["harmony"]] contains the aligned data. Store
# this as a Seurat reduction (named "scpred") in the new/test dataset
# and ask scPred to use this cell embeddings instead with
# recompute_alignment = FALSE. This will skip data scaling
# and loading projection and use the provided cell embeddings
# instead.
# First, create a Seurat dimred object via
# CreateDimReducObject() storing the independently aligned data
# and name the reduction as "scPred".
# NOTE: the name of the features should be the same as in the reference so
# the model can be applied (e.g "PC_", "harmony_"). Make sure the 
# reference@misc$scPred@reduction_key value is the same as the column
# names in the new/test aligned-embeddings columns

reference_2@misc$scPred@reduction
new_2[["harmony"]]@cell.embeddings %>% colnames()

dimred <- CreateDimReducObject(embeddings = new_2[["harmony"]]@cell.embeddings, 
                                          key = reference_2@reductions$harmony@key,
                                          assay = DefaultAssay(new_2)
                                          )
# In this case, remove original alignment as it is now stored already in "dimred".
new_2[["harmony"]] <- NULL

# Store aligned data as "scpred"
new_2[["scpred"]] <- dimred


# Make predictions using the "scpred" cell embeddings with recompute_alignment = FALSE
new_2 <- scPredict(new_2, reference_2, recompute_alignment = FALSE)


# Make umap from aligned data

new_2 <- RunUMAP(new_2, dims = 1:30, reduction = "scpred", reduction.name = "scpred_umap", reduction.key = "scpredumap_")

DimPlot(new_2, group.by = "scpred_prediction", reduction = "scpred_umap")


