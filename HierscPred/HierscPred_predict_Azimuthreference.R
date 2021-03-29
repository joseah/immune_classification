library(HierscPred)
library(Seurat)
library(SeuratDisk)
library(dplyr)

setwd('/exports/humgen/lmichielsen/sc-eQTLconsortium/datasets/pbmc_multimodal')

load(file = 'trained_ref_noRE.RData') #this loads hier and ref

#load new data
setwd('/exports/humgen/lmichielsen/DataGroningen')

newData <- LoadH5Seurat("testdata.h5Seurat")

newData <- predictTree(hier, newData, threshold = 0.75, recompute_alignment = TRUE)

save(newData, file = 'testdata_annotated_HierscPred.RData')