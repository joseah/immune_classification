library(HierscPred)
library(Seurat)
library(SeuratDisk)

setwd('/exports/humgen/lmichielsen/sc-eQTLconsortium/datasets/pbmc_multimodal')

ref <- LoadH5Seurat("pbmc_multimodal.h5seurat")

#Remove platelet, eryth, doublet 
idx = which(ref@meta.data$celltype.l2 != 'Doublet')
ref <- ref[,idx]

h <- c("pbmc/Myeloid/DC/cDC1",
       "pbmc/Myeloid/DC/cDC2",
       "pbmc/Myeloid/DC/ASDC",
       "pbmc/Myeloid/DC/pDC",
       "pbmc/Myeloid/Monocyte/CD14 Mono",
       "pbmc/Myeloid/Monocyte/CD16 Mono",
       "pbmc/HSPC",
       "pbmc/Platelet",
       "pbmc/Eryth",
       "pbmc/Lymphoid/ILC",
       "pbmc/Lymphoid/B cell/B intermediate",
       "pbmc/Lymphoid/B cell/B naive",
       "pbmc/Lymphoid/B cell/B memory",
       "pbmc/Lymphoid/B cell/Plasmablast",
       "pbmc/Lymphoid/NK cell/NK",
       "pbmc/Lymphoid/NK cell/NK_CD56bright",
       "pbmc/Lymphoid/NK cell/NK Proliferating",
       "pbmc/Lymphoid/T cell/CD8+ T cell/CD8 TCM",
       "pbmc/Lymphoid/T cell/CD8+ T cell/CD8 Naive",
       "pbmc/Lymphoid/T cell/CD8+ T cell/CD8 TEM",
       "pbmc/Lymphoid/T cell/CD8+ T cell/CD8 Proliferating",
       "pbmc/Lymphoid/T cell/CD4+ T cell/Treg",
       "pbmc/Lymphoid/T cell/CD4+ T cell/CD4 Proliferating",
       "pbmc/Lymphoid/T cell/CD4+ T cell/CD4 Naive",
       "pbmc/Lymphoid/T cell/CD4+ T cell/CD4 TCM",
       "pbmc/Lymphoid/T cell/CD4+ T cell/CD4 CTL",
       "pbmc/Lymphoid/T cell/CD4+ T cell/CD4 TEM",
       "pbmc/Lymphoid/T cell/gdT",
       "pbmc/Lymphoid/T cell/MAIT",
       "pbmc/Lymphoid/T cell/dnT")


h2 <- c("pbmc/Myeloid/DC/cDC1",
       "pbmc/Myeloid/DC/cDC2",
       "pbmc/Myeloid/DC/ASDC",
       "pbmc/Myeloid/DC/pDC",
       "pbmc/Myeloid/Monocyte/CD14 Mono",
       "pbmc/Myeloid/Monocyte/CD16 Mono",
       "pbmc/HSPC",
       "pbmc/Platelet",
       "pbmc/Eryth",
	   "pbmc/Plasmablast",
	   "pbmc/B cell/B intermediate",
       "pbmc/B cell/B naive",
       "pbmc/B cell/B memory",
	   "pbmc/proliferating/NK Proliferating",
	   "pbmc/proliferating/CD8 Proliferating",
	   "pbmc/proliferating/CD4 Proliferating",
       "pbmc/Lymphoid/MAIT",
       "pbmc/Lymphoid/Group1/NK",
       "pbmc/Lymphoid/Group1/NK_CD56bright",
	   "pbmc/Lymphoid/Group1/ILC",
       "pbmc/Lymphoid/Group1/Group1_T/gdT",
       "pbmc/Lymphoid/Group1/Group1_T/CD8 TEM",
       "pbmc/Lymphoid/Group1/Group1_T/CD4 CTL",
	   "pbmc/Lymphoid/Group2/Group2_naive/CD8 Naive",
       "pbmc/Lymphoid/Group2/Group2_naive/Treg",
       "pbmc/Lymphoid/Group2/Group2_naive/CD4 Naive",
	   "pbmc/Lymphoid/Group2/Group2_memory/CD8 TCM",
       "pbmc/Lymphoid/Group2/Group2_memory/CD4 TCM",
       "pbmc/Lymphoid/Group2/Group2_memory/CD4 TEM",
       "pbmc/Lymphoid/Group2/dnT")


hier <- create_hierarchy(h2)
hier <- train_tree(ref, hier, pvar = 'celltype.l2', reconstruction_error = FALSE)

save(hier, ref, file = 'trained_ref_noRE.RData')