suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(library(data.tree))
suppressMessages(library(Seurat))
library(scPred)
library(harmony)
library(caret)
source('misc/trainTree.R')
source('misc/predictLabels.R')
source('misc/evaluate.R')

reference <- readRDS('../../citeseq/5k_v3.RDS')
new <- readRDS('../../citeseq/5k_v3_nextgem.RDS')

xx <- reference$cell_type == 'NC/Int Monocyte'
reference$cell_type[xx] <- 'NC-Int Monocyte'

xx <- new$cell_type == 'NC/Int Monocyte'
new$cell_type[xx] <- 'NC-Int Monocyte'

xx <- reference$cell_type != 'CD3+ CD14+ cell'
reference <- reference[,xx]

xx <- new$cell_type != 'CD3+ CD14+ cell'
new <- new[,xx]

# xx <- reference$cell_type != 'gd T cell'
# reference <- reference[,xx]

hier <- c("pbmc/Myeloid/DC/cDC",
          "pbmc/Myeloid/DC/pDC",
          "pbmc/Myeloid/Monocyte/C Monocyte",
          "pbmc/Myeloid/Monocyte/NC-Int Monocyte",
          "pbmc/Lymphoid/B cell",
          "pbmc/Lymphoid/NK and T cell/T cell/CD8+ T cell",
          "pbmc/Lymphoid/NK and T cell/T cell/CD4+ T cell",
          "pbmc/Lymphoid/NK and T cell/T cell/gd T cell",
          "pbmc/Lymphoid/NK and T cell/NK cell",
          "pbmc/Lymphoid/NK and T cell/NKT cell")

# hier <- c("pbmc/Myeloid/DC/cDC",
#           "pbmc/Myeloid/DC/pDC",
#           "pbmc/Myeloid/Monocyte/C Monocyte",
#           "pbmc/Myeloid/Monocyte/NC-Int Monocyte",
#           "pbmc/Lymphoid/B cell",
#           "pbmc/Lymphoid/NK and T cell/T cell/CD8+ T cell",
#           "pbmc/Lymphoid/NK and T cell/T cell/CD4+ T cell",
#           "pbmc/Lymphoid/NK and T cell/NK cell",
#           "pbmc/Lymphoid/NK and T cell/NKT cell")


h <- createTree(hier)


# 'train' the tree 
h <- trainTree(h, reference, pVar = 'cell_type', model = "svmRadial", reconstructionError = TRUE, FN_perc = 0.05)


predictions <- predictTree(h, reference, recompute_alignment = TRUE)

res <- evaluate(new$cell_type, predictions$scpred_prediction, h)

# xx <- new$cell_type != 'gd T cell'
# res <- evaluate(new$cell_type[xx], predictions$scpred_prediction[xx], h)
# 
# xx <- new$cell_type == 'gd T cell'
# table(predictions$scpred_prediction[xx], new$cell_type[xx])

