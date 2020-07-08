suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(library(data.tree))
suppressMessages(library(Seurat))
library(scPred)

createTree <- function(stringTree){
  
  # Create dataframe
  h <- tibble(pathString = stringTree)
  
  # Convert dataframe to a tree
  h <- as.Node(h)
  
  h
}

reference <- readRDS('../immune_prediction/results/2020-06-01_5k_v3/5k_v3.RDS')
new <- readRDS('../immune_prediction/results/2020-06-01_5k_v3_nextgem/5k_v3_nextgem.RDS')


reference[['pca']] <- NULL
new[['umap']] <- NULL


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

reference$labels <- if_else(reference$cell_type == "NC/Int Monocyte", 
                            "NC-Int Monocyte", 
                            reference$cell_type)

i <- reference$labels != "CD3+ CD14+ cell"
reference <- reference[,i]
h <- createTree(hier)
DefaultAssay(reference) <- "RNA"

new$labels <- if_else(new$cell_type == "NC/Int Monocyte", "NC-Int Monocyte", new$cell_type)
i <- new$labels != "CD3+ CD14+ cell"
new <- new[,i]






refData <- CreateSeuratObject(reference@assays$RNA@counts, 
                   meta.data = reference@meta.data[, "labels", drop = FALSE])

newData <- CreateSeuratObject(new@assays$RNA@counts, 
                              meta.data = new@meta.data[, "labels", drop = FALSE])


refData$dataset <- "reference"
newData$dataset <- "new"

# tree <- h
# pVar <- "labels"
# reduction <- "pca"



# Iterate over the nodes recursively
trainNode <- function(tree, refData, newData, pVar, reduction = "pca", verbose){
  
  cat("Training parent node: ", tree$name, "\n", sep = "")
  
  cat("Integrating data...\n")
  data <- list(refData, newData)
  data <- lapply(X = data, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
  })
  features <- SelectIntegrationFeatures(object.list = data)
  data <- lapply(X = data, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  anchors <- FindIntegrationAnchors(object.list = data, 
                                    reference = 1, 
                                    reduction = "rpca", 
                                    dims = 1:30)
  
  integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, verbose = FALSE)
  integrated <- SplitObject(integrated, "dataset")
  data <- integrated$reference
  test <- integrated$new
  cat("Data succesfully integrated...\n")
  
  labels <- data@meta.data[[pVar]]
  newLabels <- labels
  
  ## First rewrite the labels 
  for(c in tree$children){
    
    namesChildren <- as.vector(c$Get('name')) 
    idxChildren <- labels %in% namesChildren
    
    newLabels[idxChildren] <- c$name
    
  }
  
  data$response = newLabels
  
  cat("Cell types in this layer of the tree:", unique(newLabels), "\n", sep = " ")

  ## get informative PCs and train classifier
  data <- getFeatureSpace(data, pvar = "response")
  data <- trainModel(data)
  
  cat("Model trained\n")
  
  
  cat("Classifying cells\n")
  
  names(test@reductions) <- "scpred"
  
  test <- scPredict(test, data, recompute_alignment = FALSE)
  test$scpred_prediction <- if_else(test$scpred_prediction == "unassigned", tree$name, test$scpred_prediction)
  

  ## Continue to children if they are not a leaf
  for(c in tree$children){
    
    if(!isLeaf(c)){
      namesChildren <- as.vector(c$Get('name')) 
      idxChildren <- labels %in% namesChildren
      
      # Get subset of the data
      dataSubset <- subset(data, cells = Cells(data)[idxChildren])
      
      # test
      idxChildren <- test$scpred_prediction %in% namesChildren
      testSubset <- subset(test, cells = Cells(test)[idxChildren])
      
      DefaultAssay(dataSubset) <- "RNA"
      DefaultAssay(testSubset) <- "RNA"
      
      i <- names(new@reductions) == "scpred"
      if(any(i)){
        testSubset@reductions[which(i)] <- NULL
      }
      
      # Do PCA on this node and continue with children
      testSubset  <- trainNode(c, dataSubset, testSubset, pVar = pVar, reduction = reduction, verbose = verbose)
      test$scpred_prediction[idxChildren] <- testSubset$scpred_prediction
      
    }
    
  }
  
  test    
  
}


pred <- trainNode(h, refData, newData, pVar = 'labels', reduction = "pca", verbose = FALSE)

res <- table(pred$scpred_prediction == pred$labels)
{res[2] / sum(res) * 100} %>% round(2)


