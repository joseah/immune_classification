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

# tree = tree we want to tree, we construct this tree before we call this function?
# data = seuratobject of the reference to use for training
# pVar = name of the column of the metadata used as labels

trainTree <- function(tree, data, pVar = 'cell_type', verbose = FALSE){
  
  labels <- data@meta.data[[pVar]]
  
  tree <-  trainNode(tree, data, pVar, verbose)
  tree
}


# Iterate over the nodes recursively
trainNode <- function(tree, data, pVar, verbose){
  
  cat("Training parent node: ", tree$name, "\n", sep = "")
  
  data <- doPCA(data)
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
  model <- data@misc$scPred
  tree$model <- model
  
  ## Continue to children if they are not a leaf
  for(c in tree$children){
    
    if(!isLeaf(c)){
      namesChildren <- as.vector(c$Get('name')) 
      idxChildren <- labels %in% namesChildren
      
      # Get subset of the data
      dataSubset <- subset(data, cells = Cells(data)[idxChildren])
      
      # Do PCA on this node and continue with children
      trainNode(c, dataSubset, pVar, verbose)
      
    }
    
  }
  
  tree    
  
}


doPCA <- function(data, verbose = FALSE){
  
  if(ncol(data) < 50){
    npcs <- ncol(data) - 2
  } else {
    npcs <- 50
  }
  
  data <- FindVariableFeatures(data, verbose = verbose)
  data <- ScaleData(data, verbose = verbose)
  data <- RunPCA(data, npcs = npcs, verbose = verbose)
  
  data
}


# small test seuratobject (downsampled 68k)
data <- readRDS('../68K/68k_corr.RDS')
# data <- readRDS("../immune_prediction/results/2020-04-14_68k_annotation/68k_annotated.RDS")
set.seed(66)
i <- sample(seq_len(nrow(data)), size = 3000)
data <- data[,i]


data[['pca']] <- NULL
data[['umap']] <- NULL

# create the tree
hier <- c("pbmc/Myeloid/DC/mDCs",
          "pbmc/Myeloid/DC/pDC",
          "pbmc/Myeloid/Monocyte/CD16+ Monocytes",
          "pbmc/Myeloid/Monocyte/CD14+ Monocytes",
          "pbmc/Lymphoid/B cells",
          "pbmc/Lymphoid/Plasma",
          "pbmc/Lymphoid/NKT cell/T cell/CD8+ T cell",
          "pbmc/Lymphoid/NKT cell/T cell/CD4+ T cell",
          "pbmc/Lymphoid/NKT cell/Natural Killer")

h <- createTree(hier)

DefaultAssay(data) <- "RNA"

# 'train' the tree 
tree <- trainTree(h, data)


new <- data
tree <- test
template <- Clone(tree)


predictTree <- function(tree, template, newData){
  
  newData <- predictNode(tree, template, newData)
  
  newData
  
}

# Iterate over the nodes recursively
predictNode <- function(tree, template, newData){
  
  # Make predictions  
  newData <- scPredict(newData, tree$model)
  
  # Assign cells to parent node if they are unassigned
  newData$scpred_prediction <- if_else(newData$scpred_prediction == "unassigned", tree$name, newData$scpred_prediction)
  
  for(c in tree$children){
    
    # If c is not a leaf node, continue with predictions
    if(!isLeaf(c)){
      # s <- SplitObject(new, "scpred_prediction")
      
      # Get subset of the data
      namesChildren <- as.vector(c$Get('name')) 
      idxChildren <- newData$scpred_prediction %in% namesChildren
      dataSubset <- subset(newData, cells = Cells(newData)[idxChildren])
      
      # Predict the labels of these cells
      dataSubset <- predictNode(c, template, dataSubset)
      
      newData$scpred_prediction[idxChildren] <- dataSubset$scpred_prediction 
      
    }
  }
  
  newData    
  
}
