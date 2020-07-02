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
  
  labels <- data@meta.data[[pVar]]
  
  data <- doPCA(data)
  
  #c <- tree$children[[1]]
  flag <- 1
  for(c in tree$children){
    
    if(flag){
      
      # Cell types
      namesChildren <- as.vector(tree$Get('name')) 
      idxChildren <- labels %in% namesChildren
      
      # idxChildren, positive samples (including node itself)
      # rest, negative samples
      # Train classifier here
      response <- vector("character", ncol(data))
      cat("Children: ", c$name, ", ", c$siblings[[1]]$name, "\n", sep = "")
      response[idxChildren] <- c$name
      response[!idxChildren] <- c$siblings[[1]]$name
      data$response <- response
      
      data <- getFeatureSpace(data, pvar = "response")
      data <- trainModel(data)
      model <- data@misc$scPred
      tree$model <- model
    }
    flag <- 0
    # If c is not a leaf node, do PCA and continue with children
    if(!isLeaf(c)){
      
      # Get subset of the data
      dataSubset <- subset(data, cells = Cells(data)[idxChildren])
      
      # Do PCA on this node and continue with children
      c <-  trainNode(c, dataSubset, pVar, verbose)
      
      
      tree$c <-  c
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
#pbmc <- readRDS('../68K/small68k.rds')
data <- readRDS("../immune_prediction/results/2020-04-14_68k_annotation/68k_annotated.RDS")
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

# Iterate over the nodes recursively
predict_tree <- function(tree, template, data){
  
  # Make predictions  
  new <- scPredict(new, tree$model)
  
  # Assign cells to parent node if they are unassigned
  new$scpred_prediction <- if_else(new$scpred_prediction == "unassigned", tree$name, new$scpred_prediction)
  
  for(c in tree$children){
    
    # If c is not a leaf node, continue with predictions
    if(!isLeaf(c)){
      s <- SplitObject(new, "scpred_prediction")
      
      # Get subset of the data
      dataSubset <- subset(data, cells = Cells(data)[idxChildren])
      
      # Do PCA on this node and continue with children
      c <-  predict_tree(c, dataSubset, pVar, verbose)
      
      
      tree$c <-  c
    }
  }
  
  tree    
  
}
