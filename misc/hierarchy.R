suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(library(data.tree))
suppressMessages(library(Seurat))

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
    
    # First we do PCA for the root node
    data <- doPCA(data, verbose)
    tree$loadings <- Loadings(data, reduction = 'pca')
    
    for(c in tree$children){
        
        # Cell types
        namesChildren <- as.vector(c$Get('name')) 
        idxChildren <- which(labels %in% namesChildren)
        
        # idxChildren, positive samples (including node itself)
        # rest, negative samples
        # Train classifier here
        
        # If c is not a leaf node, do PCA and continue with children
        if(!isLeaf(c)){
            
            # Get subset of the data
            dataSubset <- subset(data, cells = idxChildren)
            
            # Do PCA on this subset
            dataSubset <- doPCA(dataSubset, verbose)
            c$loadings <- Loadings(dataSubset, reduction = 'pca')
            
            # Continue with the children
            
        }
                
    }
    
    mylist <- list('tree' = tree, 'data' = data)
    return(mylist)
}


doPCA <- function(data, verbose){
    
    if(ncol(data) < 100){
        npcs <- ncol(data) - 2
    } else {
        npcs <- 100
    }
    
    print(npcs)
    
    data <- FindVariableFeatures(data, verbose = verbose)
    data <- ScaleData(data, verbose = verbose)
    data <- RunPCA(data, npcs = npcs, verbose = verbose)
    
    data
}

# small test seuratobject (downsampled 68k)
pbmc <- readRDS('../68K/small68k.rds')
pbmc[['pca']] <- NULL
pbmc[['umap']] <- NULL

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

# 'train' the tree 
test <- trainTree(h, pbmc)
test

t <- test$tree
t$loadings #PCA of all cells
t$Myeloid$loadings #PCA of myeloid cells