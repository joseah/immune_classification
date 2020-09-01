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

trainTree <- function(tree, data, pVar = 'cell_type', reduction = 'pca', verbose = FALSE, modelName = 'svmRadial',reconstructionError = TRUE){
  
  labels <- data@meta.data[[pVar]]
  
  tree <-  trainNode(tree, data, pVar, reduction, verbose, modelName, reconstructionError)
  tree
}


# Iterate over the nodes recursively
trainNode <- function(tree, data, pVar, reduction, verbose, modelName, reconstructionError){
  
  cat("Training parent node: ", tree$name, "\n", sep = "")
  
  labels <- data@meta.data[[pVar]]
  newLabels <- labels
  
  # if(reconstructionError){
  #   
  # } else {
  #   tree$RE <- FALSE
  # }
  # 
  data <- doPCA(data)
  
  ## First rewrite the labels 
  for(c in tree$children){
    
    namesChildren <- as.vector(c$Get('name')) 
    idxChildren <- labels %in% namesChildren
    
    newLabels[idxChildren] <- c$name
    
  }
  
  data$response = newLabels
  
  cat("Cell types in this layer of the tree:", unique(newLabels), "\n", sep = " ")
  ## get informative PCs and train classifier
  data <- getFeatureSpace(data, pvar = "response", reduction = reduction)
  data <- trainModel(data, model = modelName)
  
  cat("Model trained")
  
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
      trainNode(c, dataSubset, pVar, reduction, verbose, modelName, reconstructionError)
      
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
  data <- RunPCA(data, npcs = npcs, verbose = verbose, seed.use = 42)
  
  data
}


findThreshold <- function(data, labels, verbose = FALSE){
  
  default = data@active.assay
  
  # Add a 5fold CV loop to determine the threshold for the CV
  num_folds = 5
  folds = createFolds(labels, k = num_folds, list = FALSE)
  
  for(i in c(1:num_folds)){
    train_idx = which(folds != i)
    test_idx = which(folds == 1)
    
    train = data[,train_idx]
    test = data[,test_idx]
    test = t(as.matrix(test[[default]]))
    
    # Do pca on train 
    train = doPCA(train)
    
    train = t(as.matrix(train[[default]]@scale.data))
    means = colMeans(train)
    
    colVar <- function(x, ...) {
      sqrt(Matrix::colSums((x - means)^2, ...)/(nrow(x) - 1))
    }
    
    stdevs  <- colVar(train)
    loadings <- Loadings(train, 'pca')

    # First apply same scaling to the test data
    test <- scale(test, means, stdevs)
    
    # Transform the data
    # Order the features of the loadings and of the test data
    
    # Inverse transform the data
    
    
    
  }
  
  
}

