createTree <- function(stringTree){
  
  # Create dataframe
  h <- tibble(pathString = stringTree)
  
  # Convert dataframe to a tree
  h <- as.Node(h)
  
  h
}

# tree = tree we want to train, we construct this tree before we call this function?
# data = seuratobject of the reference to use for training
# pVar = name of the column of the metadata used as labels

trainTree <- function(tree, data, pVar = 'cell_type', reduction = 'pca', verbose = FALSE, modelName = 'svmRadial',reconstructionError = TRUE,FN_perc = 0.01){
  
  labels <- data@meta.data[[pVar]]
  
  tree <-  trainNode(tree, data, pVar, reduction, verbose, modelName, reconstructionError, FN_perc)
  tree
}


# Iterate over the nodes recursively
trainNode <- function(tree, data, pVar, reduction, verbose, modelName, reconstructionError, FN_perc){
  
  cat("Training parent node: ", tree$name, "\n", sep = "")
  
  labels <- data@meta.data[[pVar]]
  newLabels <- labels
  
  if(reconstructionError){
    tree$RE <- findThreshold(data, labels, FN_perc, verbose)
  } else {
    tree$RE <- FALSE
  }

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
      trainNode(c, dataSubset, pVar, reduction, verbose, modelName, reconstructionError, FN_perc)
      
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


findThreshold <- function(data, labels, FN_perc, verbose = FALSE){
  
  default = DefaultAssay(data)
  
  # Add a 5fold CV loop to determine the threshold for the CV
  num_folds = 5
  folds = createFolds(labels, k = num_folds, list = FALSE)
  RE_threshold <- 0
  
  for(i in c(1:num_folds)){
    train_idx = which(folds != i)
    test_idx = which(folds == i)
    
    train = data[,train_idx]
    test = data[,test_idx]

    # PCA on train, get feature loadings 
    train = doPCA(train)
    
    loadings <- Loadings(train, 'pca')
    features <- rownames(loadings)
    
    train <- as.matrix(GetAssayData(train, "data", assay = default)[features,])
    
    means = rowMeans(train)
    
    rowVar <- function(x, ...) {
      sqrt(Matrix::rowSums((x - means)^2, ...)/(ncol(x) - 1))
    }
    
    stdevs  <- rowVar(train)

    # Transform test data and calculate RE
    test <- t(as.matrix(GetAssayData(test, "data", assay = default)[features,]))
    test <- scale(test, means, stdevs)
    
    test_transformed <- test %*% loadings
    test_inverse <- test_transformed %*% t(loadings)
    
    RE = test_inverse - test
    RE = RE ^ 2
    RE = rowSums(RE)
    RE = sqrt(RE)
    
    RE_threshold[i] <- quantile(RE, 1-FN_perc)
  }
  
  RE_threshold <- median(unlist(RE_threshold))
  
  RE_threshold
  
}

