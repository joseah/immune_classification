predictTree <- function(tree, newData, threshold = 0.75, recompute_alignment = TRUE){
  
  newData <- predictNode(tree, newData, threshold, recompute_alignment)
  
  newData
  
}

# Iterate over the nodes recursively
predictNode <- function(tree, newData, threshold, recompute_alignment){
  
  if(tree$RE){
    # Look if cells are rejected based on RE
    print(cat('Threshold ',tree$RE, sep = ''))
    
    RE_rejected = reconstructionError(newData, tree$model, tree$RE)
    
    RE_rejected_idx = which(RE_rejected)
    print(cat("Rejected by RE: ", length(RE_rejected_idx), sep = ''))
    RE_notrejected_idx = which(!RE_rejected)
    print(cat("Not rejected by RE: ",length(RE_notrejected_idx), sep = ''))
    
    # Make predictions (only for cells that are not rejected!)
    newData1 <- newData[,RE_notrejected_idx]
    newData1 <- scPredict(newData1, tree$model, threshold = threshold, recompute_alignment = recompute_alignment)
    
    newData$scpred_prediction = ''
    newData$scpred_prediction[RE_notrejected_idx] = newData1$scpred_prediction
    newData$scpred_prediction[RE_rejected_idx] = 'unassigned'
    
  } else{
    
    newData <- scPredict(newData, tree$model, threshold = threshold, recompute_alignment = recompute_alignment)
    
  }
  
  
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
      dataSubset <- predictNode(c, dataSubset, threshold, recompute_alignment)
      
      newData$scpred_prediction[idxChildren] <- dataSubset$scpred_prediction
      
    }
  }
  
  newData    
  
}

reconstructionError <- function(newData, spmodel, threshold){
  
  ref_loadings <- spmodel@feature_loadings
  new_features <- rownames(newData)
  reference_features <- rownames(ref_loadings)
  
  shared_features <- intersect(reference_features, new_features)
  
  ref_loadings <- ref_loadings[shared_features,]
  
  new_data <- GetAssayData(newData, "data")[shared_features,]
  means <- spmodel@scaling$means
  stdevs  <- spmodel@scaling$stdevs
  new_data <- Matrix::t(new_data)
  
  scaled_data <- scale(new_data, means, stdevs)
  new_embeddings <- scaled_data %*% ref_loadings
  
  new_inverse = new_embeddings %*% t(ref_loadings)
  
  new_inverse = new_inverse[, order(colnames(new_inverse))]

  new_original = as.matrix(scaled_data[, order(colnames(scaled_data))])

  RE = new_inverse - new_original
  RE = RE ^ 2
  RE = rowSums(RE)
  RE = sqrt(RE)
  
  return(RE > threshold)
  
}

