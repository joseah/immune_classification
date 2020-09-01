predictTree <- function(tree, newData, recompute_alignment = TRUE){
  
  newData <- predictNode(tree, newData, recompute_alignment)
  
  newData
  
}

# Iterate over the nodes recursively
predictNode <- function(tree, newData, recompute_alignment){
  
  # Make predictions  
  newData <- scPredict(newData, tree$model, recompute_alignment = recompute_alignment)
  
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
      dataSubset <- predictNode(c, dataSubset, recompute_alignment)
      
      newData$scpred_prediction[idxChildren] <- dataSubset$scpred_prediction
      
    }
  }
  
  newData    
  
}
