## Implement metrics

# Percentage rejected
# Percentage internal
# HF1-score per cell type

evaluate <- function(y_true, y_pred, tree){
  
  # %rejected
  root = tree$name
  rejected = 100*sum(y_pred == root)/length(y_pred)
  
  # % internal
  leafs = data.frame(tree$leaves)
  internal = 100*(length(y_pred) - sum(y_pred %in% leafs) - sum(y_pred == root))/length(y_pred)
  
  # Hierarchical F1-score per cell type
  HF1_all = c()
  HF1_names <- unique(y_true)
  
  for(ct in unique(y_true)){
    
    sum_pred = 0
    sum_true = 0
    sum_overlap = 0
    
    set_true = findSet(ct, tree)
    
    if(is.null(set_true)){
      HF1_names <- HF1_names[HF1_names != ct]
      next
    }
    
    #iterate of all the cells of this cell type
    idx = which(y_true == ct)
    for(i in idx){
      
      pred = y_pred[i]
      
      if(pred == tree$name){
        set_pred = c()
      } else{
        set_pred = findSet(pred, tree)
      }
      
      overlap = length(intersect(set_true, set_pred))
      
      if((length(set_pred) > length(set_true)) & (length(set_true == overlap))){
        sum_pred = sum_pred + overlap
      } else{
        sum_pred = sum_pred + length(set_pred)
      }
      
      sum_true = sum_true + length(set_true)
      sum_overlap = sum_overlap + overlap
      
    }
    
    HP <- sum_overlap / sum_pred
    HR <- sum_overlap / sum_true
    
    HF1 <- (2*HP*HR)/(HP+HR)
    HF1_all <- append(HF1_all, HF1)
    
  }
  
  names(HF1_all) <- HF1_names
  HF1 <- median(HF1_all)
  
  return(list(HF1 = HF1, rejected = rejected, internal = internal, HF1_all = HF1_all))
  
}

findSet <- function(name_ct, tree){
  
  # print(name_ct)
  nodes_set <- c()
  nodes_all <- Traverse(tree)
  
  for(no in nodes_all){
    
    # Correct node found, now look at ancestors
    if(no$name == name_ct){
      nodes_set <- append(nodes_set, name_ct)
      while(no$parent$name != tree$name){

        nodes_set = append(nodes_set, no$parent$name)
        no <- no$parent
      }
    }
  }
  
  nodes_set
  
}
