# Script information ------------------------------------------------------

# title: 
# author: Jose Alquicira Hernandez
# date: 2020-04-23
# description: None

# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("Seurat")
library("scPred")
library('data.tree')


# Set output --------------------------------------------------------------

output <- set_output("2020-04-23", "")

# Read data ---------------------------------------------------------------

filename <- here("results", "2020-04-14_68k_annotation", "68k_annotated.RDS")
pbmc_68k <- readRDS(filename)


filename <- here("results", "2020-05-01_benchmark_pbmc_annotation_v2", "v2.RDS")
v2 <- readRDS(filename)


set.seed(66)
pbmc_68k <- pbmc_68k[, sample(1:ncol(pbmc_68k), 6000)]

pbmc_68k <- pbmc_68k[, pbmc_68k$cell_type != "Plasma"]

pbmc_68k$cell_type <- fct_drop(pbmc_68k$cell_type)

pbmc_68k$cell_type %>% unique()

# Training ----------------------------------------------------------------

hier <- c("pbmc/Myeloid/DC/mDCs",
          "pbmc/Myeloid/DC/pDC",
          "pbmc/Myeloid/Monocyte/CD16+ Monocytes",
          "pbmc/Myeloid/Monocyte/CD14+ Monocytes",
          "pbmc/Lymphoid/B cells",
          "pbmc/Lymphoid/NKT cell/T cell/CD8+ T cell",
          "pbmc/Lymphoid/NKT cell/T cell/CD4+ T cell",
          "pbmc/Lymphoid/NKT cell/Natural Killer")


cell_type <- hier %>% str_split("/") %>% map(~ .[length(.)]) %>% flatten_chr()
h <- tibble(pathString = hier, cell_type)
h <- as.Node(h)
h_test <- Clone(h)


levels <- h$Get(function(x) names(x$children))
levels <- names(levels[!is.na(levels)])
nodes <- lapply(levels, function(x) FindNode(h, x))
trainNode <- function(node, data, pVar){
  
  node$children %>% 
    map(~ .$Get("cell_type")) %>% 
    map(~ .[!is.na(.)]) %>% 
    map(as.character) -> parent_leaves
  
  
  parents <- names(parent_leaves)
  labels <- data@meta.data[, pVar] %>% as.character()
  lineages <- as.list(parent_leaves)
  
  assign_parent <- function(cell, lineages){
    lapply(parent_leaves, function(l){cell %in% l}) %>% 
      unlist() %>% 
      which() %>% 
      {if(length(.)) names(.) else NA}
  }
  
  target_labels <- data@meta.data[, pVar]
  data$scPred_pVar <- sapply(target_labels, assign_parent)
  
  data <- data[,!is.na(data$scPred_pVar)]
  if(ncol(data) < 100){
    npcs <- ncol(data) - 2
  }else{
    npcs <- 100
  }
  
  data <- FindVariableFeatures(data)
  data <- ScaleData(data)
  data <- RunPCA(data, npcs = npcs)
  data <- getFeatureSpace(data, "scPred_pVar")
  data <- trainModel(data)
  data
}
#undebug(trainNode)
#debug(trainNode)
models <- lapply(nodes, function(x) trainNode(x, pbmc_68k, pVar = "cell_type"))


names(models) <- levels
nodes_order <- h$Get(function(x) x$name) %>% as.character()
i <- if_else(nodes_order %in% names(models), TRUE, FALSE)
classifier <- as.list(i)
classifier[i] <- models
h$Set(classifier)
h
# Classification ----------------------------------------------------------
tree <- Clone(h)
tree$cells <- Cells(v2)
new_data <- v2
layers <- seq_len(tree$height)
layers <- layers[-length(layers)]

classify_children <- function(node){
  new_data <- new_data[, node$cells]
  new_data <- scPredict(node$classifier, new_data, tolerance = 0.3)
  pred <- split(new_data$scPred_prediction, new_data$scPred_prediction)
  children_names <- names(node$children)
  
  node$Do(function(x){
    j <- children_names %in% x$name
    if(any(j)){
      x$cells <- names(pred[[x$name]])
    }
  })
  node
}

for(i in layers){
  tree$Do(classify_children, 
          filterFun = function(x) x$level == layers[i] & is(x$classifier, "Seurat"))
}



pred <- tree$Get(function(x) cells = x$cells, filterFun = isLeaf, simplify = FALSE)
imap(pred, ~ data.frame(cell_type = rep(.y, length(.x)), row.names = .x)) %>% 
  Reduce(rbind, .) -> pred


v2 <- AddMetaData(v2, metadata = pred, col.name = "scPred_class")
v2$scPred_class <- if_else(is.na(v2$scPred_class), "unassigned", as.character(v2$scPred_class))

DimPlot(v2, group.by = "scPred_class", label = TRUE)


crossTab(v2, true = "cell_type", pred ="scPred_class", output = "prop")


pbmc_68k <- getFeatureSpace(pbmc_68k, "cell_type", removeOutliers = FALSE)
pbmc_68k <- trainModel(pbmc_68k)
facs <- scPredict(pbmc_68k, facs)

crossTab(facs, "cell_type", "scPred_prediction", output = "prop")

# Save data ---------------------------------------------------------------


# Session info ------------------------------------------------------------

print_session(here(output))
