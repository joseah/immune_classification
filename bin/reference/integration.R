#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Integrates azimuth CITE-seq reference
# author: Jose Alquicira Hernandez
# date: 2021-03-26
# description: Integrates reference using donor as batch variable

#   ____________________________________________________________________________
#   HPC details                                                             ####

# screen -S multimodal
# qrsh -N multimodal -l mem_requested=400G
# conda activate r-4.0.3

#   ____________________________________________________________________________
#   Import libraries                                                        ####

# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("Seurat")
library("SeuratDisk")

#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-04-12", "azimuth_reference_integration")

#   ____________________________________________________________________________
#   Import data                                                             ####

data <- LoadH5Seurat(here("data", "seurat_v4", "pbmc_multimodal.h5seurat"))

p <- DimPlot(data, reduction = "umap", group.by = "celltype.l2")
ggsave(here(output, "umap_orig.png"), p)

p <- DimPlot(data, reduction = "umap", group.by = "donor")
ggsave(here(output, "umap_donor_orig.png"), p)


#   ____________________________________________________________________________
#   Process data                                                            ####

data[["RNA"]] <- CreateAssayObject(GetAssayData(data, "counts"))
DefaultAssay(data) <- "RNA"

data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data)
data <- RunUMAP(data, dims = 1:30)


p <- DimPlot(data, reduction = "umap", group.by = "celltype.l2")
ggsave(here(output, "umap.png"), p, width = 7, height = 4)


p <- DimPlot(data, reduction = "umap", group.by = "donor", raster = FALSE)
ggsave(here(output, "umap_donor.png"), p, width = 7, height = 5)


p <- DimPlot(data, reduction = "umap", group.by = "lane", raster = FALSE)
ggsave(here(output, "umap_lane.png"), p, width = 7, height = 5)


p <- DimPlot(data, reduction = "umap", group.by = "time", raster = FALSE)
ggsave(here(output, "umap_time.png"), p, width = 7, height = 5)


#   ____________________________________________________________________________
#   Integrate data                                                          ####

data <- SplitObject(data, "donor")

for(i in seq_along(data)){
  message("Processing batch ", i, "/", length(data))
  data[[i]] <- NormalizeData(data[[i]], verbose = FALSE)
  data[[i]] <- FindVariableFeatures(data[[i]], 
                                    selection.method = "vst", 
                                    nfeatures = 2000, 
                                    verbose = FALSE)
}


inicio("Finding anchors")
anchors <- FindIntegrationAnchors(data, dims = 1:30)
fin()

inicio("Integrating data")
integrated <- IntegrateData(anchors, dims = 1:30)
fin()

#   ____________________________________________________________________________
#   Process integrated data                                                 ####


integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30)


p <- DimPlot(integrated, reduction = "umap", group.by = "celltype.l2", 
             label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
ggsave(here(output, "umap_integrated.png"), p, width = 4.5, height = 4)


p <- DimPlot(integrated, reduction = "umap", group.by = "donor", raster = FALSE)
ggsave(here(output, "umap_integrated_donor.png"), p, width = 7, height = 5)


p <- DimPlot(integrated, reduction = "umap", group.by = "lane", raster = FALSE)
ggsave(here(output, "umap_integrated_lane.png"), p, width = 7, height = 5)


p <- DimPlot(integrated, reduction = "umap", group.by = "time", raster = FALSE)
ggsave(here(output, "umap_integrated_time.png"), p, width = 7, height = 5)


#   ____________________________________________________________________________
#   Export data                                                             ####

saveRDS(integrated, here(output, "integrated.RDS"))

#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))

