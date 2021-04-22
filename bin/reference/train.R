#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Train hierarchical scPred models from CITE-seq data
# author: Lieke Michelsen, Jose Alquicira-Hernandez
# date: 2021-04-14
# description: None

#   ____________________________________________________________________________
#   Import libraries                                                        ####

# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("Seurat")
library("HierscPred")
library("doParallel")

#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-04-22", "hierchical_training")

#   ____________________________________________________________________________
#   Import data                                                             ####


inicio("Read reference")
ref <- readRDS(here("results",
                     "2021-04-12_azimuth_reference_integration",
                     "integrated.RDS"))
fin()

h <- c("pbmc/nonNKT/DC/ncDC/pDC",
       "pbmc/nonNKT/DC/ncDC/ASDC",
       "pbmc/nonNKT/DC/cDC/cDC1",
       "pbmc/nonNKT/DC/cDC/cDC2",
       "pbmc/nonNKT/HSPC-Platelet-Eryth/HSPC",
       "pbmc/nonNKT/HSPC-Platelet-Eryth/Platelet",
       "pbmc/nonNKT/HSPC-Platelet-Eryth/Eryth",
       "pbmc/nonNKT/B cell lineage/Plasmablast",
       "pbmc/nonNKT/B cell lineage/B cell/B memory",
       "pbmc/nonNKT/B cell lineage/B cell/B naive",
       "pbmc/nonNKT/B cell lineage/B cell/B intermediate",
       "pbmc/nonNKT/Monocyte/Doublet",
       "pbmc/nonNKT/Monocyte/CD14 Mono",
       "pbmc/nonNKT/Monocyte/CD16 Mono",
       "pbmc/NKT/Cycling cell/CD4 Proliferating",
       "pbmc/NKT/Cycling cell/CD8 Proliferating",
       "pbmc/NKT/Cycling cell/NK Proliferating",
       "pbmc/NKT/NK cell-ILC/ILC",
       "pbmc/NKT/NK cell-ILC/NK cell/NK",
       "pbmc/NKT/NK cell-ILC/NK cell/NK_CD56bright",
       "pbmc/NKT/CD4 CTL-CD8 TEM-gdT-MAIT/MAIT",
       "pbmc/NKT/CD4 CTL-CD8 TEM-gdT-MAIT/CD4 CTL-CD8 TEM-gdT/gdT",
       "pbmc/NKT/CD4 CTL-CD8 TEM-gdT-MAIT/CD4 CTL-CD8 TEM-gdT/CD4 CTL",
       "pbmc/NKT/CD4 CTL-CD8 TEM-gdT-MAIT/CD4 CTL-CD8 TEM-gdT/CD8 TEM",
       "pbmc/NKT/T Naive-CM-EM/CD4 TEM-CD4 TCM-CD8 TCM/CD8 TCM",
       "pbmc/NKT/T Naive-CM-EM/CD4 TEM-CD4 TCM-CD8 TCM/CD4 TEM",
       "pbmc/NKT/T Naive-CM-EM/CD4 TEM-CD4 TCM-CD8 TCM/CD4 TCM",
       "pbmc/NKT/T Naive-CM-EM/Naive-Treg-dnT/dnT",
       "pbmc/NKT/T Naive-CM-EM/Naive-Treg-dnT/CD8 Naive",
       "pbmc/NKT/T Naive-CM-EM/Naive-Treg-dnT/Treg",
       "pbmc/NKT/T Naive-CM-EM/Naive-Treg-dnT/CD4 Naive")


hier <- create_hierarchy(h)

cl <- makePSOCKcluster(6)
registerDoParallel(cl)

hier <- train_tree(ref, 
                   hier, 
                   pvar = 'celltype.l2', 
                   model = "mda", 
                   reconstruction_error = FALSE)
stopCluster(cl)



#   ____________________________________________________________________________
#   Export data                                                             ####

save(hier, ref, file = 'trained_ref_noRE.RData')

#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))

