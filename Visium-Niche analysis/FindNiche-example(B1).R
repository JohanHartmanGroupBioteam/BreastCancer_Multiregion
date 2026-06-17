#### loading PKGs ####

library(Seurat)
library(future)
library(hdf5r)
library(arrow)
# plan("multisession", workers = 10)

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(spacexr)
library(circlize)
library(ComplexHeatmap)

library(scatterpie)
library(jsonlite)

library(SPOTlight)

options(future.globals.maxSize = 50 * 1024^3)

#### paths ####
visiumpath="/Users/qiao.yang/OneDrive - Karolinska Institutet/Karolinska Ins/ProjectsAtKI/1.MultiRegin/8.visuim/visium_spatial_v4"
respath="/Users/qiao.yang/OneDrive - Karolinska Institutet/Karolinska Ins/ProjectsAtKI/1.MultiRegin/8.visuim/visium_spatial_res"
CATpath="/Users/qiao.yang/OneDrive - Karolinska Institutet/Karolinska Ins/ProjectsAtKI/1.MultiRegin/8.visuim/CTA_v4"
GTFfile="/Users/qiao.yang/OneDrive - Karolinska Institutet/Karolinska Ins/ProjectsAtKI/1.MultiRegin/1.WGS/Datasource"
refpath="/Users/qiao.yang/OneDrive - Karolinska Institutet/Karolinska Ins/ProjectsAtKI/1.MultiRegin/7.xenium/data/Wu_etal_2021_BRCA_scRNASeq"

foldername1 = "V19T26-012_A1/outs"
foldername2 = "V19T26-012_D1/outs"
patient = "BCSA1"


## loading processed obj 
# save(BCSA1.merge, file = paste(respath, patient,"Visium_processed.RData", sep = "/" ) )
load(file = paste(respath, patient,"Visium_processed.RData", sep = "/" ))

#### Find niche ####
library(nicheDE)
source("CreateNicheDEObjectFromSeurat_updated.R")

meta = BCSA1.merge[[]]

cell.pct = meta %>%
  # dplyr::filter(!is.na(celltype)) %>%
  dplyr::select(7:15)

#### Average Expression Profile From a scRNA-seq
Idents(bcsa.data) <- "celltype_major"
vignette_library_matrix  = CreateLibraryMatrixFromSeurat(bcsa.data,assay = 'RNA')
rownames(vignette_library_matrix) = str_replace( rownames(vignette_library_matrix),"-| ", ".")

#### create obj
## B1 
cells = rownames(meta)[grepl("B1", rownames(meta))]
BCSA1.merge.B1 = subset(BCSA1.merge, cells = cells )
vignette_deconv_mat = cell.pct[cells, rownames(vignette_library_matrix)]
coordinates= GetTissueCoordinates(BCSA1.merge.B1)

NDE_obj1 = CreateNicheDEObjectFromSeurat( BCSA1.merge.B1, 'SCT',coordinates = coordinates,
                                          vignette_library_matrix, as.matrix( vignette_deconv_mat),sigma = c(1,100,250))

## calculation
NDE_obj1 = CalculateEffectiveNiche(NDE_obj1)
print(NDE_obj1)

## A1 
cells = rownames(meta)[grepl("A1", rownames(meta))]
BCSA1.merge.filt = subset(BCSA1.merge, cells = cells )
vignette_deconv_mat = cell.pct[cells, rownames(vignette_library_matrix)]
coordinates= GetTissueCoordinates(BCSA1.merge.filt)

NDE_obj2 = CreateNicheDEObjectFromSeurat( BCSA1.merge.filt, 'SCT', coordinates = coordinates,
                                          vignette_library_matrix, as.matrix( vignette_deconv_mat), sigma = c(1,100,250))

## calculation
NDE_obj2 = CalculateEffectiveNiche(NDE_obj2)
print(NDE_obj2)

## B1
effective = as.data.frame( NDE_obj1@effective_niche) 

# Perform clustering
set.seed(123)
kmeans_result <- kmeans(effective, centers = 5)  # Start with 5 clusters

# Add cluster labels to your data
niche_clusters <- data.frame(visium = names(kmeans_result$cluster), Niches = kmeans_result$cluster )
niche_clusters[,1] = NULL

# If you have a Seurat object:
BCSA1.merge = AddMetaData(BCSA1.merge, niche_clusters)


## A1
effective = as.data.frame( NDE_obj2@effective_niche) 

# Perform clustering
set.seed(123)
kmeans_result <- kmeans(effective, centers = 5)  # Start with 5 clusters

# Add cluster labels to your data
niche_clusters <- data.frame(visium = names(kmeans_result$cluster), Niches = kmeans_result$cluster )
niche_clusters[,1] = NULL

# If you have a Seurat object:
BCSA1.merge = AddMetaData(BCSA1.merge, niche_clusters)

#### visualization ####

## vis
p = SpatialDimPlot(BCSA1.merge, group.by = "Niches",images = "B1" )+ ggtitle("") +
  # scale_fill_manual(values =  niches_colrs, na.value = "gray90",
  #                   name = "Niches",
  #                   limits = function(x) {
  #                     # Exclude NA from the limits (and thus from legend)
  #                     unique(x)[!is.na(unique(x))]
  #                   }
  # ) +
  guides(fill = guide_legend(
    override.aes = list(size = 4, alpha = 1)  # Increase size and set alpha to 1
  )) +
  theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.5)
  )
p
ggsave( paste(respath, patient,"visium_BCSA1_Niches_imgB1.png", sep = "/" ), p, width = 9, height = 7, dpi = 600, bg = "white")


p = SpatialDimPlot(BCSA1.merge, group.by = "Niches",images = "A1") +
  ggtitle("") +
  # scale_fill_manual(values =  niches_colrs, na.value = "gray90",
  #                   name = "Niches",
  #                   limits = function(x) {
  #                     # Exclude NA from the limits (and thus from legend)
  #                     unique(x)[!is.na(unique(x))]
  #                   }
  # ) +
  guides(fill = guide_legend(
    override.aes = list(size = 4, alpha = 1)  # Increase size and set alpha to 1
  )) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.5)
  )
p
ggsave( paste(respath, patient,"visium_BCSA1_Niches_imgA1.png", sep = "/" ), p, width = 9, height = 7, dpi = 600, bg = "white")


