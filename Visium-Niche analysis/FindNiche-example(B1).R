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

library("fastcluster")
library(dendextend)

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

chrs = paste0("chr", c(1:22, "X", "Y"))
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

source("scale_infercnv.R")

## loading processed obj 
# save(BCSA1.merge, file = paste(respath, patient,"Visium_processed.RData", sep = "/" ) )
load(file = paste(respath, patient,"Visium_processed.RData", sep = "/" ))

#### info ####
## visualization
cnv_colrs = setNames( c("#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", "#8c564b", "#e377c2", "#b5bd61"),paste0("C", c(1:8) ) )

celltype_colors = setNames(  c("#442288",  "#FED23F", brewer.pal(7, "Paired")), 
                             c("CAFs","Cancer Epithelial", "T-cells", 'Myeloid', "B-cells", "PVL", "Normal Epithelial", "Plasmablasts", "Endothelial" )
)

niches_colrs = setNames( brewer.pal(8, "Dark2") , c(1:8))
nicheHL_colrs = setNames( c("#e31a1c","#fd8d3c","#fecc5c","#ffffb2" ), c("HH", "HL", "LH", "LL"))

##### load ref for RCTD ####

# meta = read.table(file = paste( refpath,"metadata.csv", sep = "/"), sep = ",", quote = "", header = TRUE ) %>%
#   data.frame() %>%
#   magrittr::set_rownames(.$X) %>%
#   dplyr::select(-X)
# 
# bcsa.data <- Read10X(data.dir = refpath,  gene.column = 1)
# bcsa.data <- CreateSeuratObject(counts = bcsa.data, project = "TNBC", min.cells = 3, min.features = 200)
# bcsa.data@meta.data = meta
# 
# 
# bcsa.data = subset(bcsa.data, subtype == "TNBC" )
# dim(bcsa.data)
# 
# save(bcsa.data, file = paste(refpath, "bcsa.data_TNBC.RData", sep = "/"))
load(file = paste(refpath, "bcsa.data_TNBC.RData", sep = "/"))


##### load grch38 ####

# ## grch38 gene loci
load(file = paste(GTFfile, "grch38_gencode.v46.annotation.gtf.RData",sep = "/" ))
grch38_genes = grch38 %>%
  dplyr::filter( type == "gene") %>%
  dplyr::select(seqnames, start , end, width, strand, gene_id, gene_type, gene_name) %>%
  dplyr::rename( start.g = start, end.g = end  )
# rm(grch38)


#### Visium QC and RCTD: A1 (obj1) ####

##### Load #####

## obj1
visium.obj1 <- Load10X_Spatial(data.dir =  paste(visiumpath, foldername1, sep = "/"),   
                               assay = "Spatial", slice = "A1")
dim(visium.obj1)
visium.obj1[["percent.mt"]] <- PercentageFeatureSet(visium.obj1, pattern = "^MT-")

## quality control
meta <- visium.obj1@meta.data
ggplot(meta, aes( x = log2(nFeature_Spatial), log2(nCount_Spatial))) +
  geom_vline(xintercept = 9, color = "red")+
  geom_hline(yintercept = 9, color = "red")+
  geom_point(alpha = 0.5, size = 0.5)+
  theme_classic()

## subsetting
visium.obj1 <- subset(visium.obj1, subset = log2( nCount_Spatial +1) > 9 & log2(nFeature_Spatial +1) > 9 & percent.mt < 20 )
dim(visium.obj1)

## SCT
visium.obj1 <- SCTransform(visium.obj1, assay = "Spatial")


##### RCTD #####
## query data
print("query data")
query.counts <- GetAssayData(visium.obj1, assay = "Spatial", layer = "counts")
coords <- GetTissueCoordinates(visium.obj1)
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

## ref data
print("ref data")
Idents(bcsa.data) <- "celltype_major"
counts <- GetAssayData(bcsa.data, assay = "RNA", layer = "counts")
cluster <- as.factor(bcsa.data$celltype_major)
names(cluster) <- colnames(bcsa.data)
nUMI <- bcsa.data$nCount_RNA
names(nUMI) <- colnames(bcsa.data)
nUMI <- colSums(counts)
reference <- Reference(counts, cluster, nUMI, n_max_cells = 20000)

## run RCTD with many cores
print("run RCTD")
RCTD <- create.RCTD(query, reference, max_cores = 6)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

## save result
weights = RCTD@results$weights
norm_weights = normalize_weights(weights = weights)

norm_weights = data.frame(norm_weights, check.names = FALSE)
visium.obj1 <- AddMetaData(visium.obj1, metadata = norm_weights)

#### Visium QC and RCTD: B1 (obj2) ####

##### Load #####

## obj2
visium.obj2 <- Load10X_Spatial(data.dir =  paste(visiumpath, foldername2, sep = "/"),   
                               assay = "Spatial", slice = "B1")
dim(visium.obj2)
visium.obj2[["percent.mt"]] <- PercentageFeatureSet(visium.obj2, pattern = "^MT-")

## quality control
meta <- visium.obj2@meta.data
ggplot(meta, aes( x = log2(nFeature_Spatial), log2(nCount_Spatial))) +
  geom_vline(xintercept = 9, color = "red")+
  geom_hline(yintercept = 9, color = "red")+
  geom_point(alpha = 0.5, size = 0.5)+
  theme_classic()

## subsetting
visium.obj2 <- subset(visium.obj2, subset = log2( nCount_Spatial +1) > 9 & log2(nFeature_Spatial +1) > 9 & percent.mt < 20 )
dim(visium.obj2)


## SCT
visium.obj2 <- SCTransform(visium.obj2, assay = "Spatial")


##### Run RCTD #####

## query data
print("query data")
query.counts <- GetAssayData(visium.obj2, assay = "Spatial", layer = "counts")
coords <- GetTissueCoordinates(visium.obj2)
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

## ref data
print("ref data")
Idents(bcsa.data) <- "celltype_major"
counts <- GetAssayData(bcsa.data, assay = "RNA", layer = "counts")
cluster <- as.factor(bcsa.data$celltype_major)
names(cluster) <- colnames(bcsa.data)
nUMI <- bcsa.data$nCount_RNA
names(nUMI) <- colnames(bcsa.data)
nUMI <- colSums(counts)
reference <- Reference(counts, cluster, nUMI, n_max_cells = 20000)

## run RCTD with many cores
print("run RCTD")
RCTD <- create.RCTD(query, reference, max_cores = 6)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

## save result
weights = RCTD@results$weights
norm_weights = normalize_weights(weights = weights)

norm_weights = data.frame(norm_weights, check.names = FALSE)
visium.obj2 <- AddMetaData(visium.obj2, metadata = norm_weights)

#### run integrate ####

BCSA1.merge <- merge(visium.obj1, visium.obj2, add.cell.ids = c("A1", "B1"))

DefaultAssay(BCSA1.merge) <- "SCT"

VariableFeatures(BCSA1.merge) <- c(VariableFeatures(visium.obj1), VariableFeatures(visium.obj2))
BCSA1.merge <- RunPCA(BCSA1.merge, verbose = FALSE)
BCSA1.merge <- FindNeighbors(BCSA1.merge, dims = 1:30)
BCSA1.merge <- FindClusters(BCSA1.merge, verbose = FALSE)
BCSA1.merge <- RunUMAP(BCSA1.merge, dims = 1:30)

DimPlot(BCSA1.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
SpatialDimPlot(BCSA1.merge)

#### visualization RCTD ####
celltype_colors = setNames(  c("#442288",  "#FED23F", brewer.pal(7, "Paired")), 
                             c("CAFs","Cancer Epithelial", "T-cells", 'Myeloid', "B-cells", "PVL", "Normal Epithelial", "Plasmablasts", "Endothelial" )
)


## test 
metadata = BCSA1.merge[[]]
colnames(metadata)[7:15] = c("B-cells", "CAFs", "Cancer Epithelial","Endothelial", "Myeloid","Normal Epithelial","Plasmablasts", "PVL", "T-cells")

## A1
spatial_image <- GetImage(BCSA1.merge,image = "A1", mode = "raster")

test_plot <- SpatialDimPlot(BCSA1.merge, images = "A1", crop = FALSE)
plot_limits <- ggplot_build(test_plot)$layout$panel_params[[1]]
gg_data <- ggplot_build(test_plot)

point_data <- gg_data[["plot"]][[1]][["data"]]

# Combine coordinates with cell type percentages
pie_data <- cbind(point_data, metadata[rownames(point_data),] )

A1 = plotSpatialScatterpie(
  x = pie_data[,c("x","y")],
  y = pie_data[, names(celltype_colors)],
  cell_types = names(celltype_colors),
  img = spatial_image,
  pie_scale = 0.3, 
  # Rotate the image 90 degrees counterclockwise
  degrees = -90,
  # Pivot the image on its x axis
  axis = "h") +
  ggtitle("A1") +
  scale_fill_manual(values = celltype_colors, name = "RCTD anno.") + 
  theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.5)
  )
A1

# B1
spatial_image <- GetImage(BCSA1.merge,image = "B1", mode = "raster")

test_plot <- SpatialDimPlot(BCSA1.merge, images = "B1", crop = FALSE)
plot_limits <- ggplot_build(test_plot)$layout$panel_params[[1]]
gg_data <- ggplot_build(test_plot)

point_data <- gg_data[["plot"]][[1]][["data"]]

# Combine coordinates with cell type percentages
pie_data <- cbind(point_data, metadata[rownames(point_data),] )

B1 = plotSpatialScatterpie(
  x = pie_data[,c("x","y")],
  y = pie_data[, names(celltype_colors)],
  cell_types = names(celltype_colors),
  img = spatial_image,
  pie_scale = 0.3, 
  # Rotate the image 90 degrees counterclockwise
  degrees = -90,
  # Pivot the image on its x axis
  axis = "h") +
  ggtitle("B1") +
  scale_fill_manual(values = celltype_colors, name = "RCTD anno.") + 
  theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.5)
  )
B1


# combine plots
combined.plot <- (A1 | B1) + plot_layout(guides = "collect")
combined.plot

ggsave( file = paste(respath, patient, paste0( "visium_", patient, "_RCTD_plots.png"), sep = "/"), combined.plot, width = 9, height = 6, dpi = 600, bg = "white")

#### run infercnv: make obj for Bianca ####

#### GTF
print("loading GTF")
# load(paste(GTFfile, "grch38_gencode.v46.annotation.gtf.RData",sep = "/" ))

genes.grch38 = grch38 %>%
  dplyr::filter( type== "gene" ) %>%
  dplyr::select(seqnames, start, end, gene_name ) 

genes.grch38 = genes.grch38[!duplicated(genes.grch38$gene_name),]
rownames(genes.grch38) = genes.grch38$gene_name
genes.grch38[,4] = NULL
colnames(genes.grch38) = NULL

rm(grch38)

## read CTA3
A1 = read.csv(file = paste(CATpath, "BCSA1TumA1_spaceranger4.0.1_aligned.csv", sep = "/" )) %>%
  mutate(
    cellid = paste0("A1_", cellid)
  )

B1 =  read.csv(file = paste(CATpath, "BCSA1TumB2_spaceranger4.0.1_aligned.csv", sep = "/" ), header = T) %>%
  mutate(
    cellid = paste0("B1_", cellid)
  )

CTA_tumor = rbind(A1, B1) %>%
  dplyr::filter(tumor_per > 0.9 )
dim(CTA_tumor)
# 2060    9

CTA_normal = rbind(A1, B1) %>%
  dplyr::filter(tumor_per == 0 )
dim(CTA_normal)
# 943   9

CTA_all =  rbind(A1, B1[,colnames(A1)])
rownames(CTA_all) = CTA_all$cellid
CTA_all[,c(1,2)] = NULL


## adding to meta
BCSA1.merge = AddMetaData(BCSA1.merge, CTA_all)

#### high normal cell from molecular
meta = BCSA1.merge[[]]

cells_normal = rownames(meta)[meta$Cancer.Epithelial < 0.01 & !is.na(meta$Cancer.Epithelial)]
length(cells_normal)
# 26

cells_tumor = rownames(meta)[meta$Cancer.Epithelial > 0.9 & !is.na(meta$Cancer.Epithelial)]
length(cells_tumor)
# 1126

cells_tumor_all = unique(c(CTA_tumor$cellid, cells_tumor ))
length(cells_tumor_all)
# 2130

cells_normal_all = unique( c(CTA_normal$cellid, cells_normal ))
length(cells_normal_all)
# 953

## excluding ambuig cells
BC_amb = intersect(cells_normal_all, cells_tumor_all)
length(BC_amb)
# 0

cells_normal_all = cells_normal_all[! cells_normal_all %in% BC_amb]
length(cells_normal_all)
# 953

cells_tumor_all = cells_tumor_all[! cells_tumor_all %in% BC_amb]
length(cells_tumor_all)
# 2130

intersect(cells_normal_all, cells_tumor_all) %>% length()

## check tumor, filter
meta = BCSA1.merge[[]]

meta$celltype =  ifelse( rownames(meta) %in% cells_normal_all , "Normal", 
                         ifelse( rownames(meta) %in% cells_tumor_all, "Tumor", NA ) )

anncell = meta %>%
  dplyr::filter(!is.na(celltype)) %>%
  dplyr::select(celltype)
colnames(anncell) = NULL
table(anncell)

# Normal  Tumor 
# 951   2130 

#### run infercnv 

## get counts for each slice and bind by column
print("bind counts")
layers = Layers(BCSA1.merge[["Spatial"]])
layers
layersNames = c("A1", "B1")
names(layersNames) = layers


counts = NULL
for (layer in names(layersNames) ){
  
  counts.temp = SeuratObject::LayerData(object = BCSA1.merge[["Spatial"]], layer = layer )
  counts = cbind(counts, counts.temp)
}

dim(counts)
# 33538  5455

## subsetting
counts = counts[, rownames(anncell) ]
dim(counts)
# 33538  3081

save(counts, anncell, genes.grch38, file =  paste(respath, patient,"infercnvobj.RData", sep = "/" ) )


# create the infercnv object
## on bianca


save( BCSA1.merge, file = paste(respath, patient,"Visium_processed.RData", sep = "/" ))

# ## 
# test = meta %>%
#   dplyr::select(celltype)
# 
# BCSA1.merge = AddMetaData(BCSA1.merge, test)

#### CTA: visualization ####

CTA_all = rbind(A1, B1[,colnames(A1)]) %>%
  magrittr::set_rownames(.$cellid) %>%
  dplyr::select( c("tumor_per", "immune_per","stroma_per" )) 


CTA_colors = setNames(  c("#442288",  "#FED23F", brewer.pal(3, "Paired")[1]), 
                        c("Stroma","Tumor", "Immune" ))


## test 
metadata = BCSA1.merge[[]]

## A1
spatial_image <- GetImage(BCSA1.merge,image = "A1", mode = "raster")

test_plot <- SpatialDimPlot(BCSA1.merge, images = "A1", crop = FALSE)
plot_limits <- ggplot_build(test_plot)$layout$panel_params[[1]]
gg_data <- ggplot_build(test_plot)

point_data <- gg_data[["plot"]][[1]][["data"]]

# Combine coordinates with cell type percentages
pie_data <- cbind(point_data, CTA_all[rownames(point_data),]  )
colnames(pie_data)[(ncol(pie_data)-2):ncol(pie_data)] = c("Tumor", "Immune", "Stroma")

A1 = plotSpatialScatterpie(
  x = pie_data[,c("x","y")],
  y = pie_data[, names(CTA_colors)],
  cell_types = names(CTA_colors),
  img = spatial_image,
  pie_scale = 0.3, 
  # Rotate the image 90 degrees counterclockwise
  degrees = -90,
  # Pivot the image on its x axis
  axis = "h") +
  ggtitle("A1") +
  scale_fill_manual(values = CTA_colors, name = "CTA anno.") + 
  theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.5)
  )
A1


## B1
spatial_image <- GetImage(BCSA1.merge,image = "B1", mode = "raster")

test_plot <- SpatialDimPlot(BCSA1.merge, images = "B1", crop = FALSE)
plot_limits <- ggplot_build(test_plot)$layout$panel_params[[1]]
gg_data <- ggplot_build(test_plot)

point_data <- gg_data[["plot"]][[1]][["data"]]

# Combine coordinates with cell type percentages
pie_data <- cbind(point_data, CTA_all[rownames(point_data),]  )
colnames(pie_data)[(ncol(pie_data)-2):ncol(pie_data)] = c("Tumor", "Immune", "Stroma")

B1 = plotSpatialScatterpie(
  x = pie_data[,c("x","y")],
  y = pie_data[, names(CTA_colors)],
  cell_types = names(CTA_colors),
  img = spatial_image,
  pie_scale = 0.3, 
  # Rotate the image 90 degrees counterclockwise
  degrees = -90,
  # Pivot the image on its x axis
  axis = "h") +
  ggtitle("B1") +
  scale_fill_manual(values = CTA_colors, name = "CTA anno.") + 
  theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.5)
  )
B1


# combine plots
combined.plot <- (A1| B1) + plot_layout(guides = "collect")
combined.plot

ggsave( file = paste(respath, patient, paste0( "visium_", patient, "_CTAanno_plots.png"), sep = "/"), combined.plot, width = 9, height = 6, dpi = 600, bg = "white")




#### After infercnv: visualization ####

infercnv_obj = readRDS(file = paste(paste(respath, patient,"infercnvNoHMM", "run.final.infercnv_obj", sep = "/") ) )

## use scaled data to make plots
infercnv_obj = scale_infercnv(infercnv_obj) ## do scale and clean data

plot_data = infercnv_obj@expr.data %>% t() %>% data.frame()
plot_data = plot_data[ rownames(plot_data) %in% rownames(meta[meta$celltype == "Tumor",]),]

dim(plot_data)
# 2130 8885

#### order genes
# load(paste(GTFfile, "grch38_gencode.v46.annotation.gtf.RData",sep = "/" ))

genes.grch38 = grch38 %>%
  dplyr::filter( seqnames %in% chrs ) %>%
  dplyr::filter( type== "gene") %>%
  dplyr::select(seqnames, start, end, gene_name ) %>%
  mutate( seqnames = factor(seqnames, levels = chrs) ) %>%
  arrange(seqnames,start )

genes.grch38 = genes.grch38[!duplicated(genes.grch38$gene_name),]
rownames(genes.grch38) = factor(genes.grch38$gene_name, levels = genes.grch38$gene_name) 
colnames(genes.grch38) = c("chr","start", "end", "gene" )

genes.grch38.ordered = genes.grch38[genes.grch38$gene %in% colnames(plot_data),] %>%
  dplyr::filter( !grepl("ENSG00",gene ) ) %>%
  mutate(chr = factor(chr, chrs)) %>%
  arrange(chr, start)



# Cap the values
# cnv_capped <- pmin(pmax(cnv, -0.15), 0.2)
quantile(as.matrix( plot_data))

col_fun = colorRamp2(c(min(as.matrix( plot_data)), median(as.matrix( plot_data)), max(as.matrix( plot_data))), c("blue", "white", "red"))


plot_data = plot_data[, colnames(plot_data) %in% genes.grch38.ordered$gene]
dim(plot_data)
# 2130 8825
anncell = meta[ rownames(plot_data) ,] %>%
  mutate(
    Region = str_split(rownames(plot_data) , "_", simplify = T)[,1]
  )%>%
  dplyr::select(Region)
anno_cell = HeatmapAnnotation(df = anncell,which = "row",
                              col = list("Region" =  setNames( brewer.pal( 5, "Paired")[1:2] , c( "A1", "B1")) ))

dend = as.dendrogram(hclust(dist(plot_data), method = "ward.D"))
plot(dend, leaflab = "none")
#
k = 2

p = ComplexHeatmap::Heatmap(as.matrix(plot_data),
                            col = col_fun, name = "Profile",
                            
                            left_annotation = anno_cell,
                            
                            cluster_columns = FALSE,
                            
                            cluster_rows = dend, row_split = k, row_title = "C%s",
                            show_row_dend = FALSE,
                            show_row_names = FALSE,
                            show_column_names = FALSE,
                            use_raster = FALSE,
                            
                            column_split = factor(str_replace(genes.grch38.ordered$chr, "chr", "") ,levels = str_replace(chrs, "chr", "")   ) ,
                            border = TRUE,
                            
                            column_title_gp = gpar(fontsize = 10)
)


p = draw(p)

cairo_pdf(paste(respath, patient,"CNV_visium_BCSA1.pdf", sep = "/" ) , width = 13, height = 4)
p = draw(p)
dev.off()

row.dend <- row_dend(p)  #If needed, extract  dendrogram
row.list <- row_order(p)  #Extract clusters (output is a list)


cluster = data.frame(
  BC = rownames(plot_data),
  clust.cnv = NA
)
for (i in 1:length(row.list)){
  cluster[which(cluster$BC %in% rownames(plot_data)[row.list[[i]]]), "clust.cnv"] =  paste0("C", i)
}


rownames(cluster) = cluster$BC
cluster$BC = NULL

## adding
BCSA1.merge = AddMetaData(BCSA1.merge,cluster )


## visualization
cnv_colrs = setNames( c("#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", "#8c564b", "#e377c2", "#b5bd61"), paste0("C", c(1:8))  )

p = SpatialDimPlot( subset(BCSA1.merge, !is.na(clust.cnv) ) , group.by = "clust.cnv", images = "A1") + ggtitle("") +
  scale_fill_manual(values =  cnv_colrs, na.value = "gray90",
                    name = "Clone",
                    limits = function(x) {
                      # Exclude NA from the limits (and thus from legend)
                      unique(x)[!is.na(unique(x))]
                    }
  ) +
  guides(fill = guide_legend(
    override.aes = list(size = 4, alpha = 1)  # Increase size and set alpha to 1
  )) +
  theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.5)
  )
p
ggsave( paste(respath, patient,"visium_BCSA1_cnv_imgA1.png", sep = "/" ), p, width = 9, height = 7, dpi = 600, bg = "white")

p = SpatialDimPlot(subset(BCSA1.merge, !is.na(clust.cnv) ) , group.by = "clust.cnv", images = "B1") + ggtitle("") +
  scale_fill_manual(values =  cnv_colrs, na.value = "gray90",
                    name = "Clone",
                    limits = function(x) {
                      # Exclude NA from the limits (and thus from legend)
                      unique(x)[!is.na(unique(x))]
                    }
  ) +
  guides(fill = guide_legend(
    override.aes = list(size = 4, alpha = 1)  # Increase size and set alpha to 1
  )) +
  theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.5)
  )
p
ggsave( paste(respath, patient,"visium_BCSA1_cnv_imgB1.png", sep = "/" ), p, width = 9, height = 7, dpi = 600, bg = "white")


# ## save to main meta
# save( BCSA1.merge, file = paste(respath, patient,"Visium_processed.RData", sep = "/" ))


#### BCSA1 AIMS ####
library(BreastSubtypeR)

data("BreastSubtypeRobj")

AIMSgenes = BreastSubtypeRobj$genes.signature$Symbol[BreastSubtypeRobj$genes.signature$AIMS == "Yes"]
AIMSgenes.id = BreastSubtypeRobj$genes.signature$EntrezGene.ID[BreastSubtypeRobj$genes.signature$AIMS == "Yes"]
names(AIMSgenes.id)= AIMSgenes

# extract raw counts
pseudobulk1 <- GetAssayData(BCSA1.merge, assay = "Spatial", layer = "counts.1") %>% data.frame()
pseudobulk2 <- GetAssayData(BCSA1.merge, assay = "Spatial", layer = "counts.2") %>% data.frame()
pseudobulk = cbind(pseudobulk1, pseudobulk2)

## replace genes
setdiff(AIMSgenes,rownames(pseudobulk) ) %>% length()
gens_change = setdiff(AIMSgenes,rownames(pseudobulk) )
gens_change
# "KNTC2" "GARS1" "CYRIB" "MRTFB" "MACIR"
genes_indata = c("NDC80", "GARS", "FAM49B", "MKL2", "C5orf30")

rownames(pseudobulk)[ match(genes_indata, rownames(pseudobulk)) ] = gens_change

## for aims
pseudobulk = pseudobulk[rownames(pseudobulk) %in% AIMSgenes,]
dim(pseudobulk)
#1]  151 5455
colnames(pseudobulk) = str_replace_all(colnames(pseudobulk), "[.]", "-")

## for aims, tumor spot only
pseudobulk = pseudobulk[, Cells(subset(BCSA1.merge, celltype == "Tumor") )]
dim(pseudobulk)
# 151 2130

featuredata = data.frame(gene = rownames(pseudobulk))
featuredata = left_join(featuredata, grch38_genes[, c("gene_name", "width")], by = c("gene" = "gene_name"))
featuredata$width[featuredata$gene== "KNTC2"] = 45079
colnames(featuredata) = c("probe", "Length")

# FPKM from raw counts + gene lengths (in base pairs)
y_all <- edgeR::DGEList(counts = round(as.matrix(pseudobulk)))
gl <- as.numeric(featuredata$Length)
names(gl) <- rownames(pseudobulk)
gl <- gl[rownames(y_all)]
counts.fpkm <- edgeR::rpkm(y_all,
                           gene.length = gl,
                           normalized.lib.sizes = FALSE, log = FALSE
)


res = AIMS::applyAIMS(counts.fpkm, AIMSgenes.id[rownames(counts.fpkm)] )
colnames(res$cl ) = "IS"

BCSA1.merge = AddMetaData(BCSA1.merge,res$cl )

Subtype.color <- c(
  "Basal" = "red",
  "Her2" = "hotpink",
  "LumA" = "darkblue",
  "LumB" = "skyblue",
  "Normal" = "green"
)

## A1
p = SpatialDimPlot( subset(BCSA1.merge, celltype == "Tumor") , images = "A1", group.by = "IS", cols = Subtype.color) +
  scale_fill_manual(values =  Subtype.color, na.value = "gray90",
                    name = "Intrinsic sub."
  ) +
  guides(fill = guide_legend(
    override.aes = list(size = 4, alpha = 1)  # Increase size and set alpha to 1
  )) +
  theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.5)
  )
p
ggsave( paste(respath, patient,"visium_BCSA1_AIMS_imgA1.png", sep = "/" ), p, width = 9, height = 7, dpi = 600, bg = "white")


## B1
p = SpatialDimPlot( subset(BCSA1.merge, celltype == "Tumor") , images = "B1", group.by = "IS", cols = Subtype.color) +
  scale_fill_manual(values =  Subtype.color, na.value = "gray90",
                    name = "Intrinsic sub."
  ) +
  guides(fill = guide_legend(
    override.aes = list(size = 4, alpha = 1)  # Increase size and set alpha to 1
  )) +
  theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.5)
  )
p
ggsave( paste(respath, patient,"visium_BCSA1_AIMS_imgB1.png", sep = "/" ), p, width = 9, height = 7, dpi = 600, bg = "white")


## bar plot 
meta = BCSA1.merge[[]]

plot_data = meta %>%
  dplyr::filter(celltype == "Tumor") %>%
  dplyr::select(clust.cnv, IS)

p = ggplot(plot_data, aes(x = clust.cnv, fill = IS)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values =  Subtype.color, na.value = "gray90",
                    name = "Intrinsic sub.") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "",
    y = "Proportion",
    fill = "IS Type") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(face = "bold"))
p
ggsave( paste(respath, patient,"visium_BCSA1_AIMS_barplot.png", sep = "/" ), p, width = 6, height = 6, dpi = 600, bg = "white")



#### Find niche ####
library(nicheDE)
source("CreateNicheDEObjectFromSeurat.R")

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



#### calcultate immune/stroma
imme.cells = c("B.cells", "T.cells", "Myeloid", "Plasmablasts")
stroma.cells = c("CAFs", "Endothelial", "Normal.Epithelial", "PVL")


meta = BCSA1.merge[[]]

cell.pct = meta %>%
  dplyr::select(7:15, Niches) %>%
  dplyr::mutate(
    sample = str_split(rownames(.), "_", simplify = T)[,1],
    immune.pct = rowSums(.[, imme.cells], na.rm = TRUE),
    stroma.pct = rowSums(.[, stroma.cells], na.rm = TRUE),
    total.pct = immune.pct + stroma.pct
  ) %>%
  group_by(sample, Niches) %>%
  summarise(
    Immune.pct = mean(immune.pct, na.rm = TRUE),  # Average of percentages
    Stroma.pct = mean(stroma.pct, na.rm = TRUE)   # Average of percentages
  ) %>%
  ungroup()

write.table(cell.pct, file = paste(respath, patient,"sections_pct.txt", sep = "/" ), row.names = F, sep = "\t")

# save(BCSA1.merge, file = paste(respath, patient,"Visium_processed.RData", sep = "/" ))




#### Visualization niche celltype ####
imme.cells = c("B.cells", "T.cells", "Myeloid", "Plasmablasts")
stroma.cells = c("CAFs", "Endothelial", "Normal.Epithelial", "PVL")

meta = BCSA1.merge[[]]

meta.plot = meta %>%
  dplyr::mutate(
    sample = str_split(rownames(.), "_", simplify = T)[,1],
    sample = str_sub(sample,1,1)
  ) %>%
  dplyr::select(7:15, sample, Niches) %>%
  pivot_longer( cols = -c("sample", "Niches"), values_to = "pct", names_to = "celltype" ) %>%
  group_by(sample, Niches,celltype) %>%
  summarise(
    pct = mean(pct,na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    celltype = str_replace_all(celltype, "[.]", "-"),
    celltype = str_replace_all(celltype, "Normal-Epithelial", "Normal Epithelial"),
    celltype = str_replace_all(celltype, "Cancer-Epithelial", "Cancer Epithelial"),
    celltype = factor(celltype, levels = c("B-cells","Plasmablasts","T-cells", "Myeloid", "CAFs", "Endothelial", "PVL", "Normal Epithelial", "Cancer Epithelial"))
  ) 



# Basic stacked bar plot
p = ggplot(meta.plot, aes(x = factor(Niches), y = pct, fill = celltype)) +
  geom_col() +
  facet_wrap(~ sample, nrow = 1) +
  scale_fill_manual( values = celltype_colors, name = "Cell type" ) +
  labs(x = "", y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format())+
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    strip.background = element_blank()    
  )
p
ggsave( paste(respath, patient,"B1_niches_CT_barplot.png", sep = "/" ), p, width = 6, height = 6, dpi = 600, bg = "white")

## CTA
meta = BCSA1.merge[[]]

meta.plot = meta %>%
  dplyr::mutate(
    sample = str_split(rownames(.), "_", simplify = T)[,1],
    sample = str_sub(sample,1,1)
  ) %>%
  group_by(sample, Niches) %>%
  summarise(
    tumor_total = sum(tumor_per, na.rm = TRUE),
    immune_total = sum(immune_per, na.rm = TRUE),
    stroma_total = sum(stroma_per, na.rm = TRUE),
    n_spots = n()
  ) %>%
  mutate(
    tumor_per = tumor_total / (tumor_total + immune_total + stroma_total),
    immune_per = immune_total / (tumor_total + immune_total + stroma_total),
    stroma_per = stroma_total / (tumor_total + immune_total + stroma_total)
  ) %>%
  select(sample, Niches, tumor_per, immune_per, stroma_per) %>%
  pivot_longer(cols = -c("sample", "Niches"), values_to = "pct", names_to = "CTA")%>%
  mutate(
    CTA = case_when(
      CTA  == "tumor_per" ~ "Tumor",
      CTA  == "immune_per" ~ "Immune",
      CTA  == "stroma_per" ~ "Stroma",
    )
  )


# Basic stacked bar plot
p = ggplot(meta.plot, aes(x = factor(Niches), y = pct, fill = CTA)) +
  geom_col() +
  facet_wrap(~ sample, nrow = 1) +
  scale_fill_brewer(palette = "Accent") +
  labs(x = "", y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format())+
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    strip.background = element_blank()    
  )
p
ggsave( paste(respath, patient,"B9_niches_CTA_barplot.png", sep = "/" ), p, width = 6, height = 6, dpi = 600, bg = "white")

