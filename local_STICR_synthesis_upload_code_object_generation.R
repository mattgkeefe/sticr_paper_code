library(Seurat)
library(dplyr)
library(ggplot2)


ls1.data = Read10X('GW18_merge/filtered_feature_bc_matrix/')
ls1 = CreateSeuratObject(counts = ls1.data, project = "GW18_0215_3S", min.cells = 10, min.features = 200)
ls1
#15115 features across 1160 samples within 1 assay 
ls1[["percent.mt"]] = PercentageFeatureSet(ls1, pattern="^MT-")
VlnPlot(ls1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls1 <- subset(ls1, subset = 
                nFeature_RNA > 50 & nFeature_RNA < 30000 &
                nFeature_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls1
#15115 features across 1097 samples within 1 assay 


ls3.data = Read10X('GW24_merge/filtered_feature_bc_matrix/')
ls3 = CreateSeuratObject(counts = ls3.data, project = "GW24_0221_1S", min.cells = 10, min.features = 200)
ls3
#18336 features across 6843 samples within 1 assay 
ls3[["percent.mt"]] = PercentageFeatureSet(ls3, pattern="^MT-")
VlnPlot(ls3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls3 <- subset(ls3, subset = 
                nFeature_RNA > 500 & nFeature_RNA < 30000 &
                nFeature_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls3
#18336 features across 6294 samples within 1 assay 

ls11.data = Read10X('GW21_0502_sub_FBS/filtered_feature_bc_matrix')
ls11 = CreateSeuratObject(counts = ls11.data, project = "GW21_0502_1S", min.cells = 10, min.features = 200)
ls11
#14271 features across 549 samples within 1 assay 
ls11[["percent.mt"]] = PercentageFeatureSet(ls11, pattern="^MT-")
VlnPlot(ls11, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls11 <- subset(ls11, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls11, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls11
#14271 features across 430 samples within 1 assay 


ls21.data = Read10X_h5('GW23_0910_sub_FBS_S1_fixed-fastqs_cellranger/filtered_feature_bc_matrix.h5')
ls21 = CreateSeuratObject(counts = ls21.data, project = "GW23_0910_S1", min.cells = 10, min.features = 200)
ls21
#25555 features across 12136 samples within 1 assay 
ls21[["percent.mt"]] = PercentageFeatureSet(ls21, pattern="^MT-")
VlnPlot(ls21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls21 <- subset(ls21, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls21
#25555 features across 11409 samples within 1 assay

ls22.data = Read10X_h5('GW23_0910_sub_FBS_S2L2_fixed-fastqs_cellranger/filtered_feature_bc_matrix.h5')
ls22 = CreateSeuratObject(counts = ls22.data, project = "GW23_0910_S2", min.cells = 10, min.features = 200)
ls22
#21052 features across 3516 samples within 1 assay 
ls22[["percent.mt"]] = PercentageFeatureSet(ls22, pattern="^MT-")
VlnPlot(ls22, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls22 <- subset(ls22, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls22, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls22
#21052 features across 3367 samples within 1 assay 

ls23.data = Read10X_h5('GW23_0910_sub_FBS_S4_fixed-fastqs_cellranger/filtered_feature_bc_matrix.h5')
ls23 = CreateSeuratObject(counts = ls23.data, project = "GW23_0910_S4", min.cells = 10, min.features = 200)
ls23
#23357 features across 5891 samples within 1 assay
ls23[["percent.mt"]] = PercentageFeatureSet(ls23, pattern="^MT-")
VlnPlot(ls23, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls23 <- subset(ls23, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls23, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls23
#23357 features across 5519 samples within 1 assay 

ls27.data = Read10X_h5('GW22_221118_sub_FBS_S1_cellranger/filtered_feature_bc_matrix.h5')
ls27 = CreateSeuratObject(counts = ls27.data, project = "GW22_1118_S1", min.cells = 10, min.features = 200)
ls27
#21080 features across 3633 samples within 1 assay
ls27[["percent.mt"]] = PercentageFeatureSet(ls27, pattern="^MT-")
VlnPlot(ls27, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls27 <- subset(ls27, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls27, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls27
#21080 features across 3418 samples within 1 assay

ls28.data = Read10X_h5('GW22_221118_sub_FBS_S2_cellranger/filtered_feature_bc_matrix.h5')
ls28 = CreateSeuratObject(counts = ls28.data, project = "GW22_1118_S2", min.cells = 10, min.features = 200)
ls28
#20776 features across 3280 samples within 1 assay
ls28[["percent.mt"]] = PercentageFeatureSet(ls28, pattern="^MT-")
VlnPlot(ls28, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls28 <- subset(ls28, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls28, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls28
#20776 features across 3142 samples within 1 assay  

ls29.data = Read10X_h5('GW22_221118_sub_FBS_S3_cellranger/filtered_feature_bc_matrix.h5')
ls29 = CreateSeuratObject(counts = ls29.data, project = "GW22_1118_S3", min.cells = 10, min.features = 200)
ls29
#19504 features across 2480 samples within 1 assay
ls29[["percent.mt"]] = PercentageFeatureSet(ls29, pattern="^MT-")
VlnPlot(ls29, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls29 <- subset(ls29, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls29, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls29
#19504 features across 2362 samples within 1 assay

ls30.data = Read10X_h5('GW16_1220_sub_FBS_S1_cellranger/filtered_feature_bc_matrix.h5')
ls30 = CreateSeuratObject(counts = ls30.data, project = "GW16_1220_S1", min.cells = 10, min.features = 200)
ls30
#24587 features across 11123 samples within 1 assay
ls30[["percent.mt"]] = PercentageFeatureSet(ls30, pattern="^MT-")
VlnPlot(ls30, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls30 <- subset(ls30, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 1500 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls30, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls30
#24587 features across 9816 samples within 1 assay

ls31.data = Read10X_h5('GW16_1220_sub_FBS_S2_cellranger/filtered_feature_bc_matrix.h5')
ls31 = CreateSeuratObject(counts = ls31.data, project = "GW16_1220_S2", min.cells = 10, min.features = 200)
ls31
#23581 features across 5577 samples within 1 assay
ls31[["percent.mt"]] = PercentageFeatureSet(ls31, pattern="^MT-")
VlnPlot(ls31, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls31 <- subset(ls31, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 3000 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls31, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls31
#23581 features across 4250 samples within 1 assay

ls32.data = Read10X_h5('GW19_0221_sub_FBS_S1_cellranger/filtered_feature_bc_matrix.h5')
ls32 = CreateSeuratObject(counts = ls32.data, project = "GW19_0221_S1", min.cells = 10, min.features = 200)
ls32
#24862 features across 9376 samples within 1 assay 
ls32[["percent.mt"]] = PercentageFeatureSet(ls32, pattern="^MT-")
VlnPlot(ls32, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls32 <- subset(ls32, subset = 
                 nFeature_RNA > 50 & nFeature_RNA < 30000 &
                 nFeature_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls32, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls32
#24862 features across 7854 samples within 1 assay  

ls33.data = Read10X_h5('GW19_0221_sub_FBS_S2_cellranger/filtered_feature_bc_matrix.h5')
ls33 = CreateSeuratObject(counts = ls33.data, project = "GW19_0221_S2", min.cells = 10, min.features = 200)
ls33
#22614 features across 4327 samples within 1 assay 
ls33[["percent.mt"]] = PercentageFeatureSet(ls33, pattern="^MT-")
VlnPlot(ls33, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls33 <- subset(ls33, subset = 
                 nFeature_RNA > 50 & nFeature_RNA < 30000 &
                 nFeature_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls33, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls33
#22614 features across 3211 samples within 1 assay 

ls34.data = Read10X_h5('GW20_230308_sub_FBS_V1_S1_cellranger/filtered_feature_bc_matrix.h5')
ls34 = CreateSeuratObject(counts = ls34.data, project = "GW20_0308_V1_S1", min.cells = 10, min.features = 200)
ls34
#25042 features across 14694 samples within 1 assay 
ls34[["percent.mt"]] = PercentageFeatureSet(ls34, pattern="^MT-")
VlnPlot(ls34, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls34 <- subset(ls34, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls34, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls34
#25042 features across 14276 samples within 1 assay 


ls35.data = Read10X_h5('GW20_230308_sub_FBS_V1_S2_cellranger/filtered_feature_bc_matrix.h5')
ls35 = CreateSeuratObject(counts = ls35.data, project = "GW20_0308_V1_S2", min.cells = 10, min.features = 200)
ls35
#25364 features across 17590 samples within 1 assay 
ls35[["percent.mt"]] = PercentageFeatureSet(ls35, pattern="^MT-")
VlnPlot(ls35, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls35 <- subset(ls35, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls35, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls35
#25364 features across 17387 samples within 1 assay 

ls36.data = Read10X_h5('GW20_230308_sub_FBS_PFC_2S_cellranger/filtered_feature_bc_matrix.h5')
ls36 = CreateSeuratObject(counts = ls36.data, project = "GW20_0308_PFC_2S", min.cells = 10, min.features = 200)
ls36
#25260 features across 11694 samples within 1 assay 
ls36[["percent.mt"]] = PercentageFeatureSet(ls36, pattern="^MT-")
VlnPlot(ls36, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls36 <- subset(ls36, subset = 
                 nFeature_RNA > 500 & nFeature_RNA < 30000 &
                 nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(ls36, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ls36
#25260 features across 11422 samples within 1 assay 


ls_synthesis_v3 = merge(ls1, c(ls3,
                               ls11,
                               ls21,
                               ls22,
                               ls23,
                               ls27,
                               ls28,
                               ls29,
                               ls30,
                               ls31,
                               ls32,
                               ls33,
                               ls34,
                               ls35,
                               ls36
),
project="local_STICR_synthesis_v3")

ls_synthesis_v3
#27581 features across 105254 samples within 1 assay 
table(ls_synthesis_v3@meta.data$orig.ident)
#GW16_1220_S1     GW16_1220_S2     GW18_0215_3S     GW19_0221_S1     GW19_0221_S2 GW20_0308_PFC_2S 
#        9816             4250             1097             7854             3211            11422 
#GW20_0308_V1_S1  GW20_0308_V1_S2     GW21_0502_1S     GW22_1118_S1     GW22_1118_S2     GW22_1118_S3 
#          14276            17387              430             3418             3142             2362 
#GW23_0910_S1     GW23_0910_S2     GW23_0910_S4     GW24_0221_1S 
#       11409             3367             5519             6294 

# normalize data
ls_synthesis_v3 <- NormalizeData(ls_synthesis_v3, normalization.method = "LogNormalize", scale.factor = 10000)
# scale data based on normalization
all.genes <- rownames(ls_synthesis_v3)
ls_synthesis_v3 <- ScaleData(ls_synthesis_v3, features = all.genes)
#ls_synthesis_v3 <- ScaleData(ls_synthesis_v3, features = all.genes, vars.to.regress = "nCount_RNA")
# find variably expressed genes
ls_synthesis_v3 <- FindVariableFeatures(ls_synthesis_v3, selection.method = "vst", nfeatures = 2000)
# perform dimensionality reduction
ls_synthesis_v3 <- RunPCA(ls_synthesis_v3, features = VariableFeatures(object = ls_synthesis_v3))
ElbowPlot(ls_synthesis_v3)
# cluster cells
ls_synthesis_v3 <- FindNeighbors(ls_synthesis_v3, dims = 1:15)
ls_synthesis_v3 <- FindClusters(ls_synthesis_v3, resolution = .4)
# make graph
ls_synthesis_v3 <- RunUMAP(ls_synthesis_v3, dims = 1:15)

DimPlot(ls_synthesis_v3, reduction = "umap",label=T, raster=FALSE)
DimPlot(ls_synthesis_v3,group.by = "orig.ident", raster=FALSE)
FeaturePlot(ls_synthesis_v3, features=c("SATB2", "GAD2", "EOMES", "VIM"))
FeaturePlot(ls_synthesis_v3, features=c("TBR1", "BCL11B", "SATB2", "NEUROD2"))

FeaturePlot(ls_synthesis_v3, features=c("XIST"))
#FeaturePlot(ls_synthesis_v3, features=c("XIST"), split.by="orig.ident")
VlnPlot(ls_synthesis_v3, features=c("XIST"), group.by="orig.ident")

#add sex as metadata feature
dim(ls_synthesis_v3@meta.data)
head(ls_synthesis_v3@meta.data)
unique(ls_synthesis_v3@meta.data$orig.ident)
ls_synthesis_v3 = AddMetaData(ls_synthesis_v3, rep('XY', length(rownames(ls_synthesis_v3@meta.data))), 'sex')
head(ls_synthesis_v3@meta.data)
ls_synthesis_v3@meta.data$sex[ls_synthesis_v3@meta.data$orig.ident=="GW18_0215_3S"]="XX"
ls_synthesis_v3@meta.data$sex[ls_synthesis_v3@meta.data$orig.ident=="GW21_0502_1S"]="XX"
ls_synthesis_v3@meta.data$sex[ls_synthesis_v3@meta.data$orig.ident=="GW24_0221_1S"]="XX"

DimPlot(ls_synthesis_v3,group.by = "sex", raster=F)

## run harmony on orig.ident
ls_synthesis_v3_harmony_v2 <- RunHarmony(ls_synthesis_v3, "orig.ident", kmeans_init_iter_max=100, kmeans_init_nstart=50, plot_convergence = TRUE, project.dim = F)
ls_synthesis_v3_harmony_v2 <- FindNeighbors(ls_synthesis_v3_harmony_v2, reduction = "harmony", dims = 1:15)
ls_synthesis_v3_harmony_v2 <- FindClusters(ls_synthesis_v3_harmony_v2, resolution = .4)
ls_synthesis_v3_harmony_v2 <- RunUMAP(ls_synthesis_v3_harmony_v2, reduction = "harmony", dims=c(1:15))
DimPlot(ls_synthesis_v3_harmony_v2, reduction = "umap",label=T, raster=F)
## inverting the UMAP embeddings so it looks more like ls_synthesis_v3_harmony_v1
ls_synthesis_v3_harmony_v2[["umap"]]@cell.embeddings = -ls_synthesis_v3_harmony_v2[["umap"]]@cell.embeddings
DimPlot(ls_synthesis_v3_harmony_v2, reduction = "umap",label=T, raster=F)
DimPlot(ls_synthesis_v3_harmony_v2,group.by = "orig.ident", raster=F)
DimPlot(ls_synthesis_v3_harmony_v2,group.by = "sex", raster=F)

unique(ls_synthesis_v3_harmony_v2@meta.data$orig.ident)
#[1] "GW18_0215_3S"     "GW24_0221_1S"     "GW21_0502_1S"     "GW23_0910_S1"     "GW23_0910_S2"    
#[6] "GW23_0910_S4"     "GW22_1118_S1"     "GW22_1118_S2"     "GW22_1118_S3"     "GW16_1220_S1"    
#[11] "GW16_1220_S2"     "GW19_0221_S1"     "GW19_0221_S2"     "GW20_0308_V1_S1"  "GW20_0308_V1_S2" 
#[16] "GW20_0308_PFC_2S"

FeaturePlot(ls_synthesis_v3_harmony_v2, features=c("percent.mt"), raster=F)
FeaturePlot(ls_synthesis_v3_harmony_v2, features=c("nFeature_RNA"), raster=F)

FeaturePlot(ls_synthesis_v3_harmony_v2, features=c("SATB2", "OLIG2", "GAD2", "MKI67"), raster=F)

VlnPlot(ls_synthesis_v3_harmony_v2, features = c("nCount_RNA", "percent.mt"), ncol = 2, pt.size = 0)
table(ls_synthesis_v3_harmony_v2@meta.data$seurat_clusters)
#     0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16 
#16134 15102 14824 11283  8814  6760  6594  4726  4382  4239  3514  2922  2530  1168   921   888   453 

## remove cluster 10 and 15 based on nCountRNA
subset(ls_synthesis_v3_harmony_v2, seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,11,12,13,14,16))
saveRDS(ls_synthesis_v3_harmony_v2, '230402_ls_synthesis_v3_harmony_v2_obj.rds')
ls_synthesis_v3_harmony_v3= subset(ls_synthesis_v3_harmony_v2, seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,11,12,13,14,16))
ls_synthesis_v3_harmony_v3

#rerun clustering after removing those cells
ls_synthesis_v3_harmony_v3 <- FindNeighbors(ls_synthesis_v3_harmony_v3, reduction = "harmony", dims = 1:15)
ls_synthesis_v3_harmony_v3 <- FindClusters(ls_synthesis_v3_harmony_v3, resolution = .4)
ls_synthesis_v3_harmony_v3 <- RunUMAP(ls_synthesis_v3_harmony_v3, reduction = "harmony", dims=c(1:15))
DimPlot(ls_synthesis_v3_harmony_v3, reduction = "umap",label=T, raster=F)
#need to invert just UMAP_1 to restore same positioning
ls_synthesis_v3_harmony_v3[["umap"]]@cell.embeddings[,1] = -ls_synthesis_v3_harmony_v3[["umap"]]@cell.embeddings[,1]
DimPlot(ls_synthesis_v3_harmony_v3, reduction = "umap",label=T, raster=F)
VlnPlot(ls_synthesis_v3_harmony_v3, features = c("nCount_RNA", "percent.mt"), ncol = 2, pt.size = 0)
#remove low quality cells in cluster 15
table(ls_synthesis_v3_harmony_v3@meta.data$seurat_clusters)
#    0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15 
#16449 14403 13755 10345  8120  6783  6618  6026  5558  4119  3103  2793  1120   921   458   281
ls_synthesis_v3_harmony_v3= subset(ls_synthesis_v3_harmony_v3, seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14))
ls_synthesis_v3_harmony_v3 <- FindNeighbors(ls_synthesis_v3_harmony_v3, reduction = "harmony", dims = 1:15)
ls_synthesis_v3_harmony_v3 <- FindClusters(ls_synthesis_v3_harmony_v3, resolution = .4)
ls_synthesis_v3_harmony_v3 <- RunUMAP(ls_synthesis_v3_harmony_v3, reduction = "harmony", dims=c(1:15))
DimPlot(ls_synthesis_v3_harmony_v3, reduction = "umap",label=T, raster=F)
#need to invert just UMAP_2 to restore same positioning
ls_synthesis_v3_harmony_v3[["umap"]]@cell.embeddings[,2] = -ls_synthesis_v3_harmony_v3[["umap"]]@cell.embeddings[,2]
DimPlot(ls_synthesis_v3_harmony_v3, reduction = "umap",label=T, raster=F)
VlnPlot(ls_synthesis_v3_harmony_v3, features = c("nCount_RNA", "percent.mt"), ncol = 2, pt.size = 0)

# find cluster markers
ls_synthesis_v3_harmony_v3.markers <- FindAllMarkers(ls_synthesis_v3_harmony_v3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(ls_synthesis_v3_harmony_v3.markers)
ls_synthesis_v3_harmony_v3_clustermarkers_unfilt = as.data.frame(ls_synthesis_v3_harmony_v3.markers %>% group_by(cluster))
write.csv(ls_synthesis_v3_harmony_v3_clustermarkers_unfilt,'ls_synthesis_v3_harmony_v3_clustermarkers_v1_unfilt.csv',row.names=F)
ls_synthesis_v3_harmony_v3.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
ls_synthesis_v3_harmony_v3_clustermarkers = as.data.frame(ls_synthesis_v3_harmony_v3.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
write.csv(ls_synthesis_v3_harmony_v3_clustermarkers,'ls_synthesis_v3_harmony_v3_clustermarkers_v1.csv',row.names=F)


#name clusters and assign metadata for easier future parsing
clust_idents = c("EN_0","EN_early","EN_2","EN_UL","OPC","DIV","RG","IPC","IN_local","IN_CGE","IN_MGE","EN_11","OLIG","RG_unk")
full_clust_list = ls_synthesis_v3_harmony_v3@meta.data$seurat_clusters
full_clust_idents = rep("EN_DL", length(full_clust_list))
for (s in c(0:13)) {
  full_clust_idents[full_clust_list==s] = clust_idents[s+1]
}
full_clust_idents = factor(full_clust_idents, levels=c("EN_0","EN_early","EN_2","EN_UL","OPC","DIV","RG","IPC","IN_local","IN_CGE","IN_MGE","EN_11","OLIG","RG_unk"))
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, full_clust_idents, col.name = "clust_idents")
ls_synthesis_v3_harmony_v3 = SetIdent(ls_synthesis_v3_harmony_v3, value="clust_idents")
DimPlot(ls_synthesis_v3_harmony_v3, label=T, raster=F)


#add pseudoage as a metadata feature
ls_synthesis_v3 = AddMetaData(ls_synthesis_v3, rep('>GW20', length(rownames(ls_synthesis_v3@meta.data))), 'age_pseudo')
unique(ls_synthesis_v3_harmony_v3@meta.data$orig.ident)
ls_synthesis_v3_harmony_v3@meta.data$age_pseudo[ls_synthesis_v3_harmony_v3@meta.data$orig.ident=="GW18_0215_3S"]="<=GW20"
ls_synthesis_v3_harmony_v3@meta.data$age_pseudo[ls_synthesis_v3_harmony_v3@meta.data$orig.ident=="GW16_1220_S1"]="<=GW20"
ls_synthesis_v3_harmony_v3@meta.data$age_pseudo[ls_synthesis_v3_harmony_v3@meta.data$orig.ident=="GW16_1220_S2"]="<=GW20"
ls_synthesis_v3_harmony_v3@meta.data$age_pseudo[ls_synthesis_v3_harmony_v3@meta.data$orig.ident=="GW19_0221_S1"]="<=GW20"
ls_synthesis_v3_harmony_v3@meta.data$age_pseudo[ls_synthesis_v3_harmony_v3@meta.data$orig.ident=="GW19_0221_S2"]="<=GW20"
ls_synthesis_v3_harmony_v3@meta.data$age_pseudo[ls_synthesis_v3_harmony_v3@meta.data$orig.ident=="GW20_0308_V1_S1"]="<=GW20"
ls_synthesis_v3_harmony_v3@meta.data$age_pseudo[ls_synthesis_v3_harmony_v3@meta.data$orig.ident=="GW20_0308_V1_S2"]="<=GW20"
ls_synthesis_v3_harmony_v3@meta.data$age_pseudo[ls_synthesis_v3_harmony_v3@meta.data$orig.ident=="GW20_0308_PFC_2S"]="<=GW20"
table(ls_synthesis_v3_harmony_v3$age_pseudo)
#<=GW20  >GW20 
#65917  34654 

## add in STICR barcodes (with read count data)
head(ls_synthesis_v3_harmony_v3@meta.data)
unique(ls_synthesis_v3_harmony_v3@meta.data$orig.ident)
seurat_meta = ls_synthesis_v3_harmony_v3@meta.data
seurat_meta$CBC = rownames(seurat_meta)
dim(seurat_meta)
#[1] 100571     10

#9/10 S1
sample_loadpath = 'GW23_0910_sub_FBS_S1_scAR_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW23_0910_S1"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-3,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = sample_bcs
#9/10 S2
sample_loadpath = 'GW23_0910_sub_FBS_S2_scAR_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW23_0910_S2"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-3,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#9/10 S4
sample_loadpath = 'GW23_0910_sub_FBS_S4_scAR_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW23_0910_S4"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-3,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#5/02 sub FBS 
sample_loadpath = 'GW21_0502_sub_FBS_scAR_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW21_0502_1S"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-3,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#2/21 GW24 merge
sample_loadpath = 'GW24_0324_sub_FBS_S1_merge_scAR_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW24_0221_1S"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-3,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#12/20 GW16 S1
#note -- have to change the code here because there are two digits in the CBC suffix!
sample_loadpath = 'GW16_221220_S1_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW16_1220_S1"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-4,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#12/20 GW16 S2
#note -- have to change the code here because there are two digits in the CBC suffix!
sample_loadpath = 'GW16_221220_S2_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW16_1220_S2"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-4,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#11/18 GW22 S1
sample_loadpath = 'GW22_221118_S1_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW22_1118_S1"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-3,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#11/18 GW22 S2
sample_loadpath = 'GW22_221118_S2_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW22_1118_S2"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-3,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#11/18 GW22 S3
sample_loadpath = 'GW22_221118_S3_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW22_1118_S3"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-3,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#2/21/23 GW19 S1
sample_loadpath = 'GW19_0221_S1_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW19_0221_S1"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-4,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#2/21/23 GW19 S2
sample_loadpath = 'GW19_0221_S2_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW19_0221_S2"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-4,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#GW20_0308_PFC_2S
sample_loadpath = 'GW20_0308_PFC_2S_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW20_0308_PFC_2S"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-4,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#GW20_0308_V1_S1
sample_loadpath = 'GW20_0308_V1_S1_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW20_0308_V1_S1"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-4,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)
#GW20_0308_V1_S2
sample_loadpath = 'GW20_0308_V1_S2_bcs_final_metadata_UMIs_reads.csv'
sample_id = "GW20_0308_V1_S2"
sample_bcs = read.csv(sample_loadpath, row.names = 1)
head(sample_bcs)
char_len=nchar(seurat_meta$CBC[seurat_meta$orig.ident==sample_id])[1]
cbc_suffix=substr(seurat_meta$CBC[seurat_meta==sample_id][1],char_len-4,char_len)
sample_bcs$CBC = paste(sample_bcs$CBC, cbc_suffix, sep="")
head(sample_bcs)
master_sticr_umi_df = rbind(master_sticr_umi_df, sample_bcs)

seurat_meta_sticr_umi = left_join(seurat_meta, master_sticr_umi_df, by='CBC')
dim(seurat_meta_sticr_umi)
dim(ls_synthesis_v3_harmony_v3@meta.data)
rownames(seurat_meta_sticr_umi) = seurat_meta_sticr_umi$CBC
sum(is.na(seurat_meta_sticr_umi$Clone_barcodes))
table(seurat_meta_sticr_umi$orig.ident[is.na(seurat_meta_sticr_umi$Clone_barcodes)]) #check distribution of missing BCs
#GW16_1220_S1     GW16_1220_S2     GW18_0215_3S     GW19_0221_S1     GW19_0221_S2 GW20_0308_PFC_2S 
#.        325              102              955              891              547              771 
#GW20_0308_V1_S1  GW20_0308_V1_S2    GW21_0502_1S     GW22_1118_S1     GW22_1118_S2     GW22_1118_S3 
#.           676             830               63               99               86               90 
#GW23_0910_S1     GW23_0910_S2     GW23_0910_S4     GW24_0221_1S 
#        314              193              384              684 
colnames(seurat_meta_sticr_umi)
#[1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.4"
#[6] "seurat_clusters" "sex"             "age_pseudo"      "RNA_snn_res.0.6" "RNA_snn_res.0.5"
#[11] "CBC"             "Clone_barcodes"  "n_barcodes"      "Clone_size"      "Clone_IDs"      
#[16] "Clone_Index"     "UMI_Count"       "CBC_barcode"     "read_counts"     "thresh_pass" 
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, 
                                         seurat_meta_sticr_umi$Clone_barcodes, col.name='Clone_barcodes')
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, 
                                         seurat_meta_sticr_umi$Clone_size, col.name='Clone_size')
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, 
                                         seurat_meta_sticr_umi$Clone_IDs, col.name='Clone_IDs')
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, 
                                         seurat_meta_sticr_umi$Clone_Index, col.name='Clone_Index')
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, 
                                         seurat_meta_sticr_umi$UMI_Count, col.name='Barcode_UMI_Count')
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, 
                                         seurat_meta_sticr_umi$read_counts, col.name='read_counts')

#fix read_thresh_pass to make it Boolean:
read_thresh_bool = rep(FALSE, length(seurat_meta_sticr_umi$thresh_pass))
read_thresh_bool[seurat_meta_sticr_umi$thresh_pass == "True"] = TRUE
sum(read_thresh_bool)
#1] 72658
length(read_thresh_bool)
#[1] 100571
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, 
                                         read_thresh_bool, col.name='read_thresh_pass')

table(ls_synthesis_v3_harmony_v3$Barcode_UMI_Count)
#    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18 
#27359  6228  4367  3785  3526  3387  3114  2902  2581  2366  2281  2102  2083  1896  1732  1549  1530  1380
hist(log10(ls_synthesis_v3_harmony_v3$Barcode_UMI_Count))
hist(log10(ls_synthesis_v3_harmony_v3@meta.data$Barcode_UMI_Count[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass]))



ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, rep(NA, dim(ls_synthesis_v3_harmony_v3@meta.data)[1]), col.name='Clone_barcodes_filt')
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, rep(NA, dim(ls_synthesis_v3_harmony_v3@meta.data)[1]), col.name='Clone_size_filt')
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, rep(NA, dim(ls_synthesis_v3_harmony_v3@meta.data)[1]), col.name='Clone_IDs_filt')
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, rep(NA, dim(ls_synthesis_v3_harmony_v3@meta.data)[1]), col.name='Clone_Index_filt')
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, rep(NA, dim(ls_synthesis_v3_harmony_v3@meta.data)[1]), col.name='Barcode_UMI_Count_filt')
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, rep(NA, dim(ls_synthesis_v3_harmony_v3@meta.data)[1]), col.name='read_counts_filt')


ls_synthesis_v3_harmony_v3@meta.data$Clone_barcodes_filt[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass] = ls_synthesis_v3_harmony_v3@meta.data$Clone_barcodes[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass]
ls_synthesis_v3_harmony_v3@meta.data$Clone_size_filt[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass] = ls_synthesis_v3_harmony_v3@meta.data$Clone_size[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass]
ls_synthesis_v3_harmony_v3@meta.data$Clone_IDs_filt[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass] = ls_synthesis_v3_harmony_v3@meta.data$Clone_IDs[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass]
ls_synthesis_v3_harmony_v3@meta.data$Clone_Index_filt[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass] = ls_synthesis_v3_harmony_v3@meta.data$Clone_Index[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass]
ls_synthesis_v3_harmony_v3@meta.data$Barcode_UMI_Count_filt[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass] = ls_synthesis_v3_harmony_v3@meta.data$Barcode_UMI_Count[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass]
ls_synthesis_v3_harmony_v3@meta.data$read_counts_filt[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass] = ls_synthesis_v3_harmony_v3@meta.data$read_counts[ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass]

sum(ls_synthesis_v3_harmony_v3@meta.data$read_thresh_pass)
length(ls_synthesis_v3_harmony_v3@meta.data$Clone_Index)
sum(is.na(ls_synthesis_v3_harmony_v3@meta.data$Clone_Index))
sum(is.na(ls_synthesis_v3_harmony_v3@meta.data$Clone_Index_filt))

table(ls_synthesis_v3_harmony_v3@meta.data$Clone_Index)
#Index3 IndexE 
#44995  48566
table(ls_synthesis_v3_harmony_v3@meta.data$Clone_Index_filt)
#Index3 IndexE 
#34915  37743 

sum(is.na(ls_synthesis_v3_harmony_v3@meta.data$Clone_size_filt))/length(ls_synthesis_v3_harmony_v3@meta.data$Clone_size_filt)
#[1] 0.2775452
sum(!is.na(ls_synthesis_v3_harmony_v3@meta.data$Clone_size_filt))/length(ls_synthesis_v3_harmony_v3@meta.data$Clone_size_filt)
#[1] 0.7224548


## try to start subclustering to get more specific cell identities
ls_div = subset(ls_synthesis_v3_harmony_v3, clust_idents %in% "DIV")
ls_div
#27581 features across 7240 samples within 1 assay 
ls_div <- ScaleData(ls_div, features = all.genes)
ls_div = FindVariableFeatures(ls_div, selection.method="vst", nfeatures=2000)
ls_div = RunPCA(ls_div, features=VariableFeatures(object=ls_div))
ElbowPlot(ls_div)
ls_div <- RunHarmony(ls_div, c("orig.ident"))
ls_div <- FindNeighbors(ls_div, reduction = "harmony", dims = 1:15)
ls_div <- FindClusters(ls_div, resolution = .5)
ls_div <- RunUMAP(ls_div, reduction = "harmony", dims=c(1:15))
DimPlot(ls_div, reduction = "umap",label=T)
DimPlot(ls_div, reduction = "umap",group.by="orig.ident")
DimPlot(ls_div, reduction = "umap",group.by="age_pseudo")
FeaturePlot(ls_div, features=c("VIM", "HOPX", "EOMES", "DLX2", "OLIG2", "NEUROD2", "ID3", "PPP1R17", "MKI67"))
VlnPlot(ls_div, features = c("nCount_RNA", "percent.mt"), ncol = 2)
VlnPlot(ls_div, features = c("DLX2", "EOMES", "VIM", "OLIG2"), ncol=4)
#subcluster assignments:
# 0, 6, 9: IPC_EX_DIV
# 1, 2, 8: GLIA_DIV
# 7: OPC_DIV
# 3, 4, 5: IPC_IN_DIV
ls_div$subclust_idents = NA
ls_div$subclust_idents[ls_div$seurat_clusters %in% c(0,6,9)] = "IPC_EX_DIV"
ls_div$subclust_idents[ls_div$seurat_clusters %in% c(1,2,8)] = "GLIA_DIV"
ls_div$subclust_idents[ls_div$seurat_clusters %in% c(7)] = "OPC_DIV"
ls_div$subclust_idents[ls_div$seurat_clusters %in% c(3,4,5)] = "IPC_IN_DIV"
table(ls_div$subclust_idents)
#GLIA_DIV IPC_EX_DIV IPC_IN_DIV    OPC_DIV
#    2303       2060       2334        543 
ls_div = SetIdent(ls_div, value="subclust_idents")
DimPlot(ls_div, label=T)


ls_glia = subset(ls_synthesis_v3_harmony_v3, clust_idents %in% c("RG", "OPC", "OLIG", "RG_unk"))
ls_glia
#27581 features across 15778 samples within 1 assay
ls_glia <- ScaleData(ls_glia, features = all.genes)
ls_glia = FindVariableFeatures(ls_glia, selection.method="vst", nfeatures=2000)
ls_glia = RunPCA(ls_glia, features=VariableFeatures(object=ls_glia))
ElbowPlot(ls_glia)
ls_glia <- RunHarmony(ls_glia, c("orig.ident"))
ls_glia <- FindNeighbors(ls_glia, reduction = "harmony", dims = 1:15)
ls_glia <- FindClusters(ls_glia, resolution = .5)
ls_glia <- RunUMAP(ls_glia, reduction = "harmony", dims=c(1:15))
DimPlot(ls_glia, reduction = "umap",label=T, pt.size=1)
VlnPlot(ls_glia, features = c("nCount_RNA", "percent.mt"), ncol = 2)
DimPlot(ls_glia, reduction = "umap",group.by="orig.ident", pt.size=1)
DimPlot(ls_glia, reduction = "umap",group.by="age_pseudo", pt.size=1)
FeaturePlot(ls_glia, label=T, features=c("VIM", "HOPX", "ID3", "DLX5", "OLIG2", "MBP", "EOMES", "PPP1R17", "AQP4"))
FeaturePlot(ls_glia, label=T, features=c("VIM", "NEUROD2", "EOMES", "TBR1"))
FeaturePlot(ls_glia, label=T, features=c("nCount_RNA"))
VlnPlot(ls_glia, features = c("VIM", "PPP1R17", "OLIG2", "DLX5"), ncol=4)
#subcluster assignments
#0, 1, 2, 5, 7, 11, 12: GLIA
#3, 9: OPC
#8: OLIG
#6: IPC_IN
#10: IPC_EX
#4: RG_EX ##possibly doublets?
#13: NA
ls_glia$subclust_idents = NA
ls_glia$subclust_idents[ls_glia$seurat_clusters %in% c(0,1,2,5,7,11,12)] = "GLIA"
ls_glia$subclust_idents[ls_glia$seurat_clusters %in% c(3,9)] = "OPC"
ls_glia$subclust_idents[ls_glia$seurat_clusters %in% c(8)] = "OLIG"
ls_glia$subclust_idents[ls_glia$seurat_clusters %in% c(6)] = "IPC_IN"
ls_glia$subclust_idents[ls_glia$seurat_clusters %in% c(10)] = "IPC_EX"
ls_glia$subclust_idents[ls_glia$seurat_clusters %in% c(4)] = "RG_EX"
table(ls_glia$subclust_idents)
#GLIA IPC_EX IPC_IN   OLIG    OPC  RG_EX 
#9344    709   1192    755   2249   1446 
ls_glia = SetIdent(ls_glia, value="subclust_idents")
DimPlot(ls_glia, label=T)


table(ls_synthesis_v3_harmony_v3$clust_idents)
ls_in = subset(ls_synthesis_v3_harmony_v3, clust_idents %in% c("IN_CGE", "IN_MGE", "IN_local"))
ls_in
#27581 features across 12025 samples within 1 assay
ls_in <- ScaleData(ls_in, features = all.genes)
ls_in = FindVariableFeatures(ls_in, selection.method="vst", nfeatures=2000)
ls_in = RunPCA(ls_in, features=VariableFeatures(object=ls_in))
ElbowPlot(ls_in)
ls_in <- RunHarmony(ls_in, c("orig.ident"))
ls_in <- FindNeighbors(ls_in, reduction = "harmony", dims = 1:15)
ls_in <- FindClusters(ls_in, resolution = .5)
ls_in <- RunUMAP(ls_in, reduction = "harmony", dims=c(1:15))
DimPlot(ls_in, reduction = "umap",label=T, pt.size=1)
VlnPlot(ls_in, features = c("nCount_RNA", "percent.mt"), ncol = 2)
DimPlot(ls_in, reduction = "umap",group.by="orig.ident", pt.size=1)
DimPlot(ls_in, reduction = "umap",group.by="age_pseudo", pt.size=1)
FeaturePlot(ls_in, label=T, features=c("GAD2", "NEUROD2", "NR2F1", "PAX6", "MKI67", "NPY", "NR2F2", "ADARB2", "MAF"))
FeaturePlot(ls_in, label=T, features=c("nCount_RNA"))
DimPlot(ls_in, label=T, cells.highlight=rownames(ls_in@meta.data[ls_in@meta.data$"Clone_size_filt">=3,])) +ggtitle("Clone size >=3")
#subcluster assignments
# 0, 8: IN_CGE
# 2, 6: IN_MGE
# 7: IN_unk ##possibly doublets?
# 1, 3, 4, 5, 9: IN_local
# 10: IPC_IN_DIV
ls_in$subclust_idents = NA
ls_in$subclust_idents[ls_in$seurat_clusters %in% c(1, 3, 4, 5, 9)] = "IN_local"
ls_in$subclust_idents[ls_in$seurat_clusters %in% c(0, 8)] = "IN_CGE"
ls_in$subclust_idents[ls_in$seurat_clusters %in% c(2, 6)] = "IN_MGE"
ls_in$subclust_idents[ls_in$seurat_clusters %in% c(7)] = "IN_UNK"
ls_in$subclust_idents[ls_in$seurat_clusters %in% c(10)] = "IPC_IN_DIV"
ls_in = SetIdent(ls_in, value="subclust_idents")
DimPlot(ls_in, label=T)

table(ls_synthesis_v3_harmony_v3$clust_idents)
ls_ex = subset(ls_synthesis_v3_harmony_v3, clust_idents %in% c("EN_0", "EN_early", "EN_UL", "IPC", "EN_11", "EN_2"))
ls_ex
#27581 features across 65528 samples within 1 assay 
ls_ex <- ScaleData(ls_ex, features = all.genes)
ls_ex = FindVariableFeatures(ls_ex, selection.method="vst", nfeatures=2000)
ls_ex = RunPCA(ls_ex, features=VariableFeatures(object=ls_ex))
ElbowPlot(ls_ex)
ls_ex <- RunHarmony(ls_ex, c("orig.ident"))
ls_ex <- FindNeighbors(ls_ex, reduction = "harmony", dims = 1:15)
ls_ex <- FindClusters(ls_ex, resolution = .5)
ls_ex <- RunUMAP(ls_ex, reduction = "harmony", dims=c(1:15))
DimPlot(ls_ex, reduction = "umap",label=T)
VlnPlot(ls_ex, features = c("nCount_RNA", "percent.mt"), ncol = 2)
DimPlot(ls_ex, reduction = "umap",group.by="orig.ident", pt.size=0.5)
DimPlot(ls_ex, reduction = "umap",group.by="age_pseudo", pt.size=0.5)
FeaturePlot(ls_ex, label=T, features=c("EOMES", "NEUROD2", "MARCH1", "TBR1", "BCL11B", "SATB2", "CUX2", "CUX1", "TLE4"))
FeaturePlot(ls_ex, label=T, features=c("TLE4", "NR4A2", "TBR1", "OTX1", "BCL11B", "RORB", "SATB2", "POU3F2", "CUX2"))
FeaturePlot(ls_ex, label=T, features=c("VIM", "HOPX", "NEUROD1", "NEUROD2"))
DimPlot(ls_ex, reduction = "umap", group.by="Clone_Index_filt", pt.size=0.5)
DimPlot(ls_ex, reduction = "umap", split.by="Clone_Index_filt", pt.size=0.5)
#subcluster assignments
#0, 1, 2, 3, 4, 7, 9: EX
#5, 10: EX_early
#6, 8: IPC_EX
ls_ex$subclust_idents = NA
ls_ex$subclust_idents[ls_ex$seurat_clusters %in% c(0, 1, 2, 3, 4, 7, 9)] = "EX"
ls_ex$subclust_idents[ls_ex$seurat_clusters %in% c(5, 10)] = "EX_early"
ls_ex$subclust_idents[ls_ex$seurat_clusters %in% c(6, 8)] = "IPC_EX"
ls_ex = SetIdent(ls_ex, value="subclust_idents")
DimPlot(ls_ex, label=T)

ls_ex_unmerge = FindNeighbors(ls_ex, dims=1:15)
ls_ex_unmerge = FindClusters(ls_ex_unmerge, resolution=.5)
ls_ex_unmerge = RunUMAP(ls_ex_unmerge, dims=c(1:15))
DimPlot(ls_ex_unmerge, label=T)
DimPlot(ls_ex_unmerge, group.by="orig.ident")
DimPlot(ls_ex_unmerge, group.by="age_pseudo")
DimPlot(ls_ex_unmerge, group.by="sex")
DimPlot(ls_ex_unmerge, cells.highlight = rownames(ls_ex_unmerge@meta.data[ls_ex_unmerge@meta.data$Clone_size_corrected>=3,]))
FeaturePlot(ls_ex_unmerge, label=T, features=c("EOMES", "NEUROD2", "MARCH1", "TBR1", "BCL11B", "SATB2", "CUX2", "CUX1", "TLE4"))


## compile all subclust_idents and add to the main object
all_subclust_idents = c(ls_div$subclust_idents, ls_in$subclust_idents, ls_glia$subclust_idents, ls_ex$subclust_idents)
all_subclust_idents = data.frame(all_subclust_idents)
all_subclust_idents$CBC = rownames(all_subclust_idents)
colnames(all_subclust_idents) = c("subclust_idents", "CBC")
head(all_subclust_idents)
dim(all_subclust_idents)
dim(seurat_meta)
colnames(seurat_meta)
seurat_meta$CBC = rownames(seurat_meta)
head(left_join(seurat_meta, all_subclust_idents))
seurat_meta = left_join(seurat_meta, all_subclust_idents)
sum(is.na(seurat_meta$subclust_idents))
#[1] 83

ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, seurat_meta$subclust_idents, col.name = "subclust_idents")
table(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents)
#   EX   EX_early       GLIA   GLIA_DIV     IN_CGE   IN_local     IN_MGE     IN_UNK     IPC_EX 
#49657       7235       9344       2303       2860       5445       2621        828       9345 
#IPC_EX_DIV     IPC_IN IPC_IN_DIV       OLIG        OPC    OPC_DIV      RG_EX 
#      2060       1192       2605        755       2249        543       1446


# RG_EX... are these real cells or doublets?
DimPlot(ls_glia, label=T)
FeaturePlot(ls_glia, label=T, features=c("VIM", "HOPX", "NEUROD1", "NEUROD2", "MARCH1", "DLX5", "PTPRZ1", "FAM107A", "MKI67"))
FeaturePlot(ls_glia, features=c("NEUROD2"), split.by="Clone_Index", label=T)
DimPlot(ls_synthesis_v3_harmony_v3, raster=F, cells.highlight=rownames(ls_synthesis_v3_harmony_v3@meta.data[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents %in% "RG_EX",]))
DimPlot(ls_synthesis_v3_harmony_v3, raster=F, cells.highlight=rownames(ls_synthesis_v3_harmony_v3@meta.data[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents %in% "IPC_IN_DIV",]))

FeaturePlot(ls_synthesis_v3_harmony_v3, raster=F, features=c("NEUROD2", "VIM"))
FeaturePlot(ls_glia, label=T, features=c("VIM", "HOPX", "DCX", "NEUROD2"))
# RG_EX look like doublets, not going to consider for further analysis


#### -- intermediate analysis initially looking into specific glia and radial glia identities -- ####

## revisiting glial subclustering
colnames(ls_synthesis_v3_harmony_v3@meta.data)
table(ls_synthesis_v3_harmony_v3@meta.data$clust_idents)
table(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents)

ls_glia_v2 = subset(ls_synthesis_v3_harmony_v3, subclust_idents %in% c("GLIA", "GLIA_DIV"))
ls_glia_v2
#27581 features across 11647 samples within 1 assay 
ls_glia_v2 <- ScaleData(ls_glia_v2, features = all.genes)
ls_glia_v2 = FindVariableFeatures(ls_glia_v2, selection.method="vst", nfeatures=2000)
ls_glia_v2 = RunPCA(ls_glia_v2, features=VariableFeatures(object=ls_glia_v2))
ls_glia_v2 <- FindNeighbors(ls_glia_v2, dims = 1:15)
ls_glia_v2 <- FindClusters(ls_glia_v2, resolution = .5)
ls_glia_v2 <- RunUMAP(ls_glia_v2, dims=c(1:15))
DimPlot(ls_glia_v2, reduction = "umap",label=T)
DimPlot(ls_glia_v2, reduction = "umap",group.by="orig.ident")
DimPlot(ls_glia_v2, reduction = "umap",group.by="age_pseudo")
DimPlot(ls_glia_v2, reduction = "umap",group.by="sample_age")
FeaturePlot(ls_glia_v2, features=c("VIM", "HOPX", "EOMES", "DLX2", "OLIG2", "NEUROD2", "ID3", "PPP1R17", "MKI67"))
FeaturePlot(ls_glia_v2, features=c("CRYAB"), order=T)
FeaturePlot(ls_glia_v2, features=c("HOPX"), order=T)
VlnPlot(ls_glia_v2, features=c("CRYAB"))
VlnPlot(ls_glia_v2, features=c("HOPX"))
FeaturePlot(ls_glia_v2, features=c("VIM"))
FeaturePlot(ls_glia_v2, features=c("AQP4"))
FeaturePlot(ls_glia_v2, features=c("AQP4","GFAP","VIM","HOPX"))

ls_glia_v2_harmony <- RunHarmony(ls_glia_v2, c("orig.ident"))
ls_glia_v2_harmony <- FindNeighbors(ls_glia_v2_harmony, reduction = "harmony", dims = 1:15)
ls_glia_v2_harmony <- FindClusters(ls_glia_v2_harmony, resolution = .5)
ls_glia_v2_harmony <- RunUMAP(ls_glia_v2_harmony, reduction = "harmony", dims=c(1:15))
DimPlot(ls_glia_v2_harmony, reduction = "umap",label=T)
DimPlot(ls_glia_v2_harmony, reduction = "umap",group.by="orig.ident")
DimPlot(ls_glia_v2_harmony, reduction = "umap",group.by="age_pseudo")
DimPlot(ls_glia_v2_harmony, reduction = "umap",group.by="sample_age")
FeaturePlot(ls_glia_v2_harmony, features=c("VIM", "HOPX", "EOMES", "DLX2", "OLIG2", "NEUROD2", "ID3", "PPP1R17", "MKI67"))
FeaturePlot(ls_glia_v2_harmony, features=c("CRYAB"), order=T, label=T)


ls_glia_v2_harmony = SetIdent(ls_glia_v2_harmony, value="sample_age")
VlnPlot(ls_glia_v2_harmony, features="NFIA")
VlnPlot(ls_glia_v2_harmony, features="NFIX")
VlnPlot(ls_glia_v2_harmony, features="NFIB")

ls_glia_v2 = RunPCA(ls_glia_v2, features=VariableFeatures(object=ls_glia_v2), ndims.print=1:15, nfeatures.print = 10)
ls_glia_v2 = CellCycleScoring(ls_glia_v2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(ls_glia_v2, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
ls_glia_v2 = RunPCA(ls_glia_v2, features=c(s.genes, g2m.genes))
DimPlot(ls_glia_v2)
ls_glia_v2 <- ScaleData(ls_glia_v2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ls_glia_v2))
ls_glia_v2 <- RunPCA(ls_glia_v2, features = VariableFeatures(ls_glia_v2), nfeatures.print = 10)
ls_glia_v2 <- RunPCA(ls_glia_v2, features = c(s.genes, g2m.genes))
DimPlot(ls_glia_v2)
head(ls_glia_v2@meta.data)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ls_glia_v2_harmony = CellCycleScoring(ls_glia_v2_harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(ls_glia_v2_harmony[[]])
RidgePlot(ls_glia_v2_harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
ls_glia_v2_harmony = RunPCA(ls_glia_v2_harmony, features=c(s.genes, g2m.genes))
DimPlot(ls_glia_v2_harmony)
ls_glia_v2_harmony <- ScaleData(ls_glia_v2_harmony, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ls_glia_v2_harmony))
ls_glia_v2_harmony <- RunPCA(ls_glia_v2_harmony, features = VariableFeatures(ls_glia_v2_harmony), nfeatures.print = 10)
DimPlot(ls_glia_v2_harmony)
ls_glia_v2_harmony <- FindNeighbors(ls_glia_v2_harmony, reduction = "harmony", dims = 1:15)
ls_glia_v2_harmony <- FindClusters(ls_glia_v2_harmony, resolution = .5)
ls_glia_v2_harmony <- RunUMAP(ls_glia_v2_harmony, reduction = "harmony", dims=c(1:15))
DimPlot(ls_glia_v2_harmony, label=T)
DimPlot(ls_glia_v2_harmony, label=T, group.by="Phase")
#clusters 1, 4, 8, 11 are in G2M or S phase

ls_glia_v2_harmony = SetIdent(ls_glia_v2_harmony, value="seurat_clusters")
ls_glia_v2_harmony.markers <- FindAllMarkers(ls_glia_v2_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(ls_glia_v2_harmony.markers)
ls_glia_v2_harmony_clustermarkers_unfilt = as.data.frame(ls_glia_v2_harmony.markers %>% group_by(cluster))
write.csv(ls_glia_v2_harmony_clustermarkers_unfilt,'/media/chang/HDD-11/mgkeefe/R_work_in_progress_cajal/ls_glia_v2_harmony_clustermarkers_v1_unfilt.csv',row.names=F)
ls_glia_v2_harmony.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
ls_glia_v2_harmony_clustermarkers = as.data.frame(ls_glia_v2_harmony.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
write.csv(ls_glia_v2_harmony_clustermarkers,'/media/chang/HDD-11/mgkeefe/R_work_in_progress_cajal/ls_glia_v2_harmony_clustermarkers_v1.csv',row.names=F)
ls_glia_v2_harmony_clustermarkers

ls_glia_v2_harmony = SetIdent(ls_glia_v2_harmony, value="sample_age")
VlnPlot(subset(ls_glia_v2_harmony, sample_age %in% c(16,19,20,22,23,24)), features=c("EOMES"))
VlnPlot(subset(ls_glia_v2_harmony, sample_age %in% c(16,19,20,22,23,24)), features=c("HOPX"))
VlnPlot(subset(ls_glia_v2_harmony, sample_age %in% c(16,19,20,22,23,24)), features=c("CRYAB"))
VlnPlot(subset(ls_glia_v2_harmony, sample_age %in% c(16,19,20,22,23,24)), features=c("AQP4"))
VlnPlot(subset(ls_glia_v2_harmony, sample_age %in% c(16,19,20,22,23,24)), features=c("GJA1"))
ls_glia_v2_harmony = SetIdent(ls_glia_v2_harmony, value="seurat_clusters")
FeaturePlot(ls_glia_v2_harmony, features=c("GJA1", "ID3", "AQP4", "CRYAB"), label=T)
FeaturePlot(ls_glia_v2_harmony, features=c("CD44", "SPARCL1", "ITGB4", "ANGPTL4"), label=T)
FeaturePlot(ls_glia_v2_harmony, features=c("CDON", "GJA1"), label=T, pt.size=1)
FeaturePlot(ls_glia_v2_harmony, features=c("HES1", "HES5", "HES6", "EOMES"), label=T, pt.size=1)


DimPlot(ls_glia_v2_harmony, group.by="seurat_clusters", label=T)
ls_glia_v2_harmony$glial_subtype = NA
ls_glia_v2_harmony$glial_subtype[ls_glia_v2_harmony$seurat_clusters %in% c(3)] = "S_phase_div_glia"
ls_glia_v2_harmony$glial_subtype[ls_glia_v2_harmony$seurat_clusters %in% c(4,9,6)] = "G2M_phase_div_glia"
ls_glia_v2_harmony$glial_subtype[ls_glia_v2_harmony$seurat_clusters %in% c(2)] = "ASTRO"
ls_glia_v2_harmony$glial_subtype[ls_glia_v2_harmony$seurat_clusters %in% c(1)] = "Transitioning OPCs/IPCs"
ls_glia_v2_harmony$glial_subtype[ls_glia_v2_harmony$seurat_clusters %in% c(8)] = "Transitioning IPCs"
ls_glia_v2_harmony$glial_subtype[ls_glia_v2_harmony$seurat_clusters %in% c(0)] = "RG"
ls_glia_v2_harmony$glial_subtype[ls_glia_v2_harmony$seurat_clusters %in% c(5)] = "RG (putative tRG)"
ls_glia_v2_harmony$glial_subtype[ls_glia_v2_harmony$seurat_clusters %in% c(7)] = "RG (putative oRG)"

DimPlot(ls_glia_v2_harmony, group.by="glial_subtype", label=T, pt.size=0.7)
FeaturePlot(ls_glia_v2_harmony, features=c("PTPRZ1", "HOPX", "FAM107A", "EGFR"))
FeaturePlot(ls_glia_v2_harmony, features=c("VIM", "HOPX", "CRYAB", 
                                           "MKI67", "EOMES", "OLIG2",
                                           "AQP4", "GJA1", "FAM107A"))
ls_glia_v2_harmony$glial_subtype = factor(ls_glia_v2_harmony$glial_subtype, levels = c("S_phase_div_glia", "G2M_phase_div_glia", "RG", "RG (putative tRG)", "RG (putative oRG)", "Transitioning IPCs", "Transitioning OPCs/IPCs", "ASTRO"))
ls_glia_v2_harmony = SetIdent(ls_glia_v2_harmony, value="glial_subtype")


## map the new glial subcluster identitites
DimPlot(ls_synthesis_v3_harmony_v3, raster=F)
table(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents)
table(ls_div$subclust_idents)
#  GLIA_DIV IPC_EX_DIV IPC_IN_DIV    OPC_DIV 
#      2303       2060       2334        543
table(ls_in$subclust_idents)
#    IN_CGE   IN_local     IN_MGE     IN_UNK IPC_IN_DIV 
#      2860       5445       2621        828        271
table(ls_glia_v2_harmony$glial_subtype)
#S_phase_div_glia      G2M_phase_div_glia                      RG       RG (putative tRG) 
#            1359                    2374                    2412                    1068 
#RG (putative oRG)      Transitioning IPCs Transitioning OPCs/IPCs                   ASTRO 
#              802                     725                    1439                    1368
table(ls_ex$subclust_idents)
#      EX EX_early   IPC_EX 
#   49657     7235     8636 

## compile all subclust_idents and add to the main object
all_subclust_idents = data.frame(subclust_idents = ls_glia_v2_harmony$glial_subtype)
all_subclust_idents = rbind(all_subclust_idents, data.frame(subclust_idents = ls_in$subclust_idents))
all_subclust_idents = rbind(all_subclust_idents, data.frame(subclust_idents = ls_ex$subclust_idents))
all_subclust_idents = rbind(all_subclust_idents, data.frame(subclust_idents = ls_div$subclust_idents))
all_subclust_idents$CBC = rownames(all_subclust_idents)
colnames(all_subclust_idents) = c("subclust_idents_v2", "CBC")
head(all_subclust_idents)
dim(all_subclust_idents)
seurat_meta = ls_synthesis_v3_harmony_v3@meta.data
dim(seurat_meta)
length(unique(all_subclust_idents)$CBC)
table(all_subclust_idents$subclust_idents_v2)
table(seurat_meta$subclust_idents)
colnames(seurat_meta)
seurat_meta$CBC = rownames(seurat_meta)
tail(left_join(seurat_meta, all_subclust_idents))
seurat_meta = left_join(seurat_meta, all_subclust_idents)
sum(is.na(seurat_meta$subclust_idents_v2))
seurat_meta$subclust_idents_v2[seurat_meta$subclust_idents %in% 'OPC'] = 'OPC'
seurat_meta$subclust_idents_v2[seurat_meta$subclust_idents %in% 'OLIG'] = 'OLIG'
seurat_meta$subclust_idents_v2[seurat_meta$subclust_idents %in% 'RG_EX'] = 'RG_UNK'
table(seurat_meta$subclust_idents_v2)

rownames(seurat_meta) = seurat_meta$CBC
length(rownames(seurat_meta))
sum(rownames(seurat_meta) == rownames(ls_synthesis_v3_harmony_v3@meta.data))
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, seurat_meta$subclust_idents_v2, col.name = 'subclust_idents_v2')
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents == 'IPC_IN'] = 'IPC_IN'
sum(is.na(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2))
#[1] 2338
table(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents[is.na(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2)])
#GLIA_DIV       GLIA      RG_EX IPC_EX_DIV     IPC_EX   EX_EARLY         EX IPC_IN_DIV     IPC_IN 
#       1         99          0          0        709          0          0          0          0 
#IN_local    OPC_DIV        OPC       OLIG     IN_CGE     IN_MGE     IN_UNK 
#       0          0          0          0          0          0          0
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents == 'IPC_EX'] = 'IPC_EX'
sum(is.na(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2))
#[1] 1629
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2[is.na(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2)]=rep("unk", sum(is.na(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2)))
sum(is.na(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2))
#[1] 0

DimPlot(ls_synthesis_v3_harmony_v3, group.by="subclust_idents_v2", raster=F)
DimPlot(ls_synthesis_v3_harmony_v3, group.by="subclust_idents", raster=F)
table(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents)
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2[is.na(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents)] = 'unk'
DimPlot(ls_synthesis_v3_harmony_v3, group.by="subclust_idents_v2", raster=F)


#going to delete IN_UNK from the object
#these cells are almost certainly doublets and are confusing the analysis
FeaturePlot(ls_in, features="nCount_RNA", max.cutoff="q90")
FeaturePlot(ls_in, features="NEUROD2", max.cutoff="q90")
ls_in = subset(ls_in, subclust_idents %in% c("IN_local", "IN_MGE", "IN_CGE", "IPC_IN_DIV"))
library(harmony)
ls_in_new = RunPCA(ls_in, features=VariableFeatures(object=ls_in))
ls_in_new <- RunHarmony(ls_in_new, c("orig.ident"))
ls_in_new <- FindNeighbors(ls_in_new, reduction = "harmony", dims = 1:15)
ls_in_new <- FindClusters(ls_in_new, resolution = .5)
ls_in_new <- RunUMAP(ls_in_new, reduction = "harmony", dims=c(1:15))
DimPlot(ls_in_new, reduction = "umap",label=T)
#inverting the UMAP_2 embedding to make it look more like the original
ls_in_new[["umap"]]@cell.embeddings[,2] = -ls_in_new[["umap"]]@cell.embeddings[,2]


table(ls_glia_v2_harmony@meta.data$glial_subtype)
ls_glia_v2_harmony@meta.data$glial_subtype[ls_glia_v2_harmony@meta.data$glial_subtype %in% c("Transitioning OPCs/IPCs")] = "Transitioning OPCs"
ls_glia_v2_harmony@meta.data$glial_subtype[ls_glia_v2_harmony@meta.data$glial_subtype %in% c("G2M_phase_div_glia")] = "G2M dividing"
ls_glia_v2_harmony@meta.data$glial_subtype[ls_glia_v2_harmony@meta.data$glial_subtype %in% c("S_phase_div_glia")] = "S dividing"
ls_glia_v2_harmony@meta.data$glial_subtype[ls_glia_v2_harmony@meta.data$glial_subtype %in% c("ASTRO")] = "Astrocytes"
ls_glia_v2_harmony = SetIdent(ls_glia_v2_harmony, value="glial_subtype")
DimPlot(ls_glia_v2_harmony)


##### This is the start of analysis that actually makes it into final objects in the paper #####
ls_synthesis_v3_harmony_v3 = SetIdent(ls_synthesis_v3_harmony_v3, value="subclust_idents_v2")
DimPlot(ls_synthesis_v3_harmony_v3, raster=F, label=T)
DimPlot(ls_synthesis_v3_harmony_v3, raster=F, label=T)
DimPlot(subset(ls_synthesis_v3_harmony_v3, subclust_idents_v2 %in% c("IN_UNK"), invert=T), raster=F, label=T)

#generate IN subset including IN IPCs this time
DimPlot(ls_in)
ls_in_v2 = subset(ls_synthesis_v3_harmony_v3, subclust_idents_v2 %in% c("IN_CGE", "IN_MGE", "IN_local", "IPC_IN"))
ls_in_v2
#27581 features across 12118 samples within 1 assay
ls_in_v2 <- ScaleData(ls_in_v2, features = all.genes)
ls_in_v2 = FindVariableFeatures(ls_in_v2, selection.method="vst", nfeatures=2000)
ls_in_v2 = RunPCA(ls_in_v2, features=VariableFeatures(object=ls_in_v2))
ElbowPlot(ls_in_v2)
ls_in_v2 <- RunHarmony(ls_in_v2, c("orig.ident"))
ls_in_v2 <- FindNeighbors(ls_in_v2, reduction = "harmony", dims = 1:15)
ls_in_v2 <- FindClusters(ls_in_v2, resolution = .5)
ls_in_v2 <- RunUMAP(ls_in_v2, reduction = "harmony", dims=c(1:15))
DimPlot(ls_in_v2, reduction = "umap",label=T, pt.size=1)
VlnPlot(ls_in_v2, features = c("nCount_RNA", "percent.mt"), ncol = 2)
DimPlot(ls_in_v2, reduction = "umap",group.by="orig.ident", pt.size=1)
DimPlot(ls_in_v2, reduction = "umap",group.by="age_pseudo", pt.size=1)
DimPlot(subset(ls_in_v2, subset = Clone_Index %in% c("IndexE", "Index3")), label=T, cells.highlight=rownames(ls_in_v2@meta.data[ls_in_v2@meta.data$"Clone_size_corrected">=2,]), split.by="Clone_Index") +ggtitle("Clone size >1")
DimPlot(subset(ls_in_v2, subset = Clone_Index %in% c("IndexE", "Index3")), group.by="Clone_Index", split.by="age_pseudo", pt.size=1)
FeaturePlot(ls_in_v2, label=T, features=c("GAD2", "NEUROD2", "NR2F1", "PAX6", "MKI67", "NPY", "NR2F2", "ADARB2", "MAF"))
ls_in_v2[["umap"]]@cell.embeddings[,1] = -ls_in_v2[["umap"]]@cell.embeddings[,1]
ls_in_v2[["umap"]]@cell.embeddings[,2] = -ls_in_v2[["umap"]]@cell.embeddings[,2]
FeaturePlot(ls_in_v2, label=T, features=c("GAD2", "NEUROD2", "NR2F1", "PAX6", "MKI67", "NPY", "NR2F2", "ADARB2", "MAF"))
DimPlot(ls_in_v2, group.by="subclust_idents_v2", label=T)
FeaturePlot(ls_in_v2, features=c("GAD2", "NR2F2", "LHX6", "ERBB4", "PAX6", "SOX6"), ncol=3)

FeaturePlot(ls_synthesis_v3_harmony_v3, features=c("MKI67", "DLX2", "PAX6"), raster=F, ncol=3)
DotPlot(ls_synthesis_v3_harmony_v3, features=c("MKI67", "DLX2", "PAX6", "OLIG2", "EOMES", "VIM"))
ls_synthesis_v3_harmony_v3@active.ident = factor(ls_synthesis_v3_harmony_v3@active.ident, levels=c(
  "ASTRO",
  "RG",
  "RG (putative oRG)",
  "RG (putative tRG)",
  "S_phase_div_glia",
  "G2M_phase_div_glia",
  "Transitioning IPCs",
  "IPC_EX",
  "IPC_EX_DIV",
  "IPC_IN_DIV",
  "IPC_IN",
  "IN_local",
  "IN_MGE",
  "IN_CGE",
  "EX",
  "EX_early",
  "OPC_DIV",
  "Transitioning OPCs/IPCs",
  "OPC",
  "OLIG"
))

colnames(ls_in_v2@meta.data)
head(ls_in_v2@meta.data)
ls_in_v2@active.ident = SetIdent(ls_in_v2, value="subclust_idents_v2")
ls_in_v2@active.ident = factor(ls_in_v2@active.ident, levels=c("IPC_IN", "IN_local", "IN_CGE", "IN_MGE"))
VlnPlot(ls_in_v2, features=c("GAD2", "ADARB2", "LHX6", "ERBB4", "NR2F1", "SOX6", "SCGN", "PAX6"), ncol=4, pt.size=0)

real_subclusts = c("ASTRO","RG","RG (putative oRG)","RG (putative tRG)","S_phase_div_glia","G2M_phase_div_glia",
                   "Transitioning IPCs","IPC_EX","IPC_EX_DIV","IPC_IN_DIV","IPC_IN","IN_local","IN_MGE","IN_CGE",
                   "EX","EX_early","OPC_DIV","Transitioning OPCs/IPCs","OPC","OLIG")


## loading old STICR dataset to compare identities
library(Seurat)
setwd("231121_Delgado_STICR_data")
mat = Read10X(".")
meta = read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
old_sticr <- CreateSeuratObject(counts = mat, project = "STICR_v1", meta.data=meta)
old_sticr
#33538 features across 121290 samples within 1 assay 
head(old_sticr@meta.data)
old_sticr_backup = old_sticr

# scale data based on normalization
all.genes <- rownames(old_sticr)
old_sticr <- ScaleData(old_sticr, features = all.genes)
# find variably expressed genes
old_sticr <- FindVariableFeatures(old_sticr, selection.method = "vst", nfeatures = 2000)
# perform dimensionality reduction
old_sticr <- RunPCA(old_sticr, features = VariableFeatures(object = old_sticr))
ElbowPlot(old_sticr)
# cluster cells
old_sticr <- FindNeighbors(old_sticr, dims = 1:15)
old_sticr <- FindClusters(old_sticr, resolution = .5)
# make graph
old_sticr <- RunUMAP(old_sticr, dims = 1:15)
DimPlot(old_sticr, raster=F)
DimPlot(old_sticr, raster=F, group.by="Louvain_Cluster")
FeaturePlot(old_sticr, raster=F, features=c("SATB2", "DLX5", "NR2F1", "LHX6"))

original_umap = read.csv('Seurat_umap.coords.tsv.gz', sep="\t", header=F, row.names=1)
head(original_umap)
colnames(original_umap) = c("UMAP_1", "UMAP_2")
head(old_sticr[["umap"]]@cell.embeddings)
dim(original_umap)
dim(old_sticr[["umap"]]@cell.embeddings)
sum(rownames(original_umap) == rownames(old_sticr[["umap"]]@cell.embeddings))
old_sticr_new_umap = old_sticr[["umap"]]@cell.embeddings
original_umap_mat = as.matrix(original_umap)
head(original_umap_mat)[,1]
old_sticr[["umap"]]@cell.embeddings = original_umap_mat
DimPlot(old_sticr, raster=F)
DimPlot(old_sticr, raster=F, group.by="Louvain_Cluster")
FeaturePlot(old_sticr, raster=F, features=c("SATB2", "DLX5", "TSHZ1", "LHX6"))

#subset to just INs?
old_sticr = AddMetaData(old_sticr, metadata = rownames(old_sticr@meta.data), col.name = "cell_ID")
old_sticr_in_meta = read.csv('231121_Delgado_STICR_data/invitro_in_subcluster/meta.tsv', sep='\t', header=T, row.names=1)
head(old_sticr_in_meta)
old_sticr_in_meta$cell_ID = rownames(old_sticr_in_meta)
old_sticr_in = subset(old_sticr, cell_ID %in% old_sticr_in_meta$cell_ID)
old_sticr_in
#33538 features across 77748 samples within 1 assay 
# scale data based on normalization
all.genes <- rownames(old_sticr_in)
old_sticr_in <- ScaleData(old_sticr_in, features = all.genes)
# find variably expressed genes
old_sticr_in <- FindVariableFeatures(old_sticr_in, selection.method = "vst", nfeatures = 2000)
# perform dimensionality reduction
old_sticr_in <- RunPCA(old_sticr_in, features = VariableFeatures(object = old_sticr_in))
ElbowPlot(old_sticr_in)
# cluster cells
old_sticr_in <- FindNeighbors(old_sticr_in, dims = 1:15)
old_sticr_in <- FindClusters(old_sticr_in, resolution = .5)
# make graph
old_sticr_in <- RunUMAP(old_sticr_in, dims = 1:15)
DimPlot(old_sticr_in)
old_sticr_in_new_umap = old_sticr_in[["umap"]]@cell.embeddings

old_sticr_in_original_umap = read.csv('231121_Delgado_STICR_data/invitro_in_subcluster/Seurat_umap.coords.tsv.gz',
                                      sep="\t", header=F, row.names=1)
colnames(old_sticr_in_original_umap) = c("UMAP_1", "UMAP_2")
old_sticr_in_original_umap = as.matrix(old_sticr_in_original_umap)
old_sticr_in[["umap"]]@cell.embeddings = old_sticr_in_original_umap
DimPlot(old_sticr_in)
old_sticr_in = AddMetaData(old_sticr_in, old_sticr_in_meta$Leiden_Cluster_in_Paper, "IN_leiden_clust")
DimPlot(old_sticr_in, group.by = "IN_leiden_clust")

old_sticr_in = AddMetaData(old_sticr_in, NA, "IN_subtype")
old_sticr_in@meta.data$IN_subtype[old_sticr_in@meta.data$IN_leiden_clust %in% c(1,5,13,22)] = "OB"
old_sticr_in@meta.data$IN_subtype[old_sticr_in@meta.data$IN_leiden_clust %in% c(7,10,11,16,2,3,4,6,12)] = "CGE-like"
old_sticr_in@meta.data$IN_subtype[old_sticr_in@meta.data$IN_leiden_clust %in% c(21)] = "MGE"
old_sticr_in@meta.data$IN_subtype[old_sticr_in@meta.data$IN_leiden_clust %in% c(8,9,14,20,15,19,18)] = "IN_prog"
old_sticr_in@meta.data$IN_subtype[old_sticr_in@meta.data$IN_leiden_clust %in% c(17)] = "Early_IN"
table(old_sticr_in@meta.data$IN_subtype)
#CGE-like Early_IN  IN_prog      MGE       OB 
#.  37648     2940    21318     1115    14727 

#adding a new subclust metadata column for the figure
colnames(ls_synthesis_v3_harmony_v3@meta.data)
table(ls_synthesis_v3_harmony_v3$subclust_idents_v2)
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, ls_synthesis_v3_harmony_v3$subclust_idents_v2, "subclust_idents_paper")
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2 %in% 
                                                             c("RG (putative oRG)", "RG (putative tRG)")] = "RG"
DimPlot(ls_synthesis_v3_harmony_v3, group.by="subclust_idents_paper", raster=F)
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2 %in% 
                                                             c("IPC_EX_DIV", "IPC_IN_DIV", "OPC_DIV", "S_phase_div_glia", "G2M_phase_div_glia")] = "Dividing"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2 %in% 
                                                             c("ASTRO")] = "Astrocytes"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2 %in% 
                                                             c("OLIG")] = "Oligodendrocytes"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2 %in% 
                                                             c("Transitioning IPCs", "Transitioning OPCs/IPCs")] = "RG"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2 %in% 
                                                             c("RG")] = "Radial glia"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% 
                                                             c("RG")] = "Radial glia"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2 %in% 
                                                             c("IPC_EX")] = "Excitatory IPCs"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2 %in% 
                                                             c("IPC_IN")] = "Inhibitory IPCs"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2 %in% 
                                                             c("EX")] = "Excitatory neurons"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% 
                                                             c("Excitatory neurons")] = "EX"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_v2 %in% 
                                                             c("EX_early")] = "EX"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% 
                                                             c("Excitatory IPCs")] = "EX IPCs"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% 
                                                             c("Inhibitory IPCs")] = "IN IPCs"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% 
                                                             c("OPC")] = "OPCs"

DimPlot(subset(ls_synthesis_v3_harmony_v3, subclust_idents_paper %in% c("IN_UNK", "unk"), invert=T), group.by="subclust_idents_paper", raster=F)

ls_synthesis_v3_harmony_v3 = SetIdent(ls_synthesis_v3_harmony_v3, value = "subclust_idents_paper")
table(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper)
#      Astrocytes         Dividing               EX          EX IPCs          IN IPCs           IN_CGE         IN_local 
#1368             8941            56892             9345             1192             2860             5445 
#IN_MGE           IN_UNK Oligodendrocytes             OPCs      Radial glia              unk 
#2621              828              755             2249             6446             1629 
ls_synthesis_v3_harmony_v3@active.ident = factor(ls_synthesis_v3_harmony_v3@active.ident, levels=c(
  "Astrocytes",
  "Radial glia",
  "Dividing",
  "OPCs",
  "Oligodendrocytes",
  "EX IPCs",
  "EX",
  "IN IPCs",
  "IN_local",
  "IN_CGE",
  "IN_MGE",
  "IN_UNK",
  "unk"
))
DimPlot(subset(ls_synthesis_v3_harmony_v3, subclust_idents_paper %in% c("IN_UNK", "unk"), invert=T), raster=F)

table(ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper)
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, NA, "lineage_trajectories_paper")
ls_synthesis_v3_harmony_v3@meta.data$lineage_trajectories_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% c('Astrocytes', 'Radial glia')] = 'RG/Astro'
ls_synthesis_v3_harmony_v3@meta.data$lineage_trajectories_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% c('EX', 'EX IPCs')] = 'EX'
ls_synthesis_v3_harmony_v3@meta.data$lineage_trajectories_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% c('IN IPCs', 'IN_CGE', 'IN_local', 'IN_MGE')] = 'IN'
ls_synthesis_v3_harmony_v3@meta.data$lineage_trajectories_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% c('Oligodendrocytes', 'OPCs')] = 'Oligo'
ls_synthesis_v3_harmony_v3@meta.data$lineage_trajectories_paper[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% c('Dividing')] = 'Dividing'
table(ls_synthesis_v3_harmony_v3@meta.data$lineage_trajectories_paper)
#Dividing       EX       IN    Oligo RG/Astro 
#8941    66237    12118     3004     7814
div_col = "#e091e3"
ex_col = "#f2cb30"
in_col = '#fc7253'
olig_col = '#52ba86'
glia_col = '#9056d6'

DimPlot(subset(ls_synthesis_v3_harmony_v3, lineage_trajectories_paper %in% c(NA), invert=T), 
        group.by = "lineage_trajectories_paper", raster=F, cols = c('#e3b1d9', '#ebda96', '#f7a38f', '#8bad9b', '#b6a3cc')) +
  ggtitle('Lineage trajectories') +
  theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size=30, hjust=0.5))

DimPlot(subset(ls_synthesis_v3_harmony_v3, lineage_trajectories_paper %in% c(NA), invert=T), 
                group.by = "broad_celltype", raster=F, cols = c(ex_col, glia_col, in_col, olig_col)) +
  ggtitle('Broad celltypes') +
  theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size=30, hjust=0.5))

## also need to try to do the upset plot for Fig 1
library(Seurat)
library(dplyr)
library(ggplot2)
DimPlot(ls_in_v2)
ls_in_v2 = AddMetaData(ls_in_v2, NA, "subclust_idents_paper")
colnames(ls_in_v2@meta.data)
table(ls_in_v2@meta.data$subclust_idents_v2)
ls_in_v2@meta.data$subclust_idents_paper[ls_in_v2@meta.data$subclust_idents_v2 %in% "IN_CGE"] = "IN_CGE"
ls_in_v2@meta.data$subclust_idents_paper[ls_in_v2@meta.data$subclust_idents_v2 %in% "IN_MGE"] = "IN_MGE"
ls_in_v2@meta.data$subclust_idents_paper[ls_in_v2@meta.data$subclust_idents_v2 %in% "IN_local"] = "IN_local"
ls_in_v2@meta.data$subclust_idents_paper[ls_in_v2@meta.data$subclust_idents_v2 %in% "IPC_IN"] = "IN IPCs"
ls_in_v2@meta.data$subclust_idents_paper[ls_in_v2@meta.data$seurat_clusters %in% c(8)] = "IN_OB"
ls_in_v2 = SetIdent(ls_in_v2, value="subclust_idents_paper")
DimPlot(ls_in_v2)

## work on glial subtypes
table(ls_glia_v2_harmony@meta.data$glial_subtype)
DimPlot(ls_glia_v2_harmony, label=T, group.by="Phase")
ls_glia_v2_harmony = SetIdent(ls_glia_v2_harmony, value="glial_subtype")
RidgePlot(ls_glia_v2_harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

FeaturePlot(ls_glia_v2_harmony, features=c("TOP2A", "OLIG2", "EOMES", "DLX2"))
FeaturePlot(ls_synthesis_v3_harmony_v3, features=c("VIM","EOMES","AQP4"), raster=F)
DimPlot(ls_glia_v2_harmony)
VlnPlot(ls_glia_v2_harmony, features=c("EOMES"))
VlnPlot(ls_glia_v2_harmony, features=c("OLIG2"))
FeaturePlot(ls_glia_v2_harmony, features=c("OLIG2"))
FeaturePlot(ls_glia_v2_harmony, features=c("HOPX"))
FeaturePlot(ls_glia_v2_harmony, features=c("MT3"))
FeaturePlot(ls_glia_v2_harmony, features=c("GJA1"))
FeaturePlot(ls_glia_v2_harmony, features=c("CRYAB"))

FeaturePlot(ls_glia_v2_harmony, features=c("TOP2A"))
FeaturePlot(ls_glia_v2_harmony, features=c("EOMES"))
colnames(ls_glia_v2_harmony@meta.data)

#try to isolate the EOMES+ from OLIG2+ in the middle clusters
ls_glia_v2_harmony <- FindClusters(ls_glia_v2_harmony, resolution = .6)
ls_glia_v2_harmony <- RunUMAP(ls_glia_v2_harmony, reduction = "harmony", dims=c(1:15))
DimPlot(ls_glia_v2_harmony, reduction = "umap",label=T)
FeaturePlot(ls_glia_v2_harmony, features=c("EOMES"))
FeaturePlot(ls_glia_v2_harmony, features=c("OLIG2"))
VlnPlot(ls_glia_v2_harmony, features=c("EOMES", "OLIG2"))
FeaturePlot(ls_glia_v2_harmony, features=c("nCount_RNA"), max.cutoff = "q90")

ls_glia_v2_harmony <- FindClusters(ls_glia_v2_harmony, resolution = .8)
ls_glia_v2_harmony <- RunUMAP(ls_glia_v2_harmony, reduction = "harmony", dims=c(1:15))
DimPlot(ls_glia_v2_harmony, reduction = "umap",label=T)
FeaturePlot(ls_glia_v2_harmony, features=c("EOMES"))
FeaturePlot(ls_glia_v2_harmony, features=c("OLIG2"))
DimPlot(ls_glia_v2_harmony, group.by="glial_subtype",label=T)
table(ls_glia_v2_harmony@meta.data$seurat_clusters)

table(ls_glia_v2_harmony$glial_subtype)
#Astrocytes       G2M dividing                 RG  RG (putative oRG)  RG (putative tRG)         S dividing 
#1368               2374               2412                802               1068               1359 
#Transitioning IPCs Transitioning OPCs 
#725               1439
sum(ls_glia_v2_harmony$glial_subtype %in% c('Transitioning IPCs', 'Transitioning OPCs'))
#[1] 2164
ls_glia_sub_sub = subset(ls_glia_v2_harmony, glial_subtype %in% c('Transitioning IPCs', 'Transitioning OPCs'))
ls_glia_sub_sub <- ScaleData(ls_glia_sub_sub, features = all.genes)
ls_glia_sub_sub = FindVariableFeatures(ls_glia_sub_sub, selection.method="vst", nfeatures=2000)
ls_glia_sub_sub = RunPCA(ls_glia_sub_sub, features=VariableFeatures(object=ls_glia_sub_sub))
ElbowPlot(ls_glia_sub_sub)
ls_glia_sub_sub <- RunHarmony(ls_glia_sub_sub, c("orig.ident"))
ls_glia_sub_sub <- FindNeighbors(ls_glia_sub_sub, reduction = "harmony", dims = 1:15)
ls_glia_sub_sub <- FindNeighbors(ls_glia_sub_sub, dims = 1:15)
ls_glia_sub_sub <- FindClusters(ls_glia_sub_sub, resolution = .5)
#ls_glia_sub_sub <- RunUMAP(ls_glia_sub_sub, reduction = "harmony", dims=c(1:15))
ls_glia_sub_sub <- RunUMAP(ls_glia_sub_sub, dims=c(1:15))
DimPlot(ls_glia_sub_sub, reduction = "umap",label=T)
FeaturePlot(ls_glia_sub_sub, features=c("EOMES", "OLIG2"), label=T)
colnames(ls_glia_sub_sub@meta.data)
DimPlot(ls_glia_sub_sub, group.by="glial_subtype")
DimPlot(ls_glia_sub_sub, group.by="orig.ident")
table(ls_glia_sub_sub@meta.data$glial_subtype)
#Transitioning IPCs Transitioning OPCs 
#725               1439
ls_glia_sub_sub = AddMetaData(ls_glia_sub_sub, NA, "glial_subtype_v2")
ls_glia_sub_sub@meta.data$glial_subtype[ls_glia_sub_sub@meta.data$seurat_clusters %in% c(0,4,5)] = "Transitioning IPCs"
ls_glia_sub_sub@meta.data$glial_subtype[ls_glia_sub_sub@meta.data$seurat_clusters %in% c(1,2,3,6)] = "Transitioning OPCs"
table(ls_glia_sub_sub@meta.data$glial_subtype)
#Transitioning IPCs Transitioning OPCs 
#941               1223

ls_glia_v2_harmony = AddMetaData(ls_glia_v2_harmony, ls_glia_v2_harmony$glial_subtype, "glial_subtype_v2")
table(ls_glia_v2_harmony@meta.data$glial_subtype_v2)
#Astrocytes       G2M dividing                 RG  RG (putative oRG)  RG (putative tRG)         S dividing 
#1368               2374               2412                802               1068               1359 
#Transitioning IPCs Transitioning OPCs 
#725               1439 
ls_glia_v2_harmony@meta.data$glial_subtype_v2[rownames(ls_glia_v2_harmony@meta.data) %in% rownames(ls_glia_sub_sub@meta.data)[ls_glia_sub_sub@meta.data$glial_subtype %in% "Transitioning IPCs"]] = "Transitioning IPCs"
ls_glia_v2_harmony@meta.data$glial_subtype_v2[rownames(ls_glia_v2_harmony@meta.data) %in% rownames(ls_glia_sub_sub@meta.data)[ls_glia_sub_sub@meta.data$glial_subtype %in% "Transitioning OPCs"]] = "Transitioning OPCs"
table(ls_glia_v2_harmony@meta.data$glial_subtype_v2)
#Astrocytes       G2M dividing                 RG  RG (putative oRG)  RG (putative tRG)         S dividing 
#1368               2374               2412                802               1068               1359 
#Transitioning IPCs Transitioning OPCs 
#941               1223
DimPlot(ls_glia_v2_harmony, group.by="glial_subtype_v2")


ls_glia_v2_harmony = SetIdent(ls_glia_v2_harmony, value="glial_subtype_v2")
DimPlot(ls_glia_v2_harmony, label=F, raster=F)
ls_glia_v2_harmony = AddMetaData(ls_glia_v2_harmony, ls_glia_v2_harmony$glial_subtype_v2, "subclust_idents_paper")
ls_glia_v2_harmony@meta.data$subclust_idents_paper = ls_glia_v2_harmony@meta.data$glial_subtype_v2
ls_glia_v2_harmony@meta.data$subclust_idents_paper[ls_glia_v2_harmony@meta.data$glial_subtype %in% c('G2M dividing')] =
  "Dividing (G2M phase)"
ls_glia_v2_harmony@meta.data$subclust_idents_paper[ls_glia_v2_harmony@meta.data$glial_subtype %in% c('S dividing')] =
  "Dividing (S phase)"
ls_glia_v2_harmony@meta.data$subclust_idents_paper[ls_glia_v2_harmony@meta.data$glial_subtype %in% c('RG (putative tRG)')] =
  "tRG"
ls_glia_v2_harmony@meta.data$subclust_idents_paper[ls_glia_v2_harmony@meta.data$glial_subtype %in% c('RG (putative oRG)')] =
  "oRG"
ls_glia_v2_harmony@meta.data$subclust_idents_paper[ls_glia_v2_harmony@meta.data$seurat_clusters %in% c(13)] =
  NA

ls_glia_v2_harmony = SetIdent(ls_glia_v2_harmony, value="subclust_idents_paper")
DimPlot(subset(ls_glia_v2_harmony, subclust_idents_paper %in% NA, invert=T), raster=F) +
  theme(legend.text = element_text(size=16))

table(ls_glia_v2_harmony@meta.data$subclust_idents_paper)
#Astrocytes Dividing (G2M phase)   Dividing (S phase)                  oRG                   RG   Transitioning IPCs 
#1368                 2374                 1359                  769                 2412                  970 
#Transitioning OPCs                  tRG 
#1194                 1068 

ls_glia_v2_harmony@active.ident = factor(ls_glia_v2_harmony@active.ident, levels=c(
  "Dividing (S phase)",
  "Dividing (G2M phase)",
  "RG", 
  "oRG",
  "tRG",
  "Astrocytes",
  "Transitioning OPCs",
  "Transitioning IPCs"))
ls_glia_v2_harmony$subclust_idents_paper = factor(ls_glia_v2_harmony$subclust_idents_paper, levels = c(
  "Dividing (S phase)",
  "Dividing (G2M phase)",
  "RG", 
  "oRG",
  "tRG",
  "Astrocytes",
  "Transitioning OPCs",
  "Transitioning IPCs"
))

FeaturePlot(ls_glia_v2_harmony, features=c("OLIG2", "EOMES", "VIM", "HOPX"))
##going to try to get better resolution on the dividing glia
table(ls_glia_v2_harmony$subclust_idents_paper)
ls_glia_sub_div_s = subset(ls_glia_v2_harmony, subclust_idents_paper %in% c('Dividing (S phase)'))
ls_glia_sub_div_s <- ScaleData(ls_glia_sub_div_s, features = all.genes)
ls_glia_sub_div_s = FindVariableFeatures(ls_glia_sub_div_s, selection.method="vst", nfeatures=2000)
ls_glia_sub_div_s = RunPCA(ls_glia_sub_div_s, features=VariableFeatures(object=ls_glia_sub_div_s))
ElbowPlot(ls_glia_sub_div_s)
ls_glia_sub_div_s <- RunHarmony(ls_glia_sub_div_s, c("orig.ident"))
ls_glia_sub_div_s <- FindNeighbors(ls_glia_sub_div_s, reduction = "harmony", dims = 1:15)
ls_glia_sub_div_s <- FindNeighbors(ls_glia_sub_div_s, dims = 1:15)
ls_glia_sub_div_s <- FindClusters(ls_glia_sub_div_s, resolution = 0.8)
#ls_glia_sub_div_s <- RunUMAP(ls_glia_sub_div_s, reduction = "harmony", dims=c(1:15))
ls_glia_sub_div_s <- RunUMAP(ls_glia_sub_div_s, dims=c(1:15))
DimPlot(ls_glia_sub_div_s, reduction = "umap",label=T)
FeaturePlot(ls_glia_sub_div_s, features=c("EOMES", "OLIG2", "VIM", "HOPX"), label=T)
DimPlot(ls_glia_sub_div_s, group.by="orig.ident") #very interesting separation by orig.ident even after harmony...
ls_glia_sub_div_s = AddMetaData(ls_glia_sub_div_s, NA, "glia_lineage_trajectory")
ls_glia_sub_div_s@meta.data$glia_lineage_trajectory[ls_glia_sub_div_s@meta.data$seurat_clusters %in% c(7)] = "EX"
ls_glia_sub_div_s@meta.data$glia_lineage_trajectory[ls_glia_sub_div_s@meta.data$seurat_clusters %in% c(1,6)] = "OLIG"
ls_glia_sub_div_s@meta.data$glia_lineage_trajectory[ls_glia_sub_div_s@meta.data$seurat_clusters %in% c(0,2,3,4,5)] = "PROG"
DimPlot(ls_glia_sub_div_s, group.by = "glia_lineage_trajectory")
FeaturePlot(ls_glia_sub_div_s, features=c("EOMES"), label=T, split.by = "glia_lineage_trajectory")
FeaturePlot(ls_glia_sub_div_s, features=c("OLIG2"), label=T, split.by = "glia_lineage_trajectory")
FeaturePlot(ls_glia_sub_div_s, features=c("VIM"), label=T, split.by = "glia_lineage_trajectory")

ls_glia_sub_div_g2m = subset(ls_glia_v2_harmony, subclust_idents_paper %in% c('Dividing (G2M phase)'))
ls_glia_sub_div_g2m <- ScaleData(ls_glia_sub_div_g2m, features = all.genes)
ls_glia_sub_div_g2m = FindVariableFeatures(ls_glia_sub_div_g2m, selection.method="vst", nfeatures=2000)
ls_glia_sub_div_g2m = RunPCA(ls_glia_sub_div_g2m, features=VariableFeatures(object=ls_glia_sub_div_g2m))
ElbowPlot(ls_glia_sub_div_g2m)
ls_glia_sub_div_g2m <- RunHarmony(ls_glia_sub_div_g2m, c("orig.ident"))
ls_glia_sub_div_g2m <- FindNeighbors(ls_glia_sub_div_g2m, reduction = "harmony", dims = 1:15)
ls_glia_sub_div_g2m <- FindNeighbors(ls_glia_sub_div_g2m, dims = 1:15)
ls_glia_sub_div_g2m <- FindClusters(ls_glia_sub_div_g2m, resolution = 0.8)
#ls_glia_sub_div_g2m <- RunUMAP(ls_glia_sub_div_g2m, reduction = "harmony", dims=c(1:15))
ls_glia_sub_div_g2m <- RunUMAP(ls_glia_sub_div_g2m, dims=c(1:15))
DimPlot(ls_glia_sub_div_g2m, reduction = "umap",label=T)
FeaturePlot(ls_glia_sub_div_g2m, features=c("EOMES", "OLIG2", "VIM", "HOPX"), label=T)
DimPlot(ls_glia_sub_div_g2m, group.by="orig.ident") 
ls_glia_sub_div_g2m = AddMetaData(ls_glia_sub_div_g2m, "PROG", "glia_lineage_trajectory")
ls_glia_sub_div_g2m@meta.data$glia_lineage_trajectory[ls_glia_sub_div_g2m@meta.data$seurat_clusters %in% c(1)] = "EX"
ls_glia_sub_div_g2m@meta.data$glia_lineage_trajectory[ls_glia_sub_div_g2m@meta.data$seurat_clusters %in% c(2,9)] = "OLIG"
DimPlot(ls_glia_sub_div_g2m, group.by = "glia_lineage_trajectory")
FeaturePlot(ls_glia_sub_div_g2m, features=c("EOMES"), label=T, split.by = "glia_lineage_trajectory")
FeaturePlot(ls_glia_sub_div_g2m, features=c("OLIG2"), label=T, split.by = "glia_lineage_trajectory")
FeaturePlot(ls_glia_sub_div_g2m, features=c("VIM"), label=T, split.by = "glia_lineage_trajectory")


ls_glia_v2_harmony = AddMetaData(ls_glia_v2_harmony, NA, "glia_lineage_trajectory")
table(ls_glia_v2_harmony$subclust_idents_paper)
ls_glia_v2_harmony$glia_lineage_trajectory[ls_glia_v2_harmony$subclust_idents_paper %in% c("tRG", "oRG")] = "PROG"
ls_glia_v2_harmony$glia_lineage_trajectory[ls_glia_v2_harmony$subclust_idents_paper %in% c("Astrocytes", "RG")] = "ASTRO"
ls_glia_v2_harmony$glia_lineage_trajectory[ls_glia_v2_harmony$subclust_idents_paper %in% c("Transitioning OPCs")] = "OLIG"
ls_glia_v2_harmony$glia_lineage_trajectory[ls_glia_v2_harmony$subclust_idents_paper %in% c("Transitioning IPCs")] = "EX"
length(rownames(ls_glia_v2_harmony@meta.data)[rownames(ls_glia_v2_harmony@meta.data) %in% rownames(ls_glia_sub_div_g2m@meta.data)] == rownames(ls_glia_sub_div_g2m@meta.data))
sum(rownames(ls_glia_v2_harmony@meta.data)[rownames(ls_glia_v2_harmony@meta.data) %in% rownames(ls_glia_sub_div_g2m@meta.data)] == rownames(ls_glia_sub_div_g2m@meta.data))
sum(is.na(ls_glia_v2_harmony@meta.data$glia_lineage_trajectory)[rownames(ls_glia_v2_harmony@meta.data) %in% rownames(ls_glia_sub_div_g2m@meta.data)])
ls_glia_v2_harmony@meta.data$glia_lineage_trajectory[rownames(ls_glia_v2_harmony@meta.data) %in% rownames(ls_glia_sub_div_g2m@meta.data)] = ls_glia_sub_div_g2m@meta.data$glia_lineage_trajectory
ls_glia_v2_harmony@meta.data$glia_lineage_trajectory[rownames(ls_glia_v2_harmony@meta.data) %in% rownames(ls_glia_sub_div_s@meta.data)] = ls_glia_sub_div_s@meta.data$glia_lineage_trajectory
table(ls_glia_v2_harmony@meta.data$glia_lineage_trajectory)
#ASTRO    EX  OLIG  PROG 
# 3780  1303  1995  4436
DimPlot(ls_glia_v2_harmony, group.by="glia_lineage_trajectory")

colnames(ls_glia_sub_div_s@meta.data)
DimPlot(ls_glia_sub_div_s, group.by="glial_subtype")
table(ls_glia_sub_div_s@meta.data$glial_subtype)
#Transitioning IPCs Transitioning OPCs 
#725               1439
ls_glia_sub_div_s = AddMetaData(ls_glia_sub_div_s, NA, "glial_subtype_v2")
ls_glia_sub_div_s@meta.data$glial_subtype[ls_glia_sub_div_s@meta.data$seurat_clusters %in% c(0,4,5)] = "Transitioning IPCs"
ls_glia_sub_div_s@meta.data$glial_subtype[ls_glia_sub_div_s@meta.data$seurat_clusters %in% c(1,3,6)] = "Transitioning OPCs"
table(ls_glia_sub_div_s@meta.data$glial_subtype)
#Transitioning IPCs Transitioning OPCs 
#970               1194

## investigating ls_div to make sure that there are no dividing IN IPCs being labeled as dividing EX IPCs
FeaturePlot(ls_div, features=c("EOMES", "DLX2", "VIM"))
head(ls_div@meta.data)
DimPlot(ls_div, group.by="seurat_clusters")
DimPlot(ls_div, group.by="subclust_idents")
table(ls_div@meta.data$subclust_idents)

ls_div <- FindClusters(ls_div, resolution = .8)
ls_div <- RunUMAP(ls_div, reduction = "harmony", dims=c(1:15))
DimPlot(ls_div, reduction = "umap",label=T)
DimPlot(ls_div, reduction = "umap",group.by="orig.ident")
FeaturePlot(ls_div, features=c("EOMES", "DLX2", "VIM", "OLIG2"), label=T)
ls_div$subclust_idents[ls_div$seurat_clusters %in% c(1,8,7,9,4,6)] = "GLIA_DIV"
ls_div$subclust_idents[ls_div$seurat_clusters %in% c(2,3,5,12,16)] = "IPC_IN_DIV"
ls_div$subclust_idents[ls_div$seurat_clusters %in% c(0,10,11,13,14)] = "EX_IPC_DIV"
DimPlot(ls_div, reduction = "umap",group.by="RNA_snn_res.0.5")
ls_div$subclust_idents[ls_div$RNA_snn_res.0.5 %in% c(7)] = "OPC_DIV"
DimPlot(ls_div, reduction = "umap",group.by="subclust_idents")

## at this point, need to make some final "definitive" cell type calls because this is still affecting everything downstream of it
## going to make some overlapping subclustered objects so that ambiguous cells can be plotted in both and determined

#ls_prog -- glia, EX IPCS, IN IPCS, dividing, OPCs, oligos
#from ls_prog, subset down to ls_prog_div based on cell cycle scores


DimPlot(ls_synthesis_v3_harmony_v3, group.by="subclust_idents_paper", raster=F)
ls_prog = subset(ls_synthesis_v3_harmony_v3, subclust_idents_paper %in% c('Astrocytes', 'Dividing', 'EX IPCs', 'IN IPCs', 'OPCs', 'Radial glia', 'unk', 'Oligodendrocytes'))
ls_prog
#27581 features across 31925 samples within 1 assay 
ls_prog <- ScaleData(ls_prog, features = all.genes)
ls_prog = FindVariableFeatures(ls_prog, selection.method="vst", nfeatures=2000)
ls_prog = RunPCA(ls_prog, features=VariableFeatures(object=ls_prog))
ElbowPlot(ls_prog)
ls_prog <- RunHarmony(ls_prog, c("orig.ident"))
ls_prog <- FindNeighbors(ls_prog, reduction = "harmony", dims = 1:15)
ls_prog <- FindClusters(ls_prog, resolution = .5)
ls_prog <- RunUMAP(ls_prog, reduction = "harmony", dims=c(1:15))
DimPlot(ls_prog, reduction = "umap",label=T)
DimPlot(ls_prog, reduction = "umap", group.by = "orig.ident")
FeaturePlot(ls_prog, features=c("nCount_RNA"), max.cutoff="q90")
#cluster 5 is likely all doublets
FeaturePlot(ls_prog, features=c("MKI67", "TOP2A", "CENPF", "PCNA"), label=T)
FeaturePlot(ls_prog, features=c("VIM", "OLIG2", "EOMES", "DLX2"), label=T)
FeaturePlot(ls_prog, features=c("MBP", "OLIG2", "PLP1", "DLX2"), label=T)
FeaturePlot(ls_prog, features=c("VIM", "GJA1", "HOPX", "CRYAB"), label=T)

DimPlot(ls_prog, group.by="subclust_idents_paper")
ls_prog = AddMetaData(ls_prog, NA, 'subclust_idents_definitive')
ls_prog@meta.data$subclust_idents_definitive[ls_prog@meta.data$seurat_clusters %in% c(1,3,8)] = 'EX IPCs'
ls_prog@meta.data$subclust_idents_definitive[ls_prog@meta.data$seurat_clusters %in% c(12)] = 'Oligodendrocytes'

VlnPlot(ls_prog, features=c("MKI67", "PCNA"))
ls_prog_div = subset(ls_prog, seurat_clusters %in% c(2,4,6,7,9,10))
ls_prog_div <- ScaleData(ls_prog_div, features = all.genes)
ls_prog_div = FindVariableFeatures(ls_prog_div, selection.method="vst", nfeatures=2000)
ls_prog_div = RunPCA(ls_prog_div, features=VariableFeatures(object=ls_prog_div))
ElbowPlot(ls_prog_div)
ls_prog_div <- RunHarmony(ls_prog_div, c("orig.ident"))
ls_prog_div <- FindNeighbors(ls_prog_div, reduction = "harmony", dims = 1:15)
ls_prog_div <- FindClusters(ls_prog_div, resolution = .6)
ls_prog_div <- RunUMAP(ls_prog_div, reduction = "harmony", dims=c(1:15))
DimPlot(ls_prog_div, reduction = "umap",label=T)
DimPlot(ls_prog_div, reduction = "umap", group.by = "orig.ident")
FeaturePlot(ls_prog_div, features=c("nCount_RNA"), max.cutoff="q90")
FeaturePlot(ls_prog_div, features=c("MKI67", "TOP2A", "CENPF", "PCNA"), label=T)
FeaturePlot(ls_prog_div, features=c("VIM", "OLIG2", "EOMES", "DLX2"), label=T)
FeaturePlot(ls_prog_div, features=c("MBP", "OLIG2", "PLP1", "DLX2"), label=T)
FeaturePlot(ls_prog_div, features=c("VIM", "GJA1", "HOPX", "CRYAB"), label=T)


ls_prog_div@meta.data$subclust_idents_definitive[ls_prog_div@meta.data$seurat_clusters %in% c(2,9)] = 'EX IPCs'
ls_prog_div@meta.data$subclust_idents_definitive[ls_prog_div@meta.data$seurat_clusters %in% c(1,3,6,12,4)] = 'IN IPCs'
ls_prog_div@meta.data$subclust_idents_definitive[ls_prog_div@meta.data$seurat_clusters %in% c(11)] = 'OPCs'
ls_prog_div@meta.data$subclust_idents_definitive[ls_prog_div@meta.data$seurat_clusters %in% c(8,5,10)] = 'RG'
#need to split off 4 for further subclustering
#need to split off 0 for further subclustering

ls_prog_div_sub = subset(ls_prog_div, seurat_clusters %in% c(0,4,7))
ls_prog_div_sub
ls_prog_div_sub <- ScaleData(ls_prog_div_sub, features = all.genes)
ls_prog_div_sub = FindVariableFeatures(ls_prog_div_sub, selection.method="vst", nfeatures=2000)
ls_prog_div_sub = RunPCA(ls_prog_div_sub, features=VariableFeatures(object=ls_prog_div_sub))
ElbowPlot(ls_prog_div_sub)
ls_prog_div_sub <- RunHarmony(ls_prog_div_sub, c("orig.ident"))
ls_prog_div_sub <- FindNeighbors(ls_prog_div_sub, reduction = "harmony", dims = 1:15)
ls_prog_div_sub <- FindClusters(ls_prog_div_sub, resolution = .5)
ls_prog_div_sub <- RunUMAP(ls_prog_div_sub, reduction = "harmony", dims=c(1:15))
DimPlot(ls_prog_div_sub, reduction = "umap",label=T)
DimPlot(ls_prog_div_sub, reduction = "umap", group.by = "orig.ident")
FeaturePlot(ls_prog_div_sub, features=c("nCount_RNA"), max.cutoff="q90")
#cluster 5 is likely all doublets
FeaturePlot(ls_prog_div_sub, features=c("MKI67", "TOP2A", "CENPF", "PCNA"), label=T)
FeaturePlot(ls_prog_div_sub, features=c("VIM", "OLIG2", "EOMES", "DLX2"), label=T)
FeaturePlot(ls_prog_div_sub, features=c("MBP", "OLIG2", "PLP1", "DLX2"), label=T)
FeaturePlot(ls_prog_div_sub, features=c("VIM", "OLIG2", "HOPX", "PLP1"), label=T)
VlnPlot(ls_prog_div_sub, features=c("VIM", "OLIG2", "EOMES", "DLX2"))
ls_prog_div_sub@meta.data$subclust_idents_definitive[ls_prog_div_sub@meta.data$seurat_clusters %in% c(1)] = 'IN IPCs'
ls_prog_div_sub@meta.data$subclust_idents_definitive[ls_prog_div_sub@meta.data$seurat_clusters %in% c(2,5,8)] = 'OPCs'
ls_prog_div_sub@meta.data$subclust_idents_definitive[ls_prog_div_sub@meta.data$seurat_clusters %in% c(0,3,4,6,7,9,10)] = 'RG'
sum(is.na(ls_prog_div_sub@meta.data$subclust_idents_definitive))

##start to roll the identities back up the chain
ls_prog_div@meta.data$subclust_idents_definitive[rownames(ls_prog_div@meta.data) %in% rownames(ls_prog_div_sub@meta.data)] = ls_prog_div_sub@meta.data$subclust_idents_definitive
DimPlot(ls_prog_div, group.by="subclust_idents_definitive")
ls_prog_div = SetIdent(ls_prog_div, value="subclust_idents_definitive")
VlnPlot(ls_prog_div, features=c("VIM", "OLIG2", "EOMES", "DLX2"))

ls_prog@meta.data$subclust_idents_definitive[rownames(ls_prog@meta.data) %in% rownames(ls_prog_div@meta.data)] = ls_prog_div@meta.data$subclust_idents_definitive
DimPlot(ls_prog, group.by="subclust_idents_definitive")

#split off 0 and 13 for glial subclustering
ls_prog_glia = subset(ls_prog, seurat_clusters %in% c(0,13))
ls_prog_glia <- ScaleData(ls_prog_glia, features = all.genes)
ls_prog_glia = FindVariableFeatures(ls_prog_glia, selection.method="vst", nfeatures=2000)
ls_prog_glia = RunPCA(ls_prog_glia, features=VariableFeatures(object=ls_prog_glia))
ElbowPlot(ls_prog_glia)
ls_prog_glia <- RunHarmony(ls_prog_glia, c("orig.ident"))
ls_prog_glia <- FindNeighbors(ls_prog_glia, reduction = "harmony", dims = 1:15)
ls_prog_glia <- FindClusters(ls_prog_glia, resolution = .5)
ls_prog_glia <- RunUMAP(ls_prog_glia, reduction = "harmony", dims=c(1:15))
DimPlot(ls_prog_glia, reduction = "umap",label=T)
DimPlot(ls_prog_glia, reduction = "umap", group.by = "orig.ident")
FeaturePlot(ls_prog_glia, features=c("nCount_RNA"), max.cutoff="q90", label=T)
DimPlot(ls_prog_glia, reduction = "umap", group.by = "Clone_Index_filt", pt.size=1)
FeaturePlot(ls_prog_glia, features=c("VIM", "GFAP", "OLIG2", "EOMES"), label=T)
FeaturePlot(ls_prog_glia, features=c("VIM", "HOPX", "CRYAB", "GJA1"), label=T)
ls_prog_glia.markers <- FindAllMarkers(ls_prog_glia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(ls_prog_glia.markers)
ls_prog_glia_clustermarkers_unfilt = as.data.frame(ls_prog_glia.markers %>% group_by(cluster))
write.csv(ls_prog_glia_clustermarkers_unfilt,'/media/chang/HDD-11/mgkeefe/R_work_in_progress_cajal/ls_prog_glia_clustermarkers_v1_unfilt.csv',row.names=F)
ls_prog_glia.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
ls_prog_glia_clustermarkers = as.data.frame(ls_prog_glia.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
write.csv(ls_prog_glia_clustermarkers,'/media/chang/HDD-11/mgkeefe/R_work_in_progress_cajal/ls_prog_glia_clustermarkers_v1.csv',row.names=F)
ls_prog_glia_clustermarkers
FeaturePlot(ls_prog_glia, features=c("HES1", "HES5", "HES6", "NFIA"), label=T)
FeaturePlot(ls_prog_glia, features=c("OLIG2", "PDGFRA", "HOPX", "GAD1"), label=T)
FeaturePlot(ls_prog_glia, features=c("C3", "AIF1", "CCL2", "S100A11"), label=T)
FeaturePlot(ls_prog_glia, features=c("MKI67", "PCNA", "BRIP1", "LHX6"), label=T)

##have to split off 8 and and 9.....
ls_prog_glia_s_phase = subset(ls_prog_glia, seurat_clusters %in% c(8))
ls_prog_glia_s_phase <- ScaleData(ls_prog_glia_s_phase, features = all.genes)
ls_prog_glia_s_phase = FindVariableFeatures(ls_prog_glia_s_phase, selection.method="vst", nfeatures=2000)
ls_prog_glia_s_phase = RunPCA(ls_prog_glia_s_phase, features=VariableFeatures(object=ls_prog_glia_s_phase))
ElbowPlot(ls_prog_glia_s_phase)
ls_prog_glia_s_phase <- RunHarmony(ls_prog_glia_s_phase, c("orig.ident"))
ls_prog_glia_s_phase <- FindNeighbors(ls_prog_glia_s_phase, reduction = "harmony", dims = 1:15)
ls_prog_glia_s_phase <- FindClusters(ls_prog_glia_s_phase, resolution = .5)
ls_prog_glia_s_phase <- RunUMAP(ls_prog_glia_s_phase, reduction = "harmony", dims=c(1:15))
DimPlot(ls_prog_glia_s_phase, reduction = "umap",label=T)
DimPlot(ls_prog_glia_s_phase, reduction = "umap", group.by = "orig.ident")
FeaturePlot(ls_prog_glia_s_phase, features=c("nCount_RNA"), max.cutoff="q90", label=T)
FeaturePlot(ls_prog_glia_s_phase, features=c("VIM", "GAD1", "OLIG2", "HOPX"), label=T)
FeaturePlot(ls_prog_glia_s_phase, features=c("MKI67", "PCNA", "BRIP1", "HOPX"), label=T)
ls_prog_glia_s_phase.markers <- FindAllMarkers(ls_prog_glia_s_phase, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ls_prog_glia_s_phase.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
ls_prog_glia_s_phase_clustermarkers = as.data.frame(ls_prog_glia_s_phase.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
ls_prog_glia_s_phase_clustermarkers
#separation seems to be primarily on the basis of S phase or G2M phase

ls_prog_glia_g2m_phase = subset(ls_prog_glia, seurat_clusters %in% c(9))
ls_prog_glia_g2m_phase <- ScaleData(ls_prog_glia_g2m_phase, features = all.genes)
ls_prog_glia_g2m_phase = FindVariableFeatures(ls_prog_glia_g2m_phase, selection.method="vst", nfeatures=2000)
ls_prog_glia_g2m_phase = RunPCA(ls_prog_glia_g2m_phase, features=VariableFeatures(object=ls_prog_glia_g2m_phase))
ElbowPlot(ls_prog_glia_g2m_phase)
ls_prog_glia_g2m_phase <- RunHarmony(ls_prog_glia_g2m_phase, c("orig.ident"))
ls_prog_glia_g2m_phase <- FindNeighbors(ls_prog_glia_g2m_phase, reduction = "harmony", dims = 1:15)
ls_prog_glia_g2m_phase <- FindClusters(ls_prog_glia_g2m_phase, resolution = .5)
ls_prog_glia_g2m_phase <- RunUMAP(ls_prog_glia_g2m_phase, reduction = "harmony", dims=c(1:15))
DimPlot(ls_prog_glia_g2m_phase, reduction = "umap",label=T)
DimPlot(ls_prog_glia_g2m_phase, reduction = "umap", group.by = "orig.ident")
FeaturePlot(ls_prog_glia_g2m_phase, features=c("nCount_RNA"), max.cutoff="q90", label=T)
FeaturePlot(ls_prog_glia_g2m_phase, features=c("VIM", "HES5", "HES1", "HES6"), label=T)
FeaturePlot(ls_prog_glia_g2m_phase, features=c("MKI67", "PCNA", "BRIP1", "HOPX"), label=T)
DimPlot(ls_prog_glia_g2m_phase, reduction = "umap", group.by = "age_pseudo")
ls_prog_glia_g2m_phase.markers <- FindAllMarkers(ls_prog_glia_g2m_phase, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ls_prog_glia_g2m_phase.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
ls_prog_glia_g2m_phase_clustermarkers = as.data.frame(ls_prog_glia_g2m_phase.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
ls_prog_glia_g2m_phase_clustermarkers
#seems like there may be a bias for ex_neurogenic (AUTS2, TOX, ROBO2, HES6) but not strong enough to claim right now

ls_prog_glia_clust2 = subset(ls_prog_glia, seurat_clusters %in% c(2))
ls_prog_glia_clust2 <- ScaleData(ls_prog_glia_clust2, features = all.genes)
ls_prog_glia_clust2 = FindVariableFeatures(ls_prog_glia_clust2, selection.method="vst", nfeatures=2000)
ls_prog_glia_clust2 = RunPCA(ls_prog_glia_clust2, features=VariableFeatures(object=ls_prog_glia_clust2))
ElbowPlot(ls_prog_glia_clust2)
ls_prog_glia_clust2 <- RunHarmony(ls_prog_glia_clust2, c("orig.ident"))
ls_prog_glia_clust2 <- FindNeighbors(ls_prog_glia_clust2, reduction = "harmony", dims = 1:15)
ls_prog_glia_clust2 <- FindClusters(ls_prog_glia_clust2, resolution = .5)
ls_prog_glia_clust2 <- RunUMAP(ls_prog_glia_clust2, reduction = "harmony", dims=c(1:15))
DimPlot(ls_prog_glia_clust2, reduction = "umap",label=T)
DimPlot(ls_prog_glia_clust2, reduction = "umap", group.by = "orig.ident")
DimPlot(ls_prog_glia_clust2, reduction = "umap", group.by = "Clone_Index_filt")
DimPlot(ls_prog_glia_clust2, reduction = "umap", group.by = "subclust_idents_final")
FeaturePlot(ls_prog_glia_clust2, features=c("nCount_RNA"), max.cutoff="q90", label=T)
FeaturePlot(ls_prog_glia_clust2, features=c("MKI67", "PCNA", "BRIP1", "HOPX"), label=T)
FeaturePlot(ls_prog_glia_clust2, features=c("VIM", "HOPX", "MT3", "PDGFRA"), label=T)
FeaturePlot(ls_prog_glia_clust2, features=c("VIM", "OLIG2", "CLU", "PDGFRA"), label=T)
ls_prog_glia_clust2@meta.data$subclust_idents_definitive[ls_prog_glia_clust2@meta.data$seurat_clusters %in% c(0,2,5,6)] = "oRG"
ls_prog_glia_clust2@meta.data$subclust_idents_definitive[ls_prog_glia_clust2@meta.data$seurat_clusters %in% c(1,3,4)] = "OPCs"


ls_prog_glia_clust6 = subset(ls_prog_glia, seurat_clusters %in% c(6))
ls_prog_glia_clust6 <- ScaleData(ls_prog_glia_clust6, features = all.genes)
ls_prog_glia_clust6 = FindVariableFeatures(ls_prog_glia_clust6, selection.method="vst", nfeatures=2000)
ls_prog_glia_clust6 = RunPCA(ls_prog_glia_clust6, features=VariableFeatures(object=ls_prog_glia_clust6))
ElbowPlot(ls_prog_glia_clust6)
ls_prog_glia_clust6 <- RunHarmony(ls_prog_glia_clust6, c("orig.ident"))
ls_prog_glia_clust6 <- FindNeighbors(ls_prog_glia_clust6, reduction = "harmony", dims = 1:15)
ls_prog_glia_clust6 <- FindClusters(ls_prog_glia_clust6, resolution = .5)
ls_prog_glia_clust6 <- RunUMAP(ls_prog_glia_clust6, reduction = "harmony", dims=c(1:15))
DimPlot(ls_prog_glia_clust6, reduction = "umap",label=T)
DimPlot(ls_prog_glia_clust6, reduction = "umap", group.by = "orig.ident")
DimPlot(ls_prog_glia_clust6, reduction = "umap", group.by = "Clone_Index_filt")
DimPlot(ls_prog_glia_clust6, reduction = "umap", group.by = "subclust_idents_final")
FeaturePlot(ls_prog_glia_clust6, features=c("nCount_RNA"), max.cutoff="q90", label=T)
FeaturePlot(ls_prog_glia_clust6, features=c("MKI67", "PCNA", "BRIP1", "HOPX"), label=T)
FeaturePlot(ls_prog_glia_clust6, features=c("ANXA2", "CRYAB", "C3", "CTSS"), label=T)
FeaturePlot(ls_prog_glia_clust6, features=c("VIM", "HOPX", "GFAP", "AQP4"), label=T)
FeaturePlot(ls_prog_glia_clust6, features=c("VIM", "HOPX", "CLU", "MT3"), label=T)
ls_prog_glia_clust6@meta.data$subclust_idents_definitive[ls_prog_glia_clust6@meta.data$seurat_clusters %in% c(0,1,7)] = "RG"
ls_prog_glia_clust6@meta.data$subclust_idents_definitive[ls_prog_glia_clust6@meta.data$seurat_clusters %in% c(2,3,4,6)] = "tRG"

ls_prog_glia@meta.data$subclust_idents_definitive[rownames(ls_prog_glia@meta.data) %in% rownames(ls_prog_glia_clust2@meta.data)] = ls_prog_glia_clust2@meta.data$subclust_idents_definitive
ls_prog_glia@meta.data$subclust_idents_definitive[rownames(ls_prog_glia@meta.data) %in% rownames(ls_prog_glia_clust6@meta.data)] = ls_prog_glia_clust6@meta.data$subclust_idents_definitive


DimPlot(ls_prog_glia, group.by = "subclust_idents_final", pt.size=1)
DimPlot(ls_prog_glia, group.by = "subclust_idents_final", pt.size=1, split.by='seurat_clusters')
ggplot(ls_prog_glia@meta.data, aes(x=seurat_clusters, fill=subclust_idents_final)) +
  geom_bar()
ls_prog_glia@meta.data$subclust_idents_definitive[ls_prog_glia@meta.data$seurat_clusters %in% c(0, 1, 5, 7, 8, 9)] = "RG"
ls_prog_glia@meta.data$subclust_idents_definitive[ls_prog_glia@meta.data$seurat_clusters %in% c(4)] = "tRG"
ls_prog_glia@meta.data$subclust_idents_definitive[ls_prog_glia@meta.data$seurat_clusters %in% c(10)] = "EX IPCs"
ls_prog_glia@meta.data$subclust_idents_definitive[ls_prog_glia@meta.data$seurat_clusters %in% c(3)] = "Astrocytes"
DimPlot(ls_prog_glia, label=T)
DimPlot(ls_prog_glia, group.by = "subclust_idents_definitive", pt.size=1)


ggplot(ls_prog_glia@meta.data, aes(x=subclust_idents_definitive, fill=Clone_Index_filt)) +
  geom_bar()

ls_prog@meta.data$subclust_idents_definitive[rownames(ls_prog@meta.data) %in% rownames(ls_prog_glia@meta.data)] = ls_prog_glia@meta.data$subclust_idents_definitive

DimPlot(ls_prog, group.by="subclust_idents_definitive", pt.size=1)
DimPlot(ls_prog, label=T)

ls_prog_subclust_11 = subset(ls_prog, seurat_clusters %in% c(11))
ls_prog_subclust_11 <- ScaleData(ls_prog_subclust_11, features = all.genes)
ls_prog_subclust_11 = FindVariableFeatures(ls_prog_subclust_11, selection.method="vst", nfeatures=2000)
ls_prog_subclust_11 = RunPCA(ls_prog_subclust_11, features=VariableFeatures(object=ls_prog_subclust_11))
ElbowPlot(ls_prog_subclust_11)
ls_prog_subclust_11 <- RunHarmony(ls_prog_subclust_11, c("orig.ident"))
ls_prog_subclust_11 <- FindNeighbors(ls_prog_subclust_11, reduction = "harmony", dims = 1:15)
ls_prog_subclust_11 <- FindClusters(ls_prog_subclust_11, resolution = .5)
ls_prog_subclust_11 <- RunUMAP(ls_prog_subclust_11, reduction = "harmony", dims=c(1:15))
DimPlot(ls_prog_subclust_11, reduction = "umap",label=T)
DimPlot(ls_prog_subclust_11, reduction = "umap", group.by = "orig.ident")
DimPlot(ls_prog_subclust_11, reduction = "umap", group.by = "Clone_Index_filt")
DimPlot(ls_prog_subclust_11, reduction = "umap", group.by = "subclust_idents_final")
FeaturePlot(ls_prog_subclust_11, features=c("nCount_RNA"), max.cutoff="q90", label=T)
FeaturePlot(ls_prog_subclust_11, features=c("OLIG2", "PDGFRA", "GAD1", "NPY"), label=T)
#seem like pretty much all OPCs
ls_prog@meta.data$subclust_idents_definitive[ls_prog@meta.data$seurat_clusters %in% c(11)] = "OPCs"
DimPlot(ls_prog, group.by="subclust_idents_definitive", pt.size=1)


DimPlot(ls_synthesis_v3_harmony_v3, label=T, raster=F)
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, NA, "subclust_idents_definitive")
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% c('EX')] = "EX"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_paper %in% c('Oligodendrocytes')] = "Oligodendrocytes"
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive[rownames(ls_synthesis_v3_harmony_v3@meta.data) %in% rownames(ls_in_v2@meta.data)] = ls_in_v2@meta.data$subclust_idents_paper
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive[rownames(ls_synthesis_v3_harmony_v3@meta.data) %in% rownames(ls_prog@meta.data)] = ls_prog@meta.data$subclust_idents_definitive
DimPlot(ls_synthesis_v3_harmony_v3, group.by="subclust_idents_definitive", label=T, raster=F)
DimPlot(ls_synthesis_v3_harmony_v3, group.by="subclust_idents_definitive", raster=F)
## going to try to use these identities, but first going to change some of the nomenclature 
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive %in% 'EX'] = 'ENs'
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive %in% 'EX IPCs'] = 'EX_IPCs'
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive[ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive %in% 'IN IPCs'] = 'IN_IPCs'
DimPlot(ls_synthesis_v3_harmony_v3, group.by="subclust_idents_definitive", raster=F)


## rechecking subtype calls for ls_prog_glia
DimPlot(ls_prog_glia, label=T)
FeaturePlot(ls_prog_glia, features=c("MKI67", "PCNA", "PDGFRA", "EOMES"))
FeaturePlot(ls_prog_glia, features=c("OLIG2", "OLIG1", "PDGFRA", "GSX2"), label=T)
VlnPlot(ls_prog_glia, features=c("OLIG2", "PDGFRA"))
DimPlot(ls_prog_glia_clust2)
FeaturePlot(ls_prog_glia_clust2, features=c("PDGFRA", "OLIG2", "GSX2", "GFAP"))
DimPlot(ls_prog_glia, group.by="subclust_idents_definitive", label=T)


DimPlot(ls_prog_glia, group.by="subclust_idents_final")
DimPlot(ls_prog_glia_clust2, group.by="subclust_idents_final")
DimPlot(ls_prog, cells.highlight = rownames(ls_prog@meta.data)[rownames(ls_prog@meta.data) %in% rownames(ls_prog_glia_clust2@meta.data)])
# call everything in ls_prog_glia cluster2 as OPCs


#0 (or at least a subset of 0???) should definitely be called as oRG based on HOPX, GPX3, (some) MT3
#INPP1 is up in 0, looks pretty good in TN data
#cluster 5 is likely also oRG, just a lower-quality cluster, low RNA counts and RPL and RPS genes
#4 is definitely tRG based on CRYAB

FeaturePlot(ls_prog_glia, features=c("HES1", "VIM", "PAX6", "SPARCL1", "APOE", "PDGFRA", "INPP1", "GPX3", "CRYAB"), ncol=3, pt.size=0.5

ls_prog_glia = AddMetaData(ls_prog_glia, ls_prog_glia@meta.data$subclust_idents_definitive_final, col.name="subclust_names_figure")
ls_prog_glia@meta.data$subclust_names_figure[ls_prog_glia$seurat_clusters %in% c(7)] = "oRG"
ls_prog_glia@meta.data$subclust_names_figure[ls_prog_glia$seurat_clusters %in% c(9)] = "Dividing"
ls_prog_glia@meta.data$subclust_names_figure[ls_prog_glia@meta.data$subclust_names_figure %in% c("tRG")] = "tRG"
ls_prog_glia@meta.data$subclust_names_figure[ls_prog_glia@meta.data$subclust_names_figure %in% c("OPCs")] = "Early OPCs"


DimPlot(ls_prog_glia_clust6, label=T)
FeaturePlot(ls_prog_glia_clust6, features=c("AQP4", "CRYAB", "ITGB4", "SPARCL1", "C3", "IL32"), ncol=3)
FeaturePlot(ls_prog_glia_clust6, features=c("AQP4", "CRYAB", "ZEB1", "SPARCL1", "HOPX", "VIM"), ncol=3)
DimPlot(ls_prog_glia,cells.highlight = rownames(ls_prog_glia@meta.data)[rownames(ls_prog_glia@meta.data) %in% rownames(ls_prog_glia_clust6@meta.data)[ls_prog_glia_clust6@meta.data$seurat_clusters %in% c(0)]])
# ls_prog_glia_clust6 actually look like astrocytes (based on AQP4 and SPARCL1 and following DA's definition of VZ-derived astroyctes)

ls_prog_glia_clust238 = subset(ls_prog_glia, seurat_clusters %in% c(2, 3, 8))
ls_prog_glia_clust238 <- ScaleData(ls_prog_glia_clust238, features = all.genes)
ls_prog_glia_clust238 = FindVariableFeatures(ls_prog_glia_clust238, selection.method="vst", nfeatures=2000)
ls_prog_glia_clust238 = RunPCA(ls_prog_glia_clust238, features=VariableFeatures(object=ls_prog_glia_clust238))
ElbowPlot(ls_prog_glia_clust238)
ls_prog_glia_clust238 <- RunHarmony(ls_prog_glia_clust238, c("orig.ident"))
ls_prog_glia_clust238 <- FindNeighbors(ls_prog_glia_clust238, reduction = "harmony", dims = 1:15)
ls_prog_glia_clust238 <- FindClusters(ls_prog_glia_clust238, resolution = .5)
ls_prog_glia_clust238 <- RunUMAP(ls_prog_glia_clust238, reduction = "harmony", dims=c(1:15))
DimPlot(ls_prog_glia_clust238, reduction = "umap",label=T)
DimPlot(ls_prog_glia_clust238, reduction = "umap", group.by = "orig.ident")
DimPlot(ls_prog_glia_clust238, reduction = "umap", group.by = "Clone_Index_filt")
DimPlot(ls_prog_glia_clust238, reduction = "umap", group.by = "subclust_idents_final")
FeaturePlot(ls_prog_glia_clust238, features=c("nCount_RNA"), max.cutoff="q90", label=T)
FeaturePlot(ls_prog_glia_clust238, features=c("MKI67", "PCNA", "BRIP1", "HOPX"), label=T)
FeaturePlot(ls_prog_glia_clust238, features=c("VIM", "GSX2", "OLIG2", "PDGFRA"), label=T)
FeaturePlot(ls_prog_glia_clust238, features=c("GJA1", "SPARCL1", "ANGPTL4", "GPX3"), label=T)
#0 and 1 are astrocytes based on ANGPTL4, SPARCL1
#5 is OPCs based on OLIG2 and GSX2

DimPlot(ls_prog_glia, cells.highlight = rownames(ls_prog_glia@meta.data)[rownames(ls_prog_glia@meta.data) %in% rownames(ls_prog_glia_clust238@meta.data)[ls_prog_glia_clust238@meta.data$seurat_clusters %in% c(5)]])
DimPlot(ls_prog_glia, cells.highlight = rownames(ls_prog_glia@meta.data)[rownames(ls_prog_glia@meta.data) %in% rownames(ls_prog_glia_clust238@meta.data)[ls_prog_glia_clust238@meta.data$seurat_clusters %in% c(3,5)]])
DimPlot(ls_prog_glia, cells.highlight = rownames(ls_prog_glia@meta.data)[rownames(ls_prog_glia@meta.data) %in% rownames(ls_prog_glia_clust238@meta.data)[ls_prog_glia_clust238@meta.data$seurat_clusters %in% c(3,5)]])

ls_prog_glia_clust238 = AddMetaData(ls_prog_glia_clust238, NA, col.name = "subclust_idents_definitive_final")
ls_prog_glia_clust238@meta.data$subclust_idents_definitive_final[ls_prog_glia_clust238@meta.data$seurat_clusters %in% c(0)] = "Astrocytes"
ls_prog_glia_clust238@meta.data$subclust_idents_definitive_final[ls_prog_glia_clust238@meta.data$seurat_clusters %in% c(1,2,3,5)] = "OPCs"
ls_prog_glia_clust238@meta.data$subclust_idents_definitive_final[ls_prog_glia_clust238@meta.data$seurat_clusters %in% c(4,6,7)] = "oRG"

DimPlot(ls_prog_glia_clust238, group.by="subclust_idents_definitive_final")

## going to add a new metadata column to ls_prog_glia and eventually to everything....... but only going to change ls_prog_glia, will keep everything else the same
ls_prog_glia = AddMetaData(ls_prog_glia, ls_prog_glia@meta.data$subclust_idents_definitive, col.name = "subclust_idents_definitive_final")
ls_prog_glia@meta.data$subclust_idents_definitive_final[rownames(ls_prog_glia@meta.data) %in% rownames(ls_prog_glia_clust238@meta.data)] = ls_prog_glia_clust238@meta.data$subclust_idents_definitive_final
ls_prog_glia@meta.data$subclust_idents_definitive_final[ls_prog_glia@meta.data$seurat_clusters %in% c(6)] = "Astrocytes"
DimPlot(ls_prog_glia, group.by="subclust_idents_definitive_final", pt.size=1)
#after final meeting and discussion, going to use the following definitions:
#0, 1, 5, 7 are oRG
#4 is tRG
#3 and 6 are astrocytes
#2 is OPCs 
#8 is mostly OPCs and some oRG
#9 is dividing
#10 is IPCs
ls_prog_glia@meta.data$subclust_idents_definitive_final[ls_prog_glia@meta.data$seurat_clusters %in% c(0, 1, 5, 7)] = "oRG"
ls_prog_glia@meta.data$subclust_idents_definitive_final[ls_prog_glia@meta.data$seurat_clusters %in% c(2)] = "OPCs"
ls_prog_glia@meta.data$subclust_idents_definitive_final[ls_prog_glia@meta.data$seurat_clusters %in% c(3,6)] = "Astrocytes"
ls_prog_glia@meta.data$subclust_idents_definitive_final[ls_prog_glia@meta.data$seurat_clusters %in% c(4)] = "tRG"
ls_prog_glia@meta.data$subclust_idents_definitive_final[ls_prog_glia@meta.data$seurat_clusters %in% c(9)] = "Dividing"
ls_prog_glia@meta.data$subclust_idents_definitive_final[ls_prog_glia@meta.data$seurat_clusters %in% c(10)] = "EX IPCs"
DimPlot(ls_prog_glia, group.by="subclust_idents_definitive_final", pt.size=1)


## what about Denise's VZ/OSVZ astrocyte markers?
ls_prog_glia = SetIdent(ls_prog_glia, value="seurat_clusters")
FeaturePlot(ls_prog_glia, features=c("ITGB4", "ANGPTL4"))
VlnPlot(ls_prog_glia, features=c("ITGB4", "ANGPTL4"))
FeaturePlot(ls_prog_glia, features=c("FLCN", "NTRK2", "MASP1", "C2orf72", "TIMP3", "ISLR2"), ncol=3)

astro.de.markers = FindMarkers(ls_prog_glia, ident.1 = 3, ident.2 = 6)
head(astro.de.markers, n=20)
ggplot(data=astro.de.markers, aes(x=avg_log2FC, y=p_val_adj))+geom_point()
ggplot(data=astro.de.markers, aes(x=avg_log2FC, y=-log10(p_val_adj), label=rownames(astro.de.markers)))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  ggtitle("Differential expression between Cluster 6 (left) and Cluster 3 (right) -- Astroycytes")+
  theme(plot.title = element_text(hjust = 0.5))

VlnPlot(subset(ls_prog_glia, seurat_clusters %in% c(3,6)), features=c("ITGB4", "ANGPTL4"))
DotPlot(subset(ls_prog_glia, seurat_clusters %in% c(3,6)), features=c("ITGB4", "ANGPTL4"))
DotPlot(ls_prog_glia, features=c("ITGB4", "ANGPTL4"))
DotPlot(subset(ls_prog_glia, seurat_clusters %in% c(3,6)), features=c("ITGB4", "TIMP1", "ANGPTL4",  "TIMP3"))
VlnPlot(subset(ls_prog_glia, seurat_clusters %in% c(3,6)), features=c("ITGB4", "TIMP1", "ANGPTL4",  "TIMP3"))
VlnPlot(subset(ls_prog_glia, seurat_clusters %in% c(3,6)), features=c("ITGB4", "TIMP1", "ANGPTL4",  "TIMP3"))
VlnPlot(subset(ls_prog_glia, seurat_clusters %in% c(3,6)), features=c("ANGPTL2", "DPYD", "TIMP1", "TIMP3"))

ls_prog_glia = SetIdent(ls_prog_glia, value = "seurat_clusters")
DotPlot(subset(ls_prog_glia, seurat_clusters %in% c(3,6)), features=c("ITGB4", "TIMP1", "S100A11", "CRYAB", "NES", "EMP3", "MRC2", "ORMDL2", "ANGPTL4",  "TIMP3", "FLCN", "NTRK2", "NPEPL1", "MASP1", "MPP6"))
ls_prog_glia = SetIdent(ls_prog_glia, value = "Clone_Index_filt")
DotPlot(subset(ls_prog_glia, seurat_clusters %in% c(3,6) & Clone_Index_filt %in% c("IndexE", "Index3")), features=c("ITGB4", "TIMP1", "S100A11", "CRYAB", "NES", "EMP3", "MRC2", "ORMDL2", "ANGPTL4",  "TIMP3", "FLCN", "NTRK2", "NPEPL1", "MASP1", "MPP6"))
VlnPlot(subset(ls_prog_glia, seurat_clusters %in% c(3,6) & Clone_Index_filt %in% c("IndexE", "Index3")), features=c("ITGB4", "TIMP1", "S100A11", "CRYAB", "NES", "EMP3", "MRC2", "ORMDL2", "ANGPTL4",  "TIMP3", "FLCN", "NTRK2", "NPEPL1", "MASP1", "MPP6"))
DoHeatmap(subset(ls_prog_glia, seurat_clusters %in% c(3,6) & Clone_Index_filt %in% c("IndexE", "Index3")), features=c("ITGB4", "TIMP1", "S100A11", "CRYAB", "NES", "EMP3", "MRC2", "ORMDL2", "ANGPTL4",  "TIMP3", "FLCN", "NTRK2", "NPEPL1", "MASP1", "MPP6")) +
  scale_fill_gradientn(colors = c("white", "blue"))
DoHeatmap(subset(ls_prog_glia, seurat_clusters %in% c(3,6) & Clone_Index_filt %in% c("IndexE", "Index3")), features=c("ITGB4", "TIMP1", "S100A11", "CRYAB", "EMP3", "MRC2", "ANGPTL4",  "TIMP3", "FLCN", "NTRK2", "NPEPL1", "MASP1")) +
  scale_fill_gradientn(colors = c("white", "blue"))

ls_prog_glia = SetIdent(ls_prog_glia, value = "seurat_clusters")
DoHeatmap(subset(ls_prog_glia, seurat_clusters %in% c(3,6) & Clone_Index_filt %in% c("IndexE", "Index3")), features=c("ITGB4", "TIMP1", "S100A11", "CRYAB", "EMP3", "MRC2", "ANGPTL4",  "TIMP3", "FLCN", "NTRK2", "NPEPL1", "MASP1")) +
  scale_fill_gradientn(colors = c("white", "blue"))

FeaturePlot(subset(ls_prog_glia, seurat_clusters %in% c(3,6)), features=c("ITGB4", "S100A11", "ANGPTL4", "NTRK2"))


DotPlot(subset(ls_prog_glia, seurat_clusters %in% c(3,6)), features=c("ITGB4", "TIMP1", "S100A11", "CRYAB", "EMP3", "MRC2", "ANGPTL4",  "TIMP3", "NTRK2", "NPEPL1", "MASP1")) +
  ggtitle('Astrocyte gene expression')+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1)
  )


## have to roll the identities all the way up the chain now
ls_prog = AddMetaData(ls_prog, ls_prog@meta.data$subclust_idents_definitive, col.name = "subclust_idents_definitive_final")
ls_prog@meta.data$subclust_idents_definitive_final[rownames(ls_prog@meta.data) %in% rownames(ls_prog_glia@meta.data)] = ls_prog_glia@meta.data$subclust_idents_definitive_final
DimPlot(ls_prog, group.by="subclust_idents_definitive_final", pt.size=1)

ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive, col.name = "subclust_idents_definitive_final")
ls_synthesis_v3_harmony_v3@meta.data$subclust_idents_definitive_final[rownames(ls_synthesis_v3_harmony_v3@meta.data) %in% rownames(ls_prog@meta.data)] = ls_prog@meta.data$subclust_idents_definitive_final
DimPlot(ls_synthesis_v3_harmony_v3, group.by="subclust_idents_definitive_final", pt.size=1, raster=F)


## subclust_idents_definitive_final has updated calls for Astrocytes, OPCs, and RG. RG are split into tRG and oRG here.
## the tRG and oRG cells should be grouped as RG for most figures
## final update, split out the astrocyte subtypes and remake figures
colnames(ls_prog_glia@meta.data)
ls_prog_glia = AddMetaData(ls_prog_glia, ls_prog_glia@meta.data$subclust_idents_definitive_final, col.name = 'subclust_idents_revision')
ls_prog_glia@meta.data$subclust_idents_revision[ls_prog_glia@meta.data$seurat_clusters %in% c(3)] = "Astrocytes (dense smooth)"
ls_prog_glia@meta.data$subclust_idents_revision[ls_prog_glia@meta.data$seurat_clusters %in% c(6)] = "Astrocytes (dense bulbous)"

ls_prog_glia = SetIdent(ls_prog_glia, value = "subclust_idents_revision")
DimPlot(ls_prog_glia)


## final step: filter to cells that have a definitive subcluster identity or are in the IN or prog_glia sub objects
ls_synthesis_v3_harmony_drop_na = subset(ls_synthesis_v3_harmony_v3, subclust_idents_definitive_final %in% NA, invert=T)
ls_synthesis_in = ls_in_v2
ls_synthesis_prog = ls_prog_glia

ls_synthesis_v3_harmony_v3_meta = ls_synthesis_v3_harmony_v3@meta.data
ls_synthesis_v3_harmony_v3_meta = ls_synthesis_v3_harmony_v3_meta[rownames(ls_synthesis_v3_harmony_v3_meta) %in% rownames(ls_synthesis_v3_harmony_drop_na@meta.data) |
                                                                    rownames(ls_synthesis_v3_harmony_v3_meta) %in% rownames(ls_synthesis_in@meta.data) |
                                                                    rownames(ls_synthesis_v3_harmony_v3_meta) %in% rownames(ls_synthesis_prog@meta.data),]
dim(ls_synthesis_v3_harmony_v3_meta)
# [1] 97540    37

## pull subcluster identities up from the subobjects
ls_synthesis_meta_final = ls_synthesis_v3_harmony_v3_meta
ls_in_meta_final = ls_synthesis_in@meta.data
ls_prog_meta_final = ls_synthesis_prog@meta.data

ls_synthesis_meta_final$subcluster_identity = ls_synthesis_meta_final$subclust_idents_definitive_final
ls_in_meta_final$subcluster_identity = ls_in_meta_final$subclust_idents_paper
ls_prog_meta_final$subcluster_identity = ls_prog_meta_final$subclust_idents_revision

ls_synthesis_meta_final$cell = rownames(ls_synthesis_meta_final)
ls_in_meta_final$cell = rownames(ls_in_meta_final)
ls_prog_meta_final$cell = rownames(ls_prog_meta_final)

#remove intersecting cells from the IN and prog object
ls_in_cell_list = ls_in_meta_final$cell
ls_prog_meta_final = ls_prog_meta_final[!(rownames(ls_prog_meta_final) %in% ls_in_cell_list),]
dim(ls_prog_meta_final)
# [1] 6783   15
ls_in_meta_final$subclust_ident_lock = ls_in_meta_final$subcluster_identity
ls_prog_meta_final$subclust_ident_lock = ls_prog_meta_final$subcluster_identity
## for final subcluster idents, pull identities from ls_prog_meta and ls_in_meta
ls_in_meta_final$subclust_idents_from_in = ls_in_meta_final$subclust_ident_lock
ls_prog_meta_final$subclust_idents_from_prog = ls_prog_meta_final$subclust_ident_lock

ls_synthesis_meta_final = left_join(ls_synthesis_meta_final, ls_prog_meta_final[,c("cell", "subclust_idents_from_prog")], by = "cell")
ls_synthesis_meta_final = left_join(ls_synthesis_meta_final, ls_in_meta_final[,c("cell", "subclust_idents_from_in")], by = "cell")

table(ls_prog_meta_final$subclust_idents_from_prog)
# tRG                    EX IPCs Astrocytes (dense bulbous)                        oRG                   Dividing 
# 602                         52                        536                       3410                        218 
# Early OPCs  Astrocytes (dense smooth) 
# 1179                        786
table(ls_synthesis_meta_final$subclust_idents_from_prog)
#tRG                    EX IPCs Astrocytes (dense bulbous)                        oRG                   Dividing 
# 602                         52                        536                       3410                        218 
# Early OPCs  Astrocytes (dense smooth) 
# 1179                        786
table(ls_synthesis_meta_final$subclust_idents_from_in)
#  IN IPCs   IN_CGE IN_local   IN_MGE    IN_OB 
#     1192     2369     5393     2621      543

ls_synthesis_meta_final$subclust_ident_lock <- coalesce(
  ls_synthesis_meta_final$subclust_idents_from_in,
  ls_synthesis_meta_final$subclust_idents_from_prog,
  ls_synthesis_meta_final$subcluster_identity
)
tail(ls_synthesis_meta_final)
table(ls_synthesis_meta_final$subclust_ident_lock)
#Astrocytes (dense bulbous)  Astrocytes (dense smooth)                        oRG                   Dividing                 Early OPCs 
#                       536                        786                       3410                        218                       1179 
#   ENs                    EX IPCs                    EX_IPCs                    IN IPCs                     IN_CGE 
# 56582                         52                      10904                       1192                       2373 
# IN_IPCs                   IN_local                     IN_MGE                      IN_OB           Oligodendrocytes 
#    3403                       5402                       2624                        545                        772 
# OPCs                         RG                        tRG 
# 2431                       4529                        602
ls_synthesis_meta_final$subclust_ident_lock[ls_synthesis_meta_final$subclust_ident_lock == "EX IPCs"] = "EX_IPCs"
ls_synthesis_meta_final$subclust_ident_lock[ls_synthesis_meta_final$subclust_ident_lock == "IN IPCs"] = "IN_IPCs"
table(ls_synthesis_meta_final$subclust_ident_lock)
# Astrocytes (dense bulbous)  Astrocytes (dense smooth)                        oRG                   Dividing                 Early OPCs 
#                        536                        786                       3410                        218                       1179 
#   ENs                    EX_IPCs                     IN_CGE                    IN_IPCs                   IN_local 
# 56582                      10956                       2373                       4595                       5402 
# IN_MGE                      IN_OB           Oligodendrocytes                       OPCs                         RG 
#   2624                        545                        772                       2431                       4529 
# tRG 
# 602
sum(is.na(ls_synthesis_meta_final$subclust_ident_lock))
# [1] 0
table(ls_synthesis_meta_final$subcluster_identity[ls_synthesis_meta_final$subclust_ident_lock == "Dividing"])
# RG 
# 218
ls_synthesis_meta_final$subclust_ident_lock_upset = ls_synthesis_meta_final$subclust_ident_lock
ls_synthesis_meta_final$subclust_ident_lock_upset[ls_synthesis_meta_final$subclust_ident_lock_upset %in% c("Astrocytes (dense bulbous)", "Astrocytes (dense smooth)")] = "Astrocytes"
ls_synthesis_meta_final$subclust_ident_lock_upset[ls_synthesis_meta_final$subclust_ident_lock_upset %in% c("oRG", "tRG", "RG", "Dividing")] = "RG"
ls_synthesis_meta_final$subclust_ident_lock_upset[ls_synthesis_meta_final$subclust_ident_lock_upset %in% c("Early OPCs")] = "OPCs"
table(ls_synthesis_meta_final$subclust_ident_lock_upset)
# Astrocytes              ENs          EX_IPCs           IN_CGE          IN_IPCs         IN_local           IN_MGE            IN_OB 
#       1322            56582            10956             2373             4595             5402             2624              545 
# Oligodendrocytes             OPCs               RG 
#              772             3610             8759
## get the new subclust idents up onto the actual object
## make sure all cells are in the final metadata object
ls_synthesis_v3_harmony_v3 <- subset(ls_synthesis_v3_harmony_v3, cells = ls_synthesis_meta_final$cell)
ls_synthesis_v3_harmony_v3 = AddMetaData(ls_synthesis_v3_harmony_v3, ls_synthesis_meta_final$subclust_ident_lock_upset, col.name = "subclust_ident_lock_upset")
DimPlot(ls_synthesis_v3_harmony_v3, group.by="subclust_ident_lock_upset")
ls_synthesis_v3_harmony_v3@meta.data$subclust_ident_lock_upset = factor(ls_synthesis_v3_harmony_v3@meta.data$subclust_ident_lock_upset, levels = c(
  "RG",
  "Astrocytes",
  "OPCs",
  "Oligodendrocytes",
  "EX_IPCs",
  "ENs",
  "IN_IPCs",
  "IN_local",
  "IN_CGE",
  "IN_OB",
  "IN_MGE"
))