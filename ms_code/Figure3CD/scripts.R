library(Seurat)
options(Seurat.object.assay.version = "v3")
beam1 <- Read10X("beam1/umi_count/", gene.column=1)
beam1_h <- CreateSeuratObject(counts = beam1)

beam2 <- Read10X("beam2/umi_count/", gene.column=1)
beam2_h <- CreateSeuratObject(counts = beam2)

x <- c('unmapped')
beam1_h_new <- subset(beam1_h,features=setdiff(rownames(beam1_h),x))
beam2_h_new <- subset(beam2_h,features=setdiff(rownames(beam2_h),x))
beam1_h_new[["orig.ident"]] <- "BEAM_pos"
beam2_h_new[["orig.ident"]] <- "BEAM_neg"

# Merge datasets into one single seurat object
alldata_h <- merge(beam1_h_new, beam2_h_new, add.cell.ids = c("BEAM_pos",  "BEAM_neg"))
Idents(alldata_h) <- "orig.ident"
table(Idents(alldata_h))

VlnPlot(alldata_h, features="nCount_RNA")

beam1 <- readRDS("beam1.rds")
beam2 <- readRDS("beam2.rds")

#at least 200 detected genes and genes need to be expressed in at least 3 cells
beam1_c <- WhichCells(beam1, expression = nFeature_RNA > 200)
beam1_f <- rownames(beam1)[Matrix::rowSums(beam1) > 3]

beam1.filt <- subset(beam1, features = beam1_f, cells = beam1_c)
dim(beam1.filt)

#at least 200 detected genes and genes need to be expressed in at least 3 cells
beam2_c <- WhichCells(beam2, expression = nFeature_RNA > 200)
beam2_f <- rownames(beam2)[Matrix::rowSums(beam2) > 3]

beam2.filt <- subset(beam2, features = beam2_f, cells = beam2_c)
dim(beam2.filt)

# perform visualization and clustering steps
beam2.filt <- NormalizeData(beam2.filt)
beam2.filt <- FindVariableFeatures(beam2.filt)
beam2.filt <- ScaleData(beam2.filt)
beam2.filt <- RunPCA(beam2.filt, verbose = FALSE)
beam2.filt <- FindNeighbors(beam2.filt, dims = 1:30)
beam2.filt <- FindClusters(beam2.filt, resolution = 0.8, verbose = FALSE)
beam2.filt <- RunUMAP(beam2.filt, dims = 1:30)
options(repr.plot.width=10, repr.plot.height=10)
DimPlot(beam2.filt, label = TRUE) + NoLegend()

beam1.filt[["orig.ident"]] <- "BEAM_pos"
beam2.filt[["orig.ident"]] <- "BEAM_neg"


# Merge datasets into one single seurat object
alldata <- merge(beam1.filt, beam2.filt, add.cell.ids = c("BEAM_pos", "BEAM_neg"))

# Way1: Doing it using Seurat function
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")
options(repr.plot.width=18, repr.plot.height=10)
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) +
    NoLegend()

selected_mito <- WhichCells(alldata, expression = percent_mito < 5)
selected_ribo <- WhichCells(alldata, expression = percent_ribo > 10)

# and subset the object to only keep those cells
data.filt <- subset(alldata, cells = selected_mito)
data.filt <- subset(alldata, cells = selected_ribo)

dim(data.filt)

table(data.filt$orig.ident)

table(alldata$orig.ident)

options(repr.plot.width=18, repr.plot.height=10)
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) +
    NoLegend()

data.filt <- NormalizeData(data.filt)
data.filt <- FindVariableFeatures(data.filt)
data.filt <- ScaleData(data.filt)
data.filt <- RunPCA(data.filt, verbose = FALSE)
data.filt <- FindNeighbors(data.filt, dims = 1:30)
data.filt <- FindClusters(data.filt, resolution = 0.8, verbose = FALSE)
data.filt <- RunUMAP(data.filt, dims = 1:30)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(data.filt, label = TRUE) + NoLegend()

options(repr.plot.width=8, repr.plot.height=8)
ElbowPlot(data.filt, ndims=50)

options(repr.plot.width=10, repr.plot.height=8)
DimPlot(data.filt, group.by = "orig.ident", label = FALSE)

library(SingleCellExperiment)
library(Seurat)
library(SingleR)
monaco_ann1 <-function(seur_obj) {
  monaco.ref <- celldex::MonacoImmuneData()
  seur_obj<- as.SingleCellExperiment(seur_obj)

  monaco.main <- SingleR(test = seur_obj,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)

  seur_obj$monaco.main <- monaco.main$pruned.labels
}

monaco_ann2 <-function(seur_obj) {
  
  monaco.ref <- celldex::MonacoImmuneData()
  seur_obj<- as.SingleCellExperiment(seur_obj)
  
  monaco.fine <- SingleR(test = seur_obj,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
  
  singleR_labels <- monaco.fine$pruned.labels
  t <- table(singleR_labels)
  other <- names(t)[t < 5]
  singleR_labels[singleR_labels %in% other] <- NA
  
  seur_obj$monaco.fine <- singleR_labels
  
}

data.filt$mono2<-monaco_ann2(data.filt)

options(repr.plot.width=12, repr.plot.height=8)
DimPlot(data.filt, reduction = "umap", group.by = "mono2", label = FALSE, label.size = 5 ,repel = TRUE) #+ theme(legend.text=element_text(size=20)) + NoAxes()

table(data.filt$mono2)

dim(data.filt)

data.filt <- subset(data.filt, subset=mono2 %in% c("Exhausted B cells", "Naive B cells", "Non-switched memory B cells", "Switched memory B cells"))

library(dittoSeq)

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(data.filt, reduction = "umap", group.by = "mono2", label = FALSE, label.size = 5 ,repel = TRUE, cols=dittoColors()) #+ theme(legend.text=element_text(size=20)) + NoAxes()

png("Fig3C1.png",width=5.75,height=3.5,units="in",res=1200)
DimPlot(data.filt, reduction = "umap", group.by = "mono2", label = FALSE, label.size = 5 ,repel = TRUE, cols=dittoColors()) #+ theme(legend.text=element_text(size=20)) + NoAxes()
dev.off()

options(repr.plot.width=10, repr.plot.height=8)
dittoBarPlot(
    object = data.filt,
    var = "mono2",
    group.by = "orig.ident") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15))

data.filt$orig.ident_modified <- data.filt$orig.ident

data.filt$orig.ident_modified <- gsub("BEAM_pos", "BEAM pos", data.filt$orig.ident_modified)
data.filt$orig.ident_modified <- gsub("BEAM_neg", "BEAM neg", data.filt$orig.ident_modified)

png("Fig3C2.png",width=4.75,height=3.5,units="in",res=1200)
dittoBarPlot(
    object = data.filt,
    var = "mono2",
    group.by = "orig.ident_modified") + theme(legend.text=element_text(size=9), axis.text=element_text(size=9)) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

#overlapping is so low, need to input GEX barcode as whitelist and rerun Cite-Seq-Count?
joint.bcs <- intersect(colnames(data.filt), paste0(colnames(alldata_h), '-1'))
length(joint.bcs)

# Subset RNA by joint cell barcodes
data.filt.sub <- data.filt[, joint.bcs]
options(repr.plot.width=10, repr.plot.height=8)
dittoBarPlot(
    object = data.filt.sub,
    var = "mono2",
    group.by = "orig.ident") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15))

options(repr.plot.width=10, repr.plot.height=8)
joint.bcs2 <- gsub('-1', '', joint.bcs)
alldata_h.sub <- alldata_h[, joint.bcs2]
VlnPlot(alldata_h.sub, features="nCount_RNA")

options(repr.plot.width=8, repr.plot.height=7)
DimPlot(data.filt.sub, group.by = "orig.ident", label = FALSE) + theme(legend.text=element_text(size=20)) + NoAxes()

# Merge datasets into one single seurat object
beam.data <- GetAssayData(object = alldata_h.sub[['RNA']], slot='data')
colnames(beam.data) <- paste0(colnames(beam.data), '-1')
data.filt.sub <- AddMetaData(data.filt.sub, beam.data@x, col.name = 'BEAM')

options(repr.plot.width=10, repr.plot.height=8)
FeaturePlot(data.filt.sub, features = 'BEAM')

joint.bcs1 <- intersect(colnames(data.filt), paste0('BEAM_pos_', colnames(beam1_h_new), '-1'))
length(joint.bcs1)

options(repr.plot.width=10, repr.plot.height=8)
joint.bcs12 <- gsub('-1', '', joint.bcs1)
joint.bcs123 <- gsub('BEAM_pos_', '', joint.bcs12)
beam1_h_new.sub <- beam1_h_new[, joint.bcs123]
VlnPlot(beam1_h_new.sub, features="nCount_RNA")

# Subset RNA by joint cell barcodes
data.filt.sub12 <- data.filt[, joint.bcs1]
options(repr.plot.width=10, repr.plot.height=8)
dittoBarPlot(
    object = data.filt.sub12,
    var = "mono2",
    group.by = "orig.ident") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15))

str(data.filt.sub12)

colnames(beam1_h_new.sub) <- paste0("BEAM_pos_", colnames(beam1_h_new.sub), "-1")
beam1_h_new.sub <- AddMetaData(beam1_h_new.sub, data.filt.sub12$mono2, col.name = 'mono2')

options(repr.plot.width=10, repr.plot.height=8)
VlnPlot(beam1_h_new.sub, features="nCount_RNA", group.by="mono2") + stat_summary(fun.y = median, geom='point', size = 25, colour = "white", shape = 95)

joint.bcs2 <- intersect(colnames(data.filt), paste0('BEAM_neg_', colnames(beam2_h_new), '-1'))
length(joint.bcs2)

options(repr.plot.width=10, repr.plot.height=8)
joint.bcs22 <- gsub('-1', '', joint.bcs2)
joint.bcs23 <- gsub('BEAM_neg_', '', joint.bcs22)
beam2_h_new.sub <- beam2_h_new[, joint.bcs23]
VlnPlot(beam2_h_new.sub, features="nCount_RNA")

# Subset RNA by joint cell barcodes
data.filt.sub22 <- data.filt[, joint.bcs2]
options(repr.plot.width=10, repr.plot.height=8)
dittoBarPlot(
    object = data.filt.sub22,
    var = "mono2",
    group.by = "orig.ident") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15))

colnames(beam2_h_new.sub) <- paste0("BEAM_neg_", colnames(beam2_h_new.sub), "-1")
beam2_h_new.sub <- AddMetaData(beam2_h_new.sub, data.filt.sub22$mono2, col.name = 'mono2')

options(repr.plot.width=10, repr.plot.height=8)
VlnPlot(beam2_h_new.sub, features="nCount_RNA", group.by="mono2") + stat_summary(fun.y = median, geom='point', size = 25, colour = "white", shape = 95)

seurat_obj <-readRDS("cite_beam_seurat.RDS")

seurat_obj_Pos <- subset(seurat_obj, Lane == 1)
seurat_obj_Neg <- subset(seurat_obj, Lane == 2)

colnames(seurat_obj_Pos) <- paste0("BEAM_pos_", colnames(seurat_obj_Pos))
colnames(seurat_obj_Neg) <- paste0("BEAM_neg_", colnames(seurat_obj_Neg))

joint.bcs1a <- intersect(colnames(seurat_obj_Pos), colnames(beam1_h_new.sub))
length(joint.bcs1a)

beam1_h_new.suba <- beam1_h_new.sub[, joint.bcs1a]

colnames(seurat_obj_Neg) <- gsub('-2', '-1', colnames(seurat_obj_Neg))

joint.bcs2a <- intersect(colnames(seurat_obj_Neg), colnames(beam2_h_new.sub))
length(joint.bcs2a)

beam2_h_new.suba <- beam2_h_new.sub[, joint.bcs2a]

seurat_obj_Pos_sub <- seurat_obj_Pos[, joint.bcs1a]
seurat_obj_Neg_sub <- seurat_obj_Neg[, joint.bcs2a]

ADT_Pos <- GetAssay(seurat_obj_Pos_sub, assay="CITE")

ADT_Pos$counts['Ig-light-chain-kappa',] - ADT_Pos$counts['IgD',]

beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['Ig-light-chain-kappa',] - ADT_Pos$counts['IgD',], col.name = 'IgKD')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['Ig-light-chain-lambda',] - ADT_Pos$counts['IgD',], col.name = 'IgLD')

beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['Ig-light-chain-lambda',] + ADT_Pos$counts['Ig-light-chain-kappa',] - 2*ADT_Pos$counts['IgD',], col.name = 'IgLKD')

beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['Ig-light-chain-kappa',] + ADT_Pos$counts['Ig-light-chain-lambda',], col.name = 'IgKL')

beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['Ig-light-chain-kappa',], col.name = 'IgK')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['Ig-light-chain-lambda',], col.name = 'IgL')

str(beam1_h_new.suba)

options(repr.plot.width=8, repr.plot.height=8)
FeatureScatter(object = beam1_h_new.suba, feature1 = 'nCount_RNA', feature2 = 'IgKD')

options(repr.plot.width=8, repr.plot.height=8)
FeatureScatter(object = beam1_h_new.suba, feature1 = 'nCount_RNA', feature2 = 'IgLD')

options(repr.plot.width=18, repr.plot.height=6)
FeatureScatter(object = beam1_h_new.suba, feature1 = 'nCount_RNA', feature2 = 'IgKD', split.by='mono2')

options(repr.plot.width=18, repr.plot.height=6)
FeatureScatter(object = beam1_h_new.suba, feature1 = 'nCount_RNA', feature2 = 'IgK', split.by='mono2')

options(repr.plot.width=18, repr.plot.height=6)
FeatureScatter(object = beam1_h_new.suba, feature1 = 'nCount_RNA', feature2 = 'IgL', split.by='mono2')

options(repr.plot.width=18, repr.plot.height=6)
FeatureScatter(object = beam1_h_new.suba, feature1 = 'nCount_RNA', feature2 = 'IgKL', split.by='mono2',plot.cor = TRUE)

test <- beam1_h_new.suba
test <- SetIdent(test, value="mono2")
options(repr.plot.width=18, repr.plot.height=12)
g <- FeatureScatter(object = test, feature1 = 'nCount_RNA', feature2 = 'IgKL', plot.cor = TRUE)
g + facet_wrap(~colors)

library(ggpubr)

nCount_Antigen = FetchData(test,"nCount_RNA")
colnames(nCount_Antigen) = "nCount_Antigen"
testf <- data.frame(cluster = Idents(test), nCount_Antigen, IgKL = FetchData(test,"IgKL"))
testf <- subset(testf, cluster!="NA")
ggplot(testf, aes(x = nCount_Antigen , y = IgKL, col = cluster)) +
geom_point(size=1) + 
facet_wrap(~cluster)+
NoLegend() +
stat_cor(method = "pearson", size = 10) +
theme(text = element_text(size = 20))

testf <- data.frame(cluster = Idents(test), nCount_Antigen, IgKLD = FetchData(test,"IgLKD"))
testf <- subset(testf, cluster!="NA")
ggplot(testf, aes(x = nCount_Antigen , y = IgLKD, col = cluster)) +
geom_point(size=1) + 
facet_wrap(~cluster)+
NoLegend() +
stat_cor(method = "pearson", size = 10) +
theme(text = element_text(size = 20))

testf <- data.frame(cluster = Idents(test), nCount_Antigen, IgK = FetchData(test,"IgK"))
testf <- subset(testf, cluster!="NA")
ggplot(testf, aes(x = nCount_Antigen , y = IgK, col = cluster)) +
geom_point(size=1) + 
facet_wrap(~cluster)+
NoLegend() +
stat_cor(method = "pearson", size = 10) +
theme(text = element_text(size = 20))

testf <- data.frame(cluster = Idents(test), nCount_Antigen, IgKD = FetchData(test,"IgKD"))
testf <- subset(testf, cluster!="NA")
ggplot(testf, aes(x = nCount_Antigen , y = IgKD, col = cluster)) +
geom_point(size=1) + 
facet_wrap(~cluster)+
NoLegend() +
stat_cor(method = "pearson", size = 10) +
theme(text = element_text(size = 20))

testf <- data.frame(cluster = Idents(test), nCount_Antigen, IgL = FetchData(test,"IgL"))
testf <- subset(testf, cluster!="NA")
ggplot(testf, aes(x = nCount_Antigen , y = IgL, col = cluster)) +
geom_point(size=1) + 
facet_wrap(~cluster)+
NoLegend() +
stat_cor(method = "pearson", size = 10) +
theme(text = element_text(size = 20))

testf <- data.frame(cluster = Idents(test), nCount_Antigen, IgLD = FetchData(test,"IgLD"))
testf <- subset(testf, cluster!="NA")
ggplot(testf, aes(x = nCount_Antigen , y = IgLD, col = cluster)) +
geom_point(size=1) + 
facet_wrap(~cluster)+
NoLegend() +
stat_cor(method = "pearson", size = 10) +
theme(text = element_text(size = 20))

options(repr.plot.width=18, repr.plot.height=6)
FeatureScatter(object = beam1_h_new.suba, feature1 = 'nCount_RNA', feature2 = 'IgLD', split.by='mono2')

options(repr.plot.width=18, repr.plot.height=6)
FeatureScatter(object = beam1_h_new.suba, feature1 = 'nCount_RNA', feature2 = 'IgLKD', split.by='mono2')

#'CD11c''CD19.1''CD21''CD27.1''CD307d''CD307e''CD32''CD72.1''CD85j''IgD'
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['Ig-light-chain-kappa',], col.name = 'IgK')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['Ig-light-chain-lambda',], col.name = 'IgL')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['IgD',], col.name = 'IgD')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['CD85j',], col.name = 'CD85j')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['CD72.1',], col.name = 'CD72')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['CD32',], col.name = 'CD32')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['CD307e',], col.name = 'CD307e')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['CD307d',], col.name = 'CD307d')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['CD27.1',], col.name = 'CD27')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['CD21',], col.name = 'CD21')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['CD19.1',], col.name = 'CD19')
beam1_h_new.suba <- AddMetaData(beam1_h_new.suba, ADT_Pos$counts['CD11c',], col.name = 'CD11c')

seurat_obj_Pos_sub <- AddMetaData(seurat_obj_Pos_sub, beam1_h_new.suba$mono2, col.name = 'mono2')

seurat_obj_Neg_sub <- AddMetaData(seurat_obj_Neg_sub, beam2_h_new.suba$mono2, col.name = 'mono2')

seurat_obj_Pos_sub$mono2[is.na(seurat_obj_Pos_sub$mono2)] <- 'None'
seurat_obj_Pos_sub$mono2 <- gsub('Plasmacytoid dendritic cells', 'None', seurat_obj_Pos_sub$mono2)

unique(seurat_obj_Pos_sub$mono2)

seurat_obj_Neg_sub$mono2[is.na(seurat_obj_Neg_sub$mono2)] <- 'None'
seurat_obj_Neg_sub$mono2 <- gsub('Plasmacytoid dendritic cells', 'None', seurat_obj_Neg_sub$mono2)
unique(seurat_obj_Neg_sub$mono2)

seurat_obj_Pos_sub[['CITE']]@var.features

seurat_obj_Pos_sub <- ScaleData(seurat_obj_Pos_sub)

options(repr.plot.width=16, repr.plot.height=12)
DoHeatmap(seurat_obj_Pos_sub, features = seurat_obj_Pos_sub[['CITE']]@var.features, assay = "CITE", group.by = 'mono2', slot='counts', angle = 45, group.colors =dittoColors()) + 
    theme(text = element_text(size = 20)) 

png("Fig3D1.png",width=5.25,height=3.5,units="in",res=1200)
DoHeatmap(seurat_obj_Pos_sub, features = seurat_obj_Pos_sub[['CITE']]@var.features, assay = "CITE", group.by = 'mono2', slot='counts', angle = 25, size = 2, group.colors =dittoColors()) + 
    theme(text = element_text(size = 8)) + 
    theme(legend.key.size = unit(0.5, "cm"),  # Controls the size of the color bar
        legend.title = element_text(size = 8),  # Controls the size of the legend title
        legend.text = element_text(size = 6))  # Controls the size of the legend text
dev.off()

f3c2 <- merge(seurat_obj_Pos_sub, seurat_obj_Neg_sub)
f3c2$beam <- ifelse(f3c2$Lane==1, "BEAM pos", "BEAM neg")
png("Fig3C2_new.png",width=4.75,height=3.5,units="in",res=1200)
dittoBarPlot(
    object = f3c2,
    var = "mono2",
    group.by = "beam") + theme(legend.text=element_text(size=9), axis.text=element_text(size=9)) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

options(repr.plot.width=16, repr.plot.height=12)
DoHeatmap(seurat_obj_Pos_sub, features = seurat_obj_Pos_sub[['CITE']]@var.features, assay = "CITE", group.by = 'mono2', slot='data', angle = 45) + 
    theme(text = element_text(size = 20))

options(repr.plot.width=16, repr.plot.height=12)
DoHeatmap(seurat_obj_Pos_sub, features = seurat_obj_Pos_sub[['CITE']]@var.features, assay = "CITE", group.by = 'mono2', angle = 45) + 
    theme(text = element_text(size = 20))

DefaultAssay(seurat_obj_Pos_sub) <- 'CITE'

atypical <- WhichCells(seurat_obj_Pos_sub, slot='data', expression = CD19.1 > 2 & CD27.1 < 2 & IgD < 2)
seurat_obj_Pos_sub$atypical <- ifelse(colnames(seurat_obj_Pos_sub) %in% atypical, "atypical", "not")

DefaultAssay(seurat_obj_Pos_sub) <- 'RNA'
seurat_obj_Pos_sub <- NormalizeData(seurat_obj_Pos_sub)
seurat_obj_Pos_sub <- FindVariableFeatures(seurat_obj_Pos_sub)
seurat_obj_Pos_sub <- ScaleData(seurat_obj_Pos_sub)
seurat_obj_Pos_sub <- RunPCA(seurat_obj_Pos_sub, verbose = FALSE)
seurat_obj_Pos_sub <- FindNeighbors(seurat_obj_Pos_sub, dims = 1:30)
seurat_obj_Pos_sub <- FindClusters(seurat_obj_Pos_sub, resolution = 0.8, verbose = FALSE)
seurat_obj_Pos_sub <- RunUMAP(seurat_obj_Pos_sub, dims = 1:30)

options(repr.plot.width=12, repr.plot.height=12)
DimPlot(seurat_obj_Pos_sub, group.by = "atypical", label = FALSE) + theme(legend.text=element_text(size=20)) + NoAxes()

options(repr.plot.width=16, repr.plot.height=12)
DimPlot(seurat_obj_Pos_sub, group.by = "mono2", label = FALSE) + theme(legend.text=element_text(size=20)) + NoAxes()

options(repr.plot.width=10, repr.plot.height=8)
dittoBarPlot(
    object = seurat_obj_Pos_sub,
    var = "mono2",
    group.by = "atypical") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15))

png("Fig3D2.png",width=5.75,height=3.5,units="in",res=1200)
