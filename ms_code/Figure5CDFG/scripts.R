library(Seurat)
options(Seurat.object.assay.version = "v3")
 
B1_US_SeuratObj1 <- readRDS('DL1.rds')
B1_US_SeuratObj2 <- readRDS('DL2.rds')
B1_US_SeuratObj3 <- readRDS('DL3.rds')

B1_US_merge = merge(B1_US_SeuratObj1, y = c(B1_US_SeuratObj2, B1_US_SeuratObj3), add.cell.ids=c("L1", "L2", "L3"))
mito.genes = grep(pattern = "^MT-", x = rownames(B1_US_merge), value = TRUE)
B1_US_merge<-AddMetaData(object = B1_US_merge, metadata = Matrix::colSums(B1_US_merge[mito.genes,])/Matrix::colSums(B1_US_merge), col.name = "percent.mito")

percent.ribo <- PercentageFeatureSet(B1_US_merge, pattern = "^RP[SL]")
B1_US_merge <- AddMetaData(B1_US_merge, percent.ribo, col.name = "percent.ribo")

options(repr.plot.width=18, repr.plot.height=6)
VlnPlot(B1_US_merge, c("nCount_ADT", "nCount_Ag", "percent.mito","nFeature_RNA", "nCount_RNA"), group.by = "Lane", pt.size=0.1, ncol=5)

table(B1_US_merge$Best)

CV <- subset(B1_US_merge, subset = Best %in% c("CHI.017", "D001", "D012", "D024"))
CVf <- subset(CV, subset = nCount_ADT < 6000 & nCount_Ag < 3000 & nCount_ADT > 200)
CVf2 = subset(CVf, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA < 6000 & nCount_RNA > 200 & percent.mito < 0.08)

c(ncol(CVf2), ncol(CVf), ncol(CV), ncol(B1_US_merge))

options(repr.plot.width=18, repr.plot.height=4)
VlnPlot(CVf2, c("nCount_ADT", "nCount_Ag", "percent.mito","nFeature_RNA", "nCount_RNA"), group.by = "Lane", pt.size=0, ncol=5)

# perform visualization and clustering steps
SNG.US <- NormalizeData(CVf2)
SNG.US <- FindVariableFeatures(SNG.US)
SNG.US <- ScaleData(SNG.US)
SNG.US <- RunPCA(SNG.US, verbose = FALSE)
SNG.US <- FindNeighbors(SNG.US, dims = 1:30)
SNG.US <- FindClusters(SNG.US, resolution = 0.8, verbose = FALSE)
SNG.US <- RunUMAP(SNG.US, dims = 1:30)
options(repr.plot.width=10, repr.plot.height=10)
DimPlot(SNG.US, label = TRUE) + NoLegend()

library(Azimuth)
library(dittoSeq)
library(DoubletFinder)
library(SingleCellExperiment)
library(SeuratDisk)
library(patchwork)
library(dplyr)

DefaultAssay(SNG.US) <- 'RNA'
SNG.US <- RunAzimuth(SNG.US, reference = "pbmcref")

options(repr.plot.width=12, repr.plot.height=8)
DimPlot(SNG.US, reduction = "umap", group.by = 'predicted.celltype.l2', label = FALSE, label.size = 5 ,repel = TRUE) #+ theme(legend.text=element_text(size=20)) + NoAxes()

library(SingleR)

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

SNG.US$mono2<-monaco_ann2(SNG.US)

options(repr.plot.width=12, repr.plot.height=8)
DimPlot(SNG.US, reduction = "umap", group.by = 'mono2', label = FALSE, label.size = 5 ,repel = TRUE) #+ theme(legend.text=element_text(size=20)) + NoAxes()

SNG.US <- subset(SNG.US, subset=mono2 %in% c("Exhausted B cells", "Naive B cells", "Non-switched memory B cells", "Switched memory B cells"))

options(repr.plot.width=12, repr.plot.height=8)
DimPlot(SNG.US, reduction = "umap", group.by = "mono2", label = FALSE, label.size = 5 ,repel = TRUE, cols=dittoColors()) #+ theme(legend.text=element_text(size=20)) + NoAxes()

options(repr.plot.width=5, repr.plot.height=6)
dittoBarPlot(
    object = SNG.US,
    var = "mono2",
    group.by = "orig.ident") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15))

options(repr.plot.width=12, repr.plot.height=10)
FeaturePlot(SNG.US, features = c('IgD', 'CD27.1', 'HSA', 'D614G'), max.cutoff = 200)

options(repr.plot.width=12, repr.plot.height=10)
FeaturePlot(SNG.US, reduction = 'umap', features = c('DV1', 'DV2', 'DV3', 'DV4'), max.cutoff = 200)

AgData <- GetAssayData(object = SNG.US, assay = "Ag", layer = "counts")
Ag_Score <- (1 - pbeta(0.925, AgData['D614G',] + 1, AgData['HSA',] + 3)) * 100
SNG.US <- AddMetaData(SNG.US, metadata = Ag_Score, col.name = "Ag_Score")

D614G_cells <- names(which(SNG.US$Ag_Score > 90))
SNG.US$D614G <- ifelse(colnames(SNG.US) %in% D614G_cells, "Pos", "Neg")

SNG.US_D614G <- subset(SNG.US, subset = D614G == 'Pos')
options(repr.plot.width=5, repr.plot.height=6)
dittoBarPlot(
    object = SNG.US_D614G,
    var = "mono2",
    group.by = "orig.ident") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle('D614G')

options(repr.plot.width=5, repr.plot.height=6)
dittoBarPlot(
    object = SNG.US,
    var = "mono2",
    group.by = "D614G") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle('D614G')

table(SNG.US$D614G)

options(repr.plot.width=9, repr.plot.height=4)
VlnPlot(SNG.US, c("nCount_ADT", "nCount_Ag"), group.by = "mono2", pt.size=0, ncol=2)

options(repr.plot.width=6, repr.plot.height=4)
VlnPlot(SNG.US, c("nCount_ADT", "nCount_Ag"), group.by = "D614G", pt.size=0, ncol=2)

options(repr.plot.width=12, repr.plot.height=8)
DimPlot(SNG.US, reduction = "umap", group.by = "D614G", label = FALSE, label.size = 5 ,repel = TRUE, cols=dittoColors()) #+ theme(legend.text=element_text(size=20)) + NoAxes()


rownames(SNG.US[['ADT']])

options(repr.plot.width=16, repr.plot.height=12)
DoHeatmap(SNG.US, features = rownames(SNG.US[['ADT']]), assay = "ADT", group.by = 'mono2', slot='counts', angle = 45, group.colors =dittoColors()) + 
    theme(text = element_text(size = 20)) 

#data always contains the log-normed version of counts
options(repr.plot.width=16, repr.plot.height=12)
DoHeatmap(SNG.US, features = rownames(SNG.US[['ADT']]), assay = "ADT", group.by = 'mono2', slot='data', angle = 45) + 
    theme(text = element_text(size = 20))

# For DVs
ag_data <- GetAssayData(object = SNG.US, assay = "Ag", layer = "counts")

DV4 <- subset(SNG.US, cells = colnames(ag_data)[ag_data['DV4', ] > 25 & ag_data['HSA', ] < 25])
DV3 <- subset(SNG.US, cells = colnames(ag_data)[ag_data['DV3', ] > 25 & ag_data['HSA', ] < 25])
DV2 <- subset(SNG.US, cells = colnames(ag_data)[ag_data['DV2', ] > 25 & ag_data['HSA', ] < 25])
DV1 <- subset(SNG.US, cells = colnames(ag_data)[ag_data['DV1', ] > 25 & ag_data['HSA', ] < 25])
CVD <- subset(SNG.US, cells = colnames(ag_data)[ag_data['D614G', ] > 25 & ag_data['HSA', ] < 25])

# Subset RNA by joint cell barcodes
options(repr.plot.width=4, repr.plot.height=8)
p1 <- dittoBarPlot(
    object = DV1,
    var = "mono2",
    group.by = "orig.ident") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("DV1")
p2 <- dittoBarPlot(
    object = DV2,
    var = "mono2",
    group.by = "orig.ident") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("DV2")
p3 <- dittoBarPlot(
    object = DV3,
    var = "mono2",
    group.by = "orig.ident") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("DV3")
p4 <- dittoBarPlot(
    object = DV4,
    var = "mono2",
    group.by = "orig.ident") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("DV4")
p1+p2+p3+p4

options(repr.plot.width=5.4, repr.plot.height=4)
dittoBarPlot(
    object = CVD,
    var = "mono2",
    group.by = "orig.ident") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("D614G")

DV <- subset(SNG.US, cells = colnames(ag_data)[(ag_data['DV4', ] > 25 | ag_data['DV3', ] > 25 | ag_data['DV2', ] > 25 | ag_data['DV1', ] > 25) & ag_data['HSA', ] < 25 & ag_data['D614G', ] < 25])

options(repr.plot.width=5.4, repr.plot.height=4)
dittoBarPlot(
    object = DV,
    var = "mono2",
    group.by = "orig.ident") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("DV")

DV12 <- subset(DV, subset=orig.ident%in%c("L1", "L2"))

options(repr.plot.width=5.4, repr.plot.height=4)
dittoBarPlot(
    object = DV12,
    var = "mono2",
    group.by = "Best") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("DV")

table(DV$Best)
table(DV1$Best)
table(DV2$Best)
table(DV3$Best)
table(DV4$Best)

table(DV$Best, DV$Lane)
table(DV1$Best, DV1$Lane)
table(DV2$Best, DV2$Lane)
table(DV3$Best, DV3$Lane)
table(DV4$Best, DV4$Lane)

table(SNG.US$Best)

round(table(DV$Best)/table(SNG.US$Best)*1000)/10
round(table(DV1$Best)/table(SNG.US$Best)*1000)/10
round(table(DV2$Best)/table(SNG.US$Best)*1000)/10
round(table(DV3$Best)/table(SNG.US$Best)*1000)/10
round(table(DV4$Best)/table(SNG.US$Best)*1000)/10

nonDV <- subset(SNG.US, Lane==3)

options(repr.plot.width=4.8, repr.plot.height=4)
dittoBarPlot(
    object = nonDV,
    var = "mono2",
    group.by = "orig.ident") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("nonDV")

options(repr.plot.width=5.6, repr.plot.height=4)
dittoBarPlot(
    object = nonDV,
    var = "mono2",
    group.by = "Best") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("nonDV")

DV1 <- subset(SNG.US, cells = colnames(ag_data)[ag_data['DV1', ] > 25 & ag_data['HSA', ] < 25 & ag_data['D614G', ] < 25])
DV2 <- subset(SNG.US, cells = colnames(ag_data)[ag_data['DV2', ] > 25 & ag_data['HSA', ] < 25 & ag_data['D614G', ] < 25])
DV3 <- subset(SNG.US, cells = colnames(ag_data)[ag_data['DV3', ] > 25 & ag_data['HSA', ] < 25 & ag_data['D614G', ] < 25])
DV4 <- subset(SNG.US, cells = colnames(ag_data)[ag_data['DV4', ] > 25 & ag_data['HSA', ] < 25 & ag_data['D614G', ] < 25])

DV12 <- subset(DV1, subset=orig.ident%in%c("L1", "L2"))
DV22 <- subset(DV2, subset=orig.ident%in%c("L1", "L2"))
DV32 <- subset(DV3, subset=orig.ident%in%c("L1", "L2"))
DV42 <- subset(DV4, subset=orig.ident%in%c("L1", "L2"))

options(repr.plot.width=5.4, repr.plot.height=4)
dittoBarPlot(
    object = DV12,
    var = "mono2",
    group.by = "Best") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("DV1")

options(repr.plot.width=5.4, repr.plot.height=4)
dittoBarPlot(
    object = DV22,
    var = "mono2",
    group.by = "Best") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("DV2")

options(repr.plot.width=5.4, repr.plot.height=4)
dittoBarPlot(
    object = DV32,
    var = "mono2",
    group.by = "Best") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("DV3")

options(repr.plot.width=5.4, repr.plot.height=4)
dittoBarPlot(
    object = DV42,
    var = "mono2",
    group.by = "Best") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("DV4")

D614G <- subset(SNG.US, subset=D614G=="Pos")
options(repr.plot.width=5.4, repr.plot.height=4)
dittoBarPlot(
    object = D614G,
    var = "mono2",
    group.by = "Best") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("Ag+")

png("Fig5C1.png",width=5.75,height=3.5,units="in",res=1200)
dittoBarPlot(
    object = D614G,
    var = "mono2",
    group.by = "Best") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("Ag+")
dev.off()

nonD614G <- subset(SNG.US, subset=D614G=="Neg")
options(repr.plot.width=5.4, repr.plot.height=4)
dittoBarPlot(
    object = nonD614G,
    var = "mono2",
    group.by = "Best") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("Ag-")

png("Fig5C2.png",width=5.75,height=3.5,units="in",res=1200)
dittoBarPlot(
    object = nonD614G,
    var = "mono2",
    group.by = "Best") + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("Ag-")
dev.off()

options(repr.plot.width=16, repr.plot.height=12)
DoHeatmap(D614G, features = rownames(SNG.US[['ADT']]), assay = "ADT", group.by = 'mono2', slot='counts', angle = 45, group.colors =dittoColors()) + 
    theme(text = element_text(size = 20)) 

png("Fig5D.png",width=5.25,height=3.5,units="in",res=1200)
DoHeatmap(D614G, features = rownames(SNG.US[['ADT']]), assay = "ADT", group.by = 'mono2', slot='counts', angle = 25, size = 2, group.colors =dittoColors()) + 
    theme(text = element_text(size = 8)) + 
    theme(legend.key.size = unit(0.5, "cm"),  # Controls the size of the color bar
        legend.title = element_text(size = 8),  # Controls the size of the legend title
        legend.text = element_text(size = 6))  # Controls the size of the legend text
dev.off()

SNG.US0 <- subset(SNG.US, subset=Best=='D012')
Ag12 <- GetAssayData(object = SNG.US0, assay = "Ag", layer = "counts")
Agd12 <- data.frame(t(Ag12))

SNG.US0 <- subset(SNG.US, subset=Best=='CHI.017')
Ag17 <- GetAssayData(object = SNG.US0, assay = "Ag", layer = "counts")
Agd17 <- data.frame(t(Ag17))

options(repr.plot.width=18, repr.plot.height=9)
Ag <- c('DV1', 'DV2', 'DV3','DV4','BA1','HSA','D614G')
vv <- list()
for(i in 1:4){
    vv[[i]] <- ggplot(Agd12) + geom_point(aes_string(x=Ag[6], y=Ag[i]), size=0.5, shape=1) + scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 500))
    vv[[i+4]] <- ggplot(Agd17) + geom_point(aes_string(x=Ag[6], y=Ag[i]), size=0.5, shape=1) + scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 500))
}  
combined_plot <- wrap_plots(vv, ncol=4)
print(combined_plot)

png("Fig5F.png",width=5.25,height=2.5,units="in",res=1200)
Ag <- c('DV1', 'DV2', 'DV3','DV4','BA1','HSA','D614G')
vv <- list()
for(i in 1:4){
    vv[[i]] <- ggplot(Agd12) + geom_point(aes_string(x=Ag[6], y=Ag[i]), size=0.2, shape=1) + scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 500)) +
               theme(
                 axis.title.x = element_text(size = 5),  # X-axis label font size
                 axis.title.y = element_text(size = 5),  # Y-axis label font size
                 axis.text.x = element_text(size = 4),   # X-axis tick label font size
                 axis.text.y = element_text(size = 4)    # Y-axis tick label font size
               )
    vv[[i+4]] <- ggplot(Agd17) + geom_point(aes_string(x=Ag[6], y=Ag[i]), size=0.2, shape=1) + scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 500)) +
               theme(
                 axis.title.x = element_text(size = 5),  # X-axis label font size
                 axis.title.y = element_text(size = 5),  # Y-axis label font size
                 axis.text.x = element_text(size = 4),   # X-axis tick label font size
                 axis.text.y = element_text(size = 4)    # Y-axis tick label font size
               )
}  
combined_plot <- wrap_plots(vv, ncol=4)
print(combined_plot)
dev.off()

DV <- subset(SNG.US, cells = colnames(ag_data)[(ag_data['DV4', ] > 25 | ag_data['DV3', ] > 25 | ag_data['DV2', ] > 25 | ag_data['DV1', ] > 25) & ag_data['HSA', ] < 25 & ag_data['D614G', ] < 25])
SNG.US$DV <- ifelse(colnames(SNG.US) %in% colnames(DV), "DV", "nonDV")

CHI.017 <- subset(SNG.US, subset=Best=='CHI.017')
D001 <- subset(SNG.US, subset=Best=='D001')
D012 <- subset(SNG.US, subset=Best=='D012')
D024 <- subset(SNG.US, subset=Best=='D024')

options(repr.plot.width=8, repr.plot.height=4)
p1 <- dittoBarPlot(
    object = CHI.017,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("CHI.017")
p2 <- dittoBarPlot(
    object = D001,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("D001")
p3 <- dittoBarPlot(
    object = D012,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("D012")
p4 <- dittoBarPlot(
    object = D024,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("D024")


combined_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)

# Print the combined plot
print(combined_plot)

png("Fig5G.png",width=5.25,height=2.5,units="in",res=1200)
p1 <- dittoBarPlot(
    object = CHI.017,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=9), axis.text=element_text(size=8)) + ggtitle("CHI.017") 
p2 <- dittoBarPlot(
    object = D001,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=9), axis.text=element_text(size=8)) + ggtitle("D001") 
p3 <- dittoBarPlot(
    object = D012,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=9), axis.text=element_text(size=8)) + ggtitle("D012") 
p4 <- dittoBarPlot(
    object = D024,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=9), axis.text=element_text(size=8)) + ggtitle("D024") 
combined_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
print(combined_plot)
dev.off()

# Regenrate the plot with DV negative from Lane 3
nonDV$DV <- 'nonDV'
DV$DV <- 'DV'
DVnonDV <- merge(DV, nonDV) 

CHI.017 <- subset(DVnonDV, subset=Best=='CHI.017')
D001 <- subset(DVnonDV, subset=Best=='D001')
D012 <- subset(DVnonDV, subset=Best=='D012')
D024 <- subset(DVnonDV, subset=Best=='D024')

options(repr.plot.width=8, repr.plot.height=4)
p1 <- dittoBarPlot(
    object = CHI.017,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("CHI.017")
p2 <- dittoBarPlot(
    object = D001,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("D001")
p3 <- dittoBarPlot(
    object = D012,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("D012")
p4 <- dittoBarPlot(
    object = D024,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=15), axis.text=element_text(size=15)) + ggtitle("D024")
combined_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)

# Print the combined plot
print(combined_plot)

png("Fig5G.png",width=5.25,height=2.5,units="in",res=1200)
p1 <- dittoBarPlot(
    object = CHI.017,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=9), axis.text=element_text(size=8)) + ggtitle("CHI.017") 
p2 <- dittoBarPlot(
    object = D001,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=9), axis.text=element_text(size=8)) + ggtitle("D001") 
p3 <- dittoBarPlot(
    object = D012,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=9), axis.text=element_text(size=8)) + ggtitle("D012") 
p4 <- dittoBarPlot(
    object = D024,
    var = "mono2",
    group.by = "DV") + NoLegend() + theme(legend.text=element_text(size=9), axis.text=element_text(size=8)) + ggtitle("D024") 
combined_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
print(combined_plot)
dev.off()
