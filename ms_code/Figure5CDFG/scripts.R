library(Seurat)
options(Seurat.object.assay.version = "v3")
rawdir1 <- "/Volumes/chi/PROJECTS/23-303_Katzelnick/C-303-4_Pilot_7_Ag_panel_with_sequencing/Bioinformatics/Runs/240703_VH00286_150_AAF25TGM5/RUN/"

Lanes = c(3,4)
B1_US_data = list()
for(i in 1:length(Lanes)){
  B1_US_data[[i]] = Read10X_h5(paste0(rawdir1,"multi_config_",Lanes[i],"/outs/per_sample_outs/multi_config_",Lanes[i],"/count/sample_filtered_feature_bc_matrix.h5"), use.names=T)
}
rawdir2 <- "/Volumes/chi/PROJECTS/23-303_Katzelnick/C-303-4_Pilot_7_Ag_panel_with_sequencing/Bioinformatics/Runs/240607_VH00286_145_AACCFK3HV/RUN/"
B1_US_data[[3]] = Read10X_h5(paste0(rawdir2,"7multi_config_15/outs/per_sample_outs/7multi_config_15/count/sample_filtered_feature_bc_matrix.h5"), use.names=T)

B1_US_SeuratObj = list()
for(i in 1:3){
  B1_US_SeuratObj[[i]] <- CreateSeuratObject(counts = B1_US_data[[i]]$'Gene Expression', assay = "RNA", min.feature = 0, project=paste0("L",i))
  B1_US_SeuratObj[[i]][["ADT"]] <- CreateAssayObject(counts = B1_US_data[[i]]$'Antibody Capture'[1:12,])
  B1_US_SeuratObj[[i]][["Ag"]] <- CreateAssayObject(counts = B1_US_data[[i]]$'Antibody Capture'[13:19,])  
}

r1 <- read.table(gzfile("/Volumes/chi/PROJECTS/23-303_Katzelnick/C-303-4_Pilot_7_Ag_panel_with_sequencing/Bioinformatics/Runs/240703_VH00286_150_AAF25TGM5/RUN/SNP/DEMUX/multi_config_3/posterior_probabilities_refined.tsv.gz"),header=T)
rownames(r1) <- r1$BARCODE
r1$BARCODE <- NULL

# Function to extract and concatenate identifiers
extract_and_concatenate <- function(filename) {
  # Extract all matches of identifiers
  matches <- gregexpr("(CHI\\.|D0)\\d+", filename)
  identifiers <- regmatches(filename, matches)[[1]]
  # Concatenate identifiers with '+'
  paste(identifiers, collapse = "+")
}

# Apply the function to each filename
colnames(r1) <- sapply(colnames(r1), extract_and_concatenate)

r1[r1<0.95] <- 0
r1[r1>=0.95] <- 1  

r2 <- read.table(gzfile("/Volumes/chi/PROJECTS/23-303_Katzelnick/C-303-4_Pilot_7_Ag_panel_with_sequencing/Bioinformatics/Runs/240703_VH00286_150_AAF25TGM5/RUN/SNP/DEMUX/multi_config_4/posterior_probabilities_refined.tsv.gz"),header=T)
rownames(r2) <- r2$BARCODE
r2$BARCODE <- NULL

# Apply the function to each filename
colnames(r2) <- sapply(colnames(r2), extract_and_concatenate)

r2[r2<0.95] <- 0
r2[r2>=0.95] <- 1

r3 <- read.table(gzfile("/Volumes/chi/PROJECTS/23-303_Katzelnick/C-303-4_Pilot_7_Ag_panel_with_sequencing/Bioinformatics/Runs/240607_VH00286_145_AACCFK3HV/RUN/SNP/DEMUX/multi_config_15/posterior_probabilities_refined.tsv.gz"),header=T)
rownames(r3) <- r3$BARCODE
r3$BARCODE <- NULL

# Apply the function to each filename
colnames(r3) <- sapply(colnames(r3), extract_and_concatenate)

r3[r3<0.95] <- 0
r3[r3>=0.95] <- 1

#remove rows with all zeros!!!
r1 <- r1[apply(r1, 1, function(x) !all(x<1)),]
max_colnames <- apply(r1, 1, function(row) colnames(r1)[which.max(row)])
best <- as.data.frame(max_colnames)
colnames(best) <- c("Best")
B1_US_SeuratObj[[1]] <- AddMetaData(object = B1_US_SeuratObj[[1]], metadata = best)

#remove rows with all zeros!!!
r2 <- r2[apply(r2, 1, function(x) !all(x<1)),]
max_colnames <- apply(r2, 1, function(row) colnames(r2)[which.max(row)])
best <- as.data.frame(max_colnames)
colnames(best) <- c("Best")
B1_US_SeuratObj[[2]] <- AddMetaData(object = B1_US_SeuratObj[[2]], metadata = best)

#remove rows with all zeros!!!
r3 <- r3[apply(r3, 1, function(x) !all(x<1)),]
max_colnames <- apply(r3, 1, function(row) colnames(r3)[which.max(row)])
best <- as.data.frame(max_colnames)
colnames(best) <- c("Best")
B1_US_SeuratObj[[3]] <- AddMetaData(object = B1_US_SeuratObj[[3]], metadata = best)

# https://ucdavis-bioinformatics-training.github.io/2020-Advanced_Single_Cell_RNA_Seq/data_analysis/VDJ_Analysis_fixed
add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
    tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))

    # Remove the -1 at the end of each barcode.
    # Subsets so only the first line of each barcode is kept,
    # as each entry for given barcode will have same clonotype.
    tcr <- tcr[!duplicated(tcr$barcode), ]

    # Only keep the barcode and clonotype columns. 
    # We'll get additional clonotype info from the clonotype table.
    tcr <- tcr[,c("barcode", "raw_clonotype_id")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

    # Clonotype-centric info.
    clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep=""))

    # Slap the AA sequences onto our original table by clonotype_id.
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
    names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"

    # Reorder so barcodes are first column and set them as rownames.
    tcr <- tcr[, c(2,1,3)]
    rownames(tcr) <- tcr[,1]
    tcr[,1] <- NULL
    colnames(tcr) <- paste(type, colnames(tcr), sep="_")
    # Add to the Seurat object's metadata.
    clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
    return(clono_seurat)
}

for(i in 1:length(Lanes)){
    B1_US_SeuratObj[[i]] <- add_clonotype(paste0(rawdir1,"multi_config_",Lanes[i],"/outs/per_sample_outs/multi_config_",Lanes[i],"/vdj_b/"), B1_US_SeuratObj[[i]], "b")
}
B1_US_SeuratObj[[3]] <- add_clonotype(paste0(rawdir2,"7multi_config_",15,"/outs/per_sample_outs/7multi_config_",15,"/vdj_b/"), B1_US_SeuratObj[[3]], "b")

for(i in 1:3){
  B1_US_SeuratObj[[i]] <- RenameCells(B1_US_SeuratObj[[i]], new.names = paste0(substr(colnames(B1_US_SeuratObj[[i]]), start = 1, stop = 17),i))
  B1_US_SeuratObj[[i]]$Lane  <- rep(i, length(colnames(B1_US_SeuratObj[[i]])))
}  

B1_US_merge = merge(B1_US_SeuratObj[[1]], B1_US_SeuratObj[2:3], add.cell.ids=c("L1", "L2", "L3"))
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
source('/Users/wangl43/code/functions.R')
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

# Switch to BCR analysis
library(scRepertoire)

S1 <- read.csv(paste0(rawdir1,"multi_config_",Lanes[1],"/outs/per_sample_outs/multi_config_",Lanes[1],"/vdj_b/filtered_contig_annotations.csv"))
S2 <- read.csv(paste0(rawdir1,"multi_config_",Lanes[2],"/outs/per_sample_outs/multi_config_",Lanes[2],"/vdj_b/filtered_contig_annotations.csv"))
S3 <- read.csv(paste0(rawdir2,"7multi_config_15/outs/per_sample_outs/7multi_config_15/vdj_b/filtered_contig_annotations.csv"))

contig_list <- list(S1, S2, S3)
combined <- combineBCR(contig_list, samples=c("L1", "L2", "L3"), threshold=0.85)

table(SNG.US$Lane)

CHI <- subset(SNG.US, subset = Best %in% c("CHI.017", "D001", "D012", "D024"))
pos_barcodes <- rownames(CHI@meta.data[CHI@meta.data$Lane == 1 | CHI@meta.data$Lane == 2, ])
neg_barcodes <- rownames(CHI@meta.data[CHI@meta.data$Lane == 3, ])

length(pos_barcodes)
length(neg_barcodes)
head(neg_barcodes)

head(combined$L1$barcode)
head(combined$L2$barcode)
head(combined$L3$barcode)

pos_barcodes <- gsub('-2', '-1', pos_barcodes)
neg_barcodes <- gsub('-3', '-1', neg_barcodes)

combined$L1$sample <- ifelse(combined$L1$barcode %in% pos_barcodes, "Pos", combined$L1$sample)
combined$L1$sample <- ifelse(combined$L1$barcode %in% neg_barcodes, "Neg", combined$L1$sample)
combined$L2$sample <- ifelse(combined$L2$barcode %in% pos_barcodes, "Pos", combined$L2$sample)
combined$L2$sample <- ifelse(combined$L2$barcode %in% neg_barcodes, "Neg", combined$L2$sample)
combined$L3$sample <- ifelse(combined$L3$barcode %in% pos_barcodes, "Pos", combined$L3$sample)
combined$L3$sample <- ifelse(combined$L3$barcode %in% neg_barcodes, "Neg", combined$L3$sample)

table(combined$L1$sample)
table(combined$L2$sample)
table(combined$L3$sample)

# Combine all data frames into one
combined_data <- bind_rows(combined)

# Remove DV cells
combined_data <- subset(combined_data, subset=sample %in% c("Pos", "Neg"))

# Split into a list based on the sample column
reorganized_list <- combined_data %>%
  group_split(sample, .keep = TRUE) %>%
  setNames(unique(combined_data$sample))

options(repr.plot.width=4, repr.plot.height=6)
clonalHomeostasis(reorganized_list, cloneCall = "gene")

options(repr.plot.width=3, repr.plot.height=6)
clonalProportion(reorganized_list, cloneCall = "gene") 

options(repr.plot.width=3, repr.plot.height=6)
clonalOverlap(reorganized_list, 
              cloneCall = "gene+nt", 
              method = "morisita")

head(S1$barcode)
head(S2$barcode)
head(S3$barcode)

# Need to start with S1,S2 since BCRcombine dropped columns from 31 to 12, including the chain column used in plots below
pos_barcodes <- gsub('L1_', '', pos_barcodes)
pos_barcodes <- gsub('L2_', '', pos_barcodes)
pos_barcodes <- gsub('L3_', '', pos_barcodes)
neg_barcodes <- gsub('L1_', '', neg_barcodes)
neg_barcodes <- gsub('L2_', '', neg_barcodes)
neg_barcodes <- gsub('L3_', '', neg_barcodes)

S1 <- subset(S1, barcode %in% pos_barcodes | barcode %in% neg_barcodes)
S2 <- subset(S2, barcode %in% c(pos_barcodes, neg_barcodes))
S3 <- subset(S3, barcode %in% c(pos_barcodes, neg_barcodes))

#S1$sample <-ifelse(S1$barcode %in% pos_barcodes, "Pos", "Neg")
#S2$sample <-ifelse(S2$barcode %in% pos_barcodes, "Pos", "Neg")
#S3$sample <-ifelse(S3$barcode %in% pos_barcodes, "Pos", "Neg")
S1$sample <- 'Pos'
S2$sample <- 'Pos'
S3$sample <- 'Neg'

contig_list <- list(S1, S2, S3)

# Combine all data frames into one
combined_S <- bind_rows(contig_list)

# Split into a list based on the sample column
reorganized_S <- combined_S %>%
  group_split(sample, .keep = TRUE) %>%
  setNames(unique(combined_data$sample))

library(tm)
myTdm <- as.matrix(TermDocumentMatrix(reorganized_S$Pos$chain))
FreqMat <- data.frame(chain = rownames(myTdm), 
                      Freq = rowSums(myTdm), 
                      row.names = NULL)
myTdm <- as.matrix(TermDocumentMatrix(reorganized_S$Neg$chain))
FreqMat2 <- data.frame(chain = rownames(myTdm), 
                      Freq = rowSums(myTdm), 
                      row.names = NULL)
FM <- rbind(data.frame(id = "D614G_pos", FreqMat),
      data.frame(id = "D614G_neg", FreqMat2))

options(repr.plot.width=5, repr.plot.height=6)
FM%>%
  ggplot(aes(x = id, y = Freq, fill = chain)) +
  geom_bar(stat = "identity",position = "fill")+
  scale_y_continuous(labels = scales::percent) + theme(legend.text=element_text(size=15), axis.text=element_text(size=15))

myTdm <- as.matrix(TermDocumentMatrix(reorganized_S$Pos$c_gene))
FreqMat <- data.frame(chain = toupper(rownames(myTdm)), 
                      Freq = rowSums(myTdm), 
                      row.names = NULL)
myTdm <- as.matrix(TermDocumentMatrix(reorganized_S$Neg$c_gene))
FreqMat2 <- data.frame(chain = toupper(rownames(myTdm)), 
                      Freq = rowSums(myTdm), 
                      row.names = NULL)
FM <- rbind(data.frame(id = "D614G_pos", FreqMat),
      data.frame(id = "D614G_neg", FreqMat2))

options(repr.plot.width=5, repr.plot.height=6)
FM%>%
  ggplot(aes(x = id, y = Freq, fill = chain)) +
  geom_bar(stat = "identity",position = "fill")+
  scale_y_continuous(labels = scales::percent) + theme(legend.text=element_text(size=15), axis.text=element_text(size=15))

# Going back to reorganized_list
options(repr.plot.width=4, repr.plot.height=6)
clonalQuant(reorganized_list, cloneCall="gene+nt", chain = "IGH")

options(repr.plot.width=8, repr.plot.height=6)
clonalCompare(reorganized_list, 
                  top.clones = 10, 
                  samples = c("Pos", "Neg"), 
                  cloneCall="aa", 
                  graph = "alluvial")

#https://opig.stats.ox.ac.uk/webapps/covabdab/
clonalCompare(reorganized_list, 
                  top.clones = 10, 
                  samples = c("Pos", "Neg"), 
                  cloneCall="aa", 
                  exportTable = T)

library(Biostrings)
library(dplyr)

similarity <- function(target_seqs) {
    # Load and preprocess the CoV-AbDab data
    covid <- read.csv("/Users/wangl43/work/test/306_4/CoV-AbDab_080224.csv") %>%
        filter(grepl("uman", Heavy.V.Gene)) %>%
        mutate(CDRH3 = gsub("\\s+", "", CDRH3)) %>% # Remove whitespaces
        filter(nchar(CDRH3) > 6)

    # Function to compute the best match for a single target sequence
    compute_best_match <- function(target_seq) {
        # Align all CDRH3 sequences to the target sequence
        alignments <- lapply(covid$CDRH3, function(cdrh3) {
            pairwiseAlignment(
                cdrh3, target_seq,
                substitutionMatrix = "BLOSUM62",
                gapOpening = -10,
                gapExtension = -0.5
            )
        })

        # Extract alignment scores and percent identities
        scores <- sapply(alignments, score)
        pids <- sapply(alignments, pid)

        # Find the best alignment
        best_idx <- which.max(scores)
        best <- covid[best_idx, ]

        # Return results for this target sequence
        list(
            Sequence = target_seq,
            BestIdentity = sprintf("%.2f%%", pids[best_idx]),
            Name = best$Name,
            BindsTo = best[["Binds.to"]]
        )
    }

    # Apply the best match computation to all target sequences
    results <- lapply(target_seqs, compute_best_match)

    # Convert the list of results to a data frame
    results_df <- do.call(rbind, lapply(results, as.data.frame))
    return(results_df)
}

target_seqs <- clonalCompare(reorganized_list, 
                  top.clones = 10,            
                  samples = c("Pos"), 
                  cloneCall="aa", 
                  exportTable = T) %>%
          mutate(clones = gsub("^C(.*)W_.*", "\\1", clones)) %>%
          mutate(clones = gsub("(.*)_.*", "\\1", clones)) %>%
          filter(nchar(clones) > 6)
target_seqs

result <- similarity(target_seqs$clones)
result

result[order(result$BestIdentity), ]

cdr3 <- FetchData(SNG.US, vars = c("b_cdr3s_aa", 'Best'))
cdr3[grepl("ARDLIDYGMDV", cdr3$b_cdr3s_aa), , drop = FALSE]

cdr3 <- FetchData(SNG.US, vars = c("b_cdr3s_aa", 'Best'))
cdr3[grepl("ARDLMVYGMDV", cdr3$b_cdr3s_aa), , drop = FALSE]

cdr3 <- FetchData(SNG.US, vars = c("b_cdr3s_aa", 'Best'))
cdr3[grepl("ARDRIDYGMDV", cdr3$b_cdr3s_aa), , drop = FALSE]

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
