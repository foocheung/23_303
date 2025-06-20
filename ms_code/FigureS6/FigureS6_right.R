# --- Required Libraries ---
library(Seurat)
library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)
library(Matrix)
library(ggpubr)

# --- Helper function ---
get_matched_umi <- function(ag6_df, bcr_path, label, prefix = "") {
  bcr <- read_csv(bcr_path, show_col_types = FALSE)
  ag6_df %>%
    filter(barcode %in% paste0(prefix, bcr$Clean_Clone_ID)) %>%
    mutate(Group = label)
}

# ==========================
# --- BEAM Main Data Load ---
# ==========================
beam_data <- Read10X_h5("../../data/multi_config_10_sample_filtered_feature_bc_matrix.h5")
seurat_beam <- CreateSeuratObject(counts = beam_data[["Gene Expression"]])
seurat_beam[["ADT"]] <- CreateAssayObject(counts = beam_data[["Antibody Capture"]])
seurat_beam[["Antigen"]] <- CreateAssayObject(counts = beam_data[["Antigen Capture"]])
seurat_beam[["percent.mt"]] <- PercentageFeatureSet(seurat_beam, pattern = "^MT-")
seurat_beam <- subset(seurat_beam, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & percent.mt < 10)

# Light chain
DefaultAssay(seurat_beam) <- "Antigen"
beam_ag6 <- data.frame(barcode = colnames(seurat_beam), ag6_umi = GetAssayData(seurat_beam)["ag6", ])

beam_df_light <- bind_rows(
  get_matched_umi(beam_ag6, "../../data/intersection_sets_light/intersection_6_BEAM_Positive_AND_DualLabeling_Positive.csv", "BEAM_AND_Dual"),
  get_matched_umi(beam_ag6, "../../data/intersection_sets_light/intersection_1_BEAM_Positive.csv", "BEAM_Only")
)

# Heavy chain
beam_d614g <- data.frame(barcode = colnames(seurat_beam), ag6_umi = GetAssayData(seurat_beam)["ag6", ])
beam_df_heavy <- bind_rows(
  get_matched_umi(beam_d614g, "../../data/intersection_sets/intersection_6_BEAM_Positive_AND_DualLabeling_Positive.csv", "BEAM_AND_Dual"),
  get_matched_umi(beam_d614g, "../../data/intersection_sets/intersection_1_BEAM_Positive.csv", "BEAM_Only")
)

# ==========================
# --- BEAM Negative Data Load ---
# ==========================
beam_neg_data <- Read10X_h5("../../data/multi_config_11_sample_filtered_feature_bc_matrix.h5")
seurat_beam_neg <- CreateSeuratObject(counts = beam_neg_data[["Gene Expression"]])
seurat_beam_neg[["ADT"]] <- CreateAssayObject(counts = beam_neg_data[["Antibody Capture"]])
seurat_beam_neg[["Antigen"]] <- CreateAssayObject(counts = beam_neg_data[["Antigen Capture"]])
seurat_beam_neg[["percent.mt"]] <- PercentageFeatureSet(seurat_beam_neg, pattern = "^MT-")
seurat_beam_neg <- subset(seurat_beam_neg, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & percent.mt < 10)

# Light chain
DefaultAssay(seurat_beam_neg) <- "Antigen"
beam_neg_ag6 <- data.frame(barcode = colnames(seurat_beam_neg), ag6_umi = GetAssayData(seurat_beam_neg)["ag6", ])
beam_df_light <- bind_rows(
  beam_df_light,
  get_matched_umi(beam_neg_ag6, "../../data/intersection_sets_light/intersection_2_BEAM_Negative.csv", "BEAM_Neg_Only")
)

# Heavy chain
beam_neg_d614g <- data.frame(barcode = colnames(seurat_beam_neg), ag6_umi = GetAssayData(seurat_beam_neg)["ag6", ])
beam_df_heavy <- bind_rows(
  beam_df_heavy,
  get_matched_umi(beam_neg_d614g, "../../data/intersection_sets/intersection_2_BEAM_Negative.csv", "BEAM_Neg_Only")
)

# ==========================
# --- DUAL Data Load ---
# ==========================
dual1 <- Read10X_h5("../../data/multi_config3_matrix.h5")
dual2 <- Read10X_h5("../../data/multi_config4_matrix.h5")

rna1 <- dual1[["Gene Expression"]]
rna2 <- dual2[["Gene Expression"]]
colnames(rna1) <- paste0("MC3_", colnames(rna1))
colnames(rna2) <- paste0("MC4_", colnames(rna2))
dual_rna <- cbind(rna1, rna2)

adt1 <- dual1[["Antibody Capture"]]
adt2 <- dual2[["Antibody Capture"]]
colnames(adt1) <- paste0("MC3_", colnames(adt1))
colnames(adt2) <- paste0("MC4_", colnames(adt2))
dual_adt <- cbind(adt1, adt2)

seurat_dual <- CreateSeuratObject(counts = dual_rna)
seurat_dual[["ADT"]] <- CreateAssayObject(counts = dual_adt)
seurat_dual[["percent.mt"]] <- PercentageFeatureSet(seurat_dual, pattern = "^MT-")
seurat_dual <- subset(seurat_dual, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & percent.mt < 10)

# Light chain
DefaultAssay(seurat_dual) <- "ADT"
dual_ag6 <- data.frame(barcode = colnames(seurat_dual), ag6_umi = GetAssayData(seurat_dual)["D614G", ])
dual_df_light <- bind_rows(
  get_matched_umi(dual_ag6, "../../data/intersection_sets_light/intersection_6_BEAM_Positive_AND_DualLabeling_Positive.csv", "BEAM_AND_Dual", prefix = "MC3_"),
  get_matched_umi(dual_ag6, "../../data/intersection_sets_light/intersection_3_DualLabeling_Positive.csv", "Dual_Only", prefix = "MC3_")
)

# Heavy chain
dual_d614g <- data.frame(barcode = colnames(seurat_dual), ag6_umi = GetAssayData(seurat_dual)["D614G", ])
dual_df_heavy <- bind_rows(
  get_matched_umi(dual_d614g, "../../data/intersection_sets/intersection_6_BEAM_Positive_AND_DualLabeling_Positive.csv", "BEAM_AND_Dual", prefix = "MC3_"),
  get_matched_umi(dual_d614g, "../../data/intersection_sets/intersection_3_DualLabeling_Positive.csv", "Dual_Only", prefix = "MC3_")
)

# ==========================
# --- DUAL Negative Data Load ---
# ==========================
dual_neg_data <- Read10X_h5("../../data/7multi_config15_matrix.h5")
dual_neg_rna <- dual_neg_data[["Gene Expression"]]
dual_neg_adt <- dual_neg_data[["Antibody Capture"]]
colnames(dual_neg_rna) <- paste0("MC3_", colnames(dual_neg_rna))
colnames(dual_neg_adt) <- paste0("MC3_", colnames(dual_neg_adt))

seurat_dual_neg <- CreateSeuratObject(counts = dual_neg_rna)
seurat_dual_neg[["ADT"]] <- CreateAssayObject(counts = dual_neg_adt)
seurat_dual_neg[["percent.mt"]] <- PercentageFeatureSet(seurat_dual_neg, pattern = "^MT-")
seurat_dual_neg <- subset(seurat_dual_neg, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & percent.mt < 10)

# Light chain
DefaultAssay(seurat_dual_neg) <- "ADT"
dual_neg_ag6 <- data.frame(barcode = colnames(seurat_dual_neg), ag6_umi = GetAssayData(seurat_dual_neg)["D614G", ])
dual_df_light <- bind_rows(
  dual_df_light,
  get_matched_umi(dual_neg_ag6, "../../data/intersection_sets_light/intersection_4_DualLabeling_Negative.csv", "Dual_Neg_Only", prefix = "MC3_")
)

# Heavy chain
dual_neg_d614g <- data.frame(barcode = colnames(seurat_dual_neg), ag6_umi = GetAssayData(seurat_dual_neg)["D614G", ])
dual_df_heavy <- bind_rows(
  dual_df_heavy,
  get_matched_umi(dual_neg_d614g, "../../data/intersection_sets/intersection_4_DualLabeling_Negative.csv", "Dual_Neg_Only", prefix = "MC3_")
)

# ==========================
# --- Plotting Functions ---
# ==========================
make_log_plot <- function(df, title, comparisons = NULL) {
  p <- ggplot(df, aes(x = Group, y = ag6_umi, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
    scale_y_log10() +
    labs(title = title, y = "UMI (log10)", x = "") +
    theme_bw() + theme(legend.position = "none")
  
  if (!is.null(comparisons)) {
    p <- p + stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif")
  }
  
  return(p)
}

# ==========================
# --- Generate & Save Plots ---
# ==========================
beam_comparisons <- list(
  c("BEAM_AND_Dual", "BEAM_Only"),
  c("BEAM_AND_Dual", "BEAM_Neg_Only"),
  c("BEAM_Only", "BEAM_Neg_Only")
)

dual_comparisons <- list(
  c("BEAM_AND_Dual", "Dual_Only"),
  c("BEAM_AND_Dual", "Dual_Neg_Only"),
  c("Dual_Only", "Dual_Neg_Only")
)

p1 <- make_log_plot(beam_df_light, "BEAM Light Chain: D614G UMI", beam_comparisons)
p2 <- make_log_plot(beam_df_heavy, "BEAM Heavy Chain: D614G UMI", beam_comparisons)
p3 <- make_log_plot(dual_df_light, "DUAL Light Chain: D614G UMI", dual_comparisons)
p4 <- make_log_plot(dual_df_heavy, "DUAL Heavy Chain: D614G UMI", dual_comparisons)



ggsave("FigureS6_right.pdf")

