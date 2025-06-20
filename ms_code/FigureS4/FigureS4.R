library(dplyr)
library(readr)
library(tidyr)
library(ComplexHeatmap)
library(ggplot2)
library(grid)
library(gridExtra)
library(stringdist)
library(ape)
library(igraph)

# -------- Tree builder function ---------
make_colored_tree_by_edit_distance <- function(clone_data, label, outdir = ".", cache_prefix = "tree_cache", force_rebuild = TRUE) {
  message("Building or loading tree for: ", label)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(outdir, paste0(cache_prefix, "_", gsub("[^A-Za-z0-9]", "_", label), ".rds"))
  
  group_colors <- c("BEAM_Only" = "steelblue", "Dual_Only" = "darkorange", "Shared" = "purple")
  
  if (file.exists(cache_file) && !force_rebuild) {
    message("Loading cached tree from: ", cache_file)
    cache_data <- readRDS(cache_file)
    tree <- cache_data$tree
    tip_colors <- cache_data$tip_colors
    tip_cdr3s <- cache_data$tip_cdr3s
  } else {
    message("Building new tree and computing clusters...")
    clone_data <- clone_data %>% filter(!duplicated(junction_aa))
    if (nrow(clone_data) < 2) return(NULL)
    
    aa_seqs <- setNames(clone_data$junction_aa, clone_data$Clone_ID)
    n <- length(aa_seqs)
    
    dist_matrix <- matrix(1, n, n)
    rownames(dist_matrix) <- colnames(dist_matrix) <- names(aa_seqs)
    
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        pid <- 1 - stringdist::stringdist(aa_seqs[i], aa_seqs[j], method = "lv") / max(nchar(aa_seqs[i]), nchar(aa_seqs[j]))
        coverage <- min(nchar(aa_seqs[i]), nchar(aa_seqs[j])) / max(nchar(aa_seqs[i]), nchar(aa_seqs[j]))
        if (pid >= 0.9 && coverage >= 0.9) dist_matrix[i, j] <- 0
      }
    }
    
    tree <- nj(as.dist(stringdistmatrix(aa_seqs, aa_seqs, method = "lv")))
    tree$tip.label <- names(aa_seqs)[as.integer(tree$tip.label)]
    tip_cdr3s <- aa_seqs[tree$tip.label]
    
    graph <- graph_from_adjacency_matrix(dist_matrix == 0, mode = "undirected", diag = FALSE)
    cluster_membership <- components(graph)$membership
    clone_data$Cluster <- paste0("Cluster_", cluster_membership[match(clone_data$Clone_ID, names(aa_seqs))])
    
    cluster_counts <- table(clone_data$Cluster)
    clone_data <- clone_data %>%
      mutate(Isolated = Cluster %in% names(cluster_counts[cluster_counts == 1]))
    
    label_map <- clone_data %>%
      filter(junction_aa %in% tip_cdr3s) %>%
      distinct(junction_aa, Group)
    tip_colors <- group_colors[label_map$Group[match(tip_cdr3s, label_map$junction_aa)]]
    
    saveRDS(list(tree = tree, tip_colors = tip_colors, tip_cdr3s = tip_cdr3s), file = cache_file)
    message("Tree cached to: ", cache_file)
  }
  
  # ── Save CDR3-labeled tree ──
  labeled_pdf <- "FigureS4_CDR3Labels.pdf"
  tree_labeled <- tree
  tree_labeled$tip.label <- tip_cdr3s
  
  pdf(file = labeled_pdf, width = 30, height = 26)
  plot(tree_labeled, show.tip.label = FALSE, main = paste("CDR3-Labeled Tree -", label), cex = 0.6)
  tree_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  text(
    x = tree_plot$xx[1:tree_plot$Ntip],
    y = tree_plot$yy[1:tree_plot$Ntip],
    labels = tree_labeled$tip.label,
    col = tip_colors,
    pos = 4,
    offset = 0.5,
    cex = 0.6
  )
  legend("topright", legend = names(group_colors), col = group_colors, pch = 19, cex = 0.8)
  dev.off()
  message("CDR3-labeled tree saved to: ", labeled_pdf)
  
  # ── Save colored-tip-only tree ──
  color_pdf <-  "FigureS4.pdf"
  pdf(file = color_pdf, width = 30, height = 26)
  plot(tree, show.tip.label = FALSE, main = paste("Colored CDR3 Tree -", label), cex = 0.6)
  tiplabels(pch = 19, col = tip_colors, cex = 0.6)
  legend("topright", legend = names(group_colors), col = group_colors, pch = 19, cex = 0.8)
  dev.off()
  message("Tree plot saved to: ", color_pdf)
  

}

# --------- Wrapper function ---------
process_upset_merge_tree <- function(matrix_csv, label_tag, test_mode = FALSE) {
  df <- read_csv(matrix_csv, col_types = cols()) %>%
    mutate(code = paste0(BEAM_Positive, BEAM_Negative, DualLabeling_Positive, DualLabeling_Negative))
  
  beam_only <- df %>% filter(code == "1000") %>% mutate(Group = "BEAM_Only")
  dual_only <- df %>% filter(code == "0010") %>% mutate(Group = "Dual_Only")
  shared    <- df %>% filter(code == "1010") %>% mutate(Group = "Shared")
  
  merged <- bind_rows(beam_only, dual_only, shared) %>%
    filter(!duplicated(junction_aa)) %>%
    mutate(Clone_ID = paste0("Clone", row_number()))
  
  if (test_mode) {
    set.seed(123)
    merged <- merged %>% sample_frac(0.2)
    message("Test mode: sampled ", nrow(merged), " rows.")
  }
  
  out_dir <- dirname(matrix_csv)
  make_colored_tree_by_edit_distance(merged, paste0(label_tag, "_v3_Merged"), outdir = out_dir, force_rebuild = TRUE)
}

# -------- Process individual set ---------
dirs <- "../../data"

for (dir in dirs) {
  for (type in c("Light")) {
    matrix_file <- list.files(dir, pattern = paste0("JunctionAA_Matrix_", type, "\\.csv"), full.names = TRUE)
    if (length(matrix_file) != 1) next
    label_tag <- paste0(gsub("UpSet_Exported_Sets_", "", basename(dir)), "_", type)
    message("\nMerging sequences for tree: ", label_tag)
    process_upset_merge_tree(matrix_file, label_tag)
  }
}
