# ------------------------------------------------------------------------------
# Title: covabdab-clone-scatter.R
# Description: analyzes BEAM heavy and light chain sequences, identifies hits 
# against CoV-AbDab, and visualizes clone_count vs ag6_umi with hit highlighting and log scaling
# Created by: Foo Cheung (foo.cheung@nih.gov)
# Date: 2025
# ------------------------------------------------------------------------------


library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggrepel)

# === Function to extract AA sequences (every 2nd line) ===
extract_junction_aa <- function(x) {
  x <- x$copies
  aa_lines <- x[seq(2, length(x), 2)]  # lines 2, 4, 6, ...
  return(aa_lines)
}

# === Parse beam_heavy and beam_light FASTA files ===
beam_heavy_raw <- read_csv("../../data/v2_grep_matched_v3_BEAM_Heavy_to_CoVAbDab.fasta", col_names = "copies")
beam_light_raw <- read_csv("../../data/v2_grep_matched_v3_BEAM_Light_to_CoVAbDab.fasta", col_names = "copies")

beam_heavy_hits <- extract_junction_aa(beam_heavy_raw)
beam_light_hits <- extract_junction_aa(beam_light_raw)

# === Load metadata and filter groups ===
heavy_df <- read_csv("../../data/v3_BEAM_Heavy_UMI_Metadata.csv", show_col_types = FALSE) %>%
  filter(Group %in% c("BEAM_AND_Dual", "BEAM_Only"))
light_df <- read_csv("../../data/v3_BEAM_light_UMI_Metadata.csv", show_col_types = FALSE) %>%
  filter(Group %in% c("BEAM_AND_Dual", "BEAM_Only"))

# === Summarize data and add hit label (using raw values) ===
heavy_summary <- heavy_df %>%
  group_by(junction_aa) %>%
  summarise(
    ag6_umi = first(ag6_umi),
    umi_count = first(umi_count),
    Group = first(Group),
    clone_count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    is_hit = if_else(junction_aa %in% beam_heavy_hits, "Hit", "No Hit")
  )

light_summary <- light_df %>%
  group_by(junction_aa) %>%
  summarise(
    ag6_umi = first(ag6_umi),
    umi_count = first(umi_count),
    Group = first(Group),
    clone_count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    is_hit = if_else(junction_aa %in% beam_light_hits, "Hit", "No Hit")
  )



library(ggplot2)
library(gridExtra)
library(scales)
library(ggplot2)
library(gridExtra)
library(scales)

# === Adjusted values ===
heavy_summary <- heavy_summary %>%
  mutate(ag6_umi_adj = ag6_umi + 1)

light_summary <- light_summary %>%
  mutate(ag6_umi_adj = ag6_umi + 1)

# === Determine range for x-axis ticks ===
x_min_h <- min(heavy_summary$clone_count, na.rm = TRUE)
x_max_h <- max(heavy_summary$clone_count, na.rm = TRUE)

x_min_l <- min(light_summary$clone_count, na.rm = TRUE)
x_max_l <- max(light_summary$clone_count, na.rm = TRUE)

# === Heavy Plot ===
p1 <- ggplot(heavy_summary, aes(x = clone_count, y = ag6_umi_adj)) +
  geom_point(data = subset(heavy_summary, is_hit == "No Hit"),
             color = "gray", alpha = 0.5, size = 1.5) +
  geom_point(data = subset(heavy_summary, is_hit == "Hit"),
             color = "red", alpha = 1, size = 3) +
  scale_x_continuous(breaks = seq(x_min_h, x_max_h, by = 1),
                     labels = label_number(accuracy = 1)) +
  scale_y_continuous(trans = "log10", labels = label_number(accuracy = 1)) +
  labs(title = "Heavy Chain: Clone Count vs log10(ag6_umi + 1)",
       x = "Clone Count (raw, 1-unit steps)", y = "ag6_umi (log10 +1)")

# === Light Plot ===
p2 <- ggplot(light_summary, aes(x = clone_count, y = ag6_umi_adj)) +
  geom_point(data = subset(light_summary, is_hit == "No Hit"),
             color = "gray", alpha = 0.5, size = 1.5) +
  geom_point(data = subset(light_summary, is_hit == "Hit"),
             color = "blue", alpha = 1, size = 3) +
  scale_x_continuous(breaks = seq(x_min_l, x_max_l, by = 1),
                     labels = label_number(accuracy = 1)) +
  scale_y_continuous(trans = "log10", labels = label_number(accuracy = 1)) +
  labs(title = "Light Chain: Clone Count vs log10(ag6_umi + 1)",
       x = "Clone Count (raw, 1-unit steps)", y = "ag6_umi (log10 +1)")

# === Save to png ===
pdf("Figure3E.pdf")
p1
dev.off()
