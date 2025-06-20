library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)


df<-read.csv("../../data/heavy_cdr3_fasta_with_annotation_hits.csv")

df <- df %>%
  mutate(
    SimplifiedGroup = case_when(
      grepl("BEAM_Only", Group) ~ "BEAM_Only",
      grepl("Dual_Only", Group) ~ "Dual_Only",
      grepl("Shared", Group) ~ "Shared",
      TRUE ~ "Other"
    )
  )

# Step 1: Calculate frequency table
freq_df <- df %>%
  filter(SimplifiedGroup %in% c("BEAM_Only", "Dual_Only", "Shared")) %>%
  mutate(CopyCount = BEAM_Count + Dual_Count) %>%
  filter(CopyCount > 1, CopyCount <= 12) %>%
  group_by(SimplifiedGroup, CopyCount) %>%
  summarise(Count = n(), .groups = "drop")

# Step 2: Add total per group to convert to frequency
freq_df <- freq_df %>%
  group_by(SimplifiedGroup) %>%
  mutate(Frequency = Count / sum(Count))





freq_plot3 <- ggplot(freq_df, aes(x = CopyCount, y = Frequency, color = SimplifiedGroup)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  # geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -1, size = 3) +
  scale_x_continuous(breaks = 2:12) +
  labs(
    title = "Proportion of Clones (Light Chain) at Each Copy Count",
    x = "Clone Copy Count",
    y = "Proportion per Group"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())+ theme(plot.margin = margin(10, 20, 10, 10),
                                               clip = "off")

df<-read.csv("../../data/heavy_cdr3_fasta_with_annotation_hits.csv")

df <- df %>%
  mutate(
    SimplifiedGroup = case_when(
      grepl("BEAM_Only", Group) ~ "BEAM_Only",
      grepl("Dual_Only", Group) ~ "Dual_Only",
      grepl("Shared", Group) ~ "Shared",
      TRUE ~ "Other"
    )
  )

# Step 1: Calculate frequency table
freq_df <- df %>%
  filter(SimplifiedGroup %in% c("BEAM_Only", "Dual_Only", "Shared")) %>%
  mutate(CopyCount = BEAM_Count + Dual_Count) %>%
  filter(CopyCount > 1, CopyCount <= 12) %>%
  group_by(SimplifiedGroup, CopyCount) %>%
  summarise(Count = n(), .groups = "drop")

# Step 2: Add total per group to convert to frequency
freq_df <- freq_df %>%
  group_by(SimplifiedGroup) %>%
  mutate(Frequency = Count / sum(Count))


freq_plot1 <- ggplot(freq_df, aes(x = CopyCount, y = Frequency, color = SimplifiedGroup)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  # geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -1, size = 3) +
  scale_x_continuous(breaks = 2:12) +
  labs(
    title = "Proportion of Clones (Heavy Chain) at Each Copy Count",
    x = "Clone Copy Count",
    y = "Proportion per Group"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())+ theme(plot.margin = margin(10, 20, 10, 10),
                                               clip = "off") #+ geom_text(aes(label = Count),
#          position = position_dodge(width = 0.9),
#           vjust = -0.2, size = 3)

pdf("FigureS6_center.pdf", width = 15)
freq_plot3
dev.off()