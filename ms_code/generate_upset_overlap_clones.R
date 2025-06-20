# --- Required Libraries ---
library(dplyr)
library(readr)
library(tidyr)

# -------- Paths to Source TSVs ---------
beam_pos_path <- "./s_10/filtered_contig_db-pass.tsv"
beam_neg_path <- "./s_11/filtered_contig_db-pass.tsv"
dual_pos1_path <- "./S3_multi_config_3.joined.tsv"
dual_pos2_path <- "./S3_multi_config_4.joined.tsv"
dual_neg_path  <- "./S3_7multi_config_15.joined.tsv"

# -------- Read function for contig/joined files ---------
read_clone_data <- function(file_path, dataset_name, filter_guess = NULL, exclude = FALSE, has_best_guess = FALSE) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }

  df <- read_tsv(file_path, col_types = cols())
  if (has_best_guess) {
    df <- df %>% filter(DROPLET.TYPE == "SNG")
    if (!is.null(filter_guess)) {
      if (exclude) {
        df <- df %>% filter(!grepl(filter_guess, BEST.GUESS))
      } else {
        df <- df %>% filter(grepl(filter_guess, BEST.GUESS))
      }
    }
  }

  df %>%
    select(sequence_id, junction_aa, v_call) %>%
    mutate(Dataset = dataset_name) %>%
    distinct()
}

# -------- Load all contig/joined data for chain mapping ---------
load_all_chain_map <- function() {
  beam_pos <- read_clone_data(beam_pos_path, "BEAM_Positive")
  beam_neg <- read_clone_data(beam_neg_path, "BEAM_Negative")
  dual_pos1 <- read_clone_data(dual_pos1_path, "DualLabeling_Positive1", has_best_guess = TRUE)
  dual_pos2 <- read_clone_data(dual_pos2_path, "DualLabeling_Positive2", has_best_guess = TRUE)
  dual_neg <- read_clone_data(dual_neg_path, "DualLabeling_Negative", has_best_guess = TRUE)

  bind_rows(beam_pos, beam_neg, dual_pos1, dual_pos2, dual_neg) %>%
    mutate(Chain_Type = case_when(
      grepl("IGH|TRB", v_call, ignore.case = TRUE) ~ "Heavy",
      grepl("IG[KL]|TRA|TRG", v_call, ignore.case = TRUE) ~ "Light",
      TRUE ~ "Unknown"
    )) %>%
    select(junction_aa, Chain_Type) %>%
    distinct()
}

# -------- Function to write FASTA file --------
write_junction_fasta <- function(df, file_path) {
  fasta_lines <- df %>%
    filter(!is.na(junction_aa)) %>%
    mutate(
      header = paste0(
        ">Clone_", Clone_ID,
        " | Chain=", Chain_Type,
        " | Group=", Group
      ),
      sequence = junction_aa
    ) %>%
    select(header, sequence) %>%
    pivot_longer(cols = everything(), values_to = "line") %>%
    pull(line)

  writeLines(fasta_lines, file_path)
  return(nrow(df))
}

# -------- Process matrix to FASTA ---------
process_fasta_only <- function(matrix_csv, label_tag, chain_map) {
  df <- read_csv(matrix_csv, col_types = cols()) %>%
    mutate(code = paste0(BEAM_Positive, BEAM_Negative, DualLabeling_Positive, DualLabeling_Negative))

  beam_only <- df %>% filter(code == "1000") %>% mutate(Group = "BEAM_Only")
  dual_only <- df %>% filter(code == "0010") %>% mutate(Group = "Dual_Only")
  shared    <- df %>% filter(code == "1010") %>% mutate(Group = "Shared")

  merged <- bind_rows(beam_only, dual_only, shared) %>%
    filter(!duplicated(junction_aa)) %>%
    mutate(Clone_ID = paste0(row_number())) %>%
    left_join(chain_map, by = "junction_aa")

  out_dir <- dirname(matrix_csv)
  fasta_path <- file.path(out_dir, paste0(label_tag, "_JunctionAA.fasta"))
  n_written <- write_junction_fasta(merged, fasta_path)

  return(tibble(Matrix = label_tag, Sequences = n_written))
}

# -------- Main FASTA-only loop ---------
chain_map <- load_all_chain_map()
dirs <- list.dirs(".", recursive = FALSE, full.names = TRUE)
dirs <- dirs[grepl("UpSet_Exported_Sets_PosNegFiltered", dirs)]

summary_table <- list()

for (dir in dirs) {
  for (type in c("All", "Heavy", "Light")) {
    matrix_file <- list.files(dir, pattern = paste0("JunctionAA_Matrix_", type, "\\.csv"), full.names = TRUE)
    if (length(matrix_file) != 1) next
    label_tag <- paste0(gsub("UpSet_Exported_Sets_", "", basename(dir)), "_", type)
    message("Generating FASTA for: ", label_tag)
    summary_table[[label_tag]] <- process_fasta_only(matrix_file, label_tag, chain_map)
  }
}

# -------- Summary of written FASTA sequences ---------
summary_df <- bind_rows(summary_table)
write_csv(summary_df, "JunctionAA_FASTA_Sequence_Counts.csv")
print(summary_df)

# Save annotated matrix to CSV for visual inspection
write_csv(data_all,  "F2_Annotated_All_Data.csv")
write_csv(data_pos,  "F2_Annotated_PosFiltered.csv")
write_csv(data_both, "F2_Annotated_PosNegFiltered.csv")
