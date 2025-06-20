library(readr)
library(dplyr)
library(stringr)

# === List all UMI metadata files ===
input_files <- list.files("../../data", pattern = "BEAM_(Heavy|Light)_UMI_Metadata.csv$", full.names = TRUE)

for (file in input_files) {
  cat("\U0001F4C2 Processing:", file, "\n")
  flush.console()
  
  base <- str_remove(basename(file), "_UMI_Metadata.csv")
  csv_path <- paste0("string_matched_", base, "_to_CoVAbDab_details.csv")
  
  # Read and clean BCR metadata
  bcr <- read_csv(file, show_col_types = FALSE) %>%
    mutate(junction_aa = trimws(toupper(junction_aa))) %>%
    filter(Group != "BEAM_Neg_Only")
  
  bcr_summary <- bcr %>%
    group_by(junction_aa) %>%
    mutate(clone_count = n()) %>%
    ungroup()
  
  unique_junctions <- unique(bcr_summary$junction_aa)
  fasta_lines <- c()
  detail_rows <- list()
  match_counter <- 0
  
  for (jaa in unique_junctions) {
    bcr_row <- bcr_summary %>% filter(junction_aa == jaa) %>% slice(1)
    
    grep_cmd <- paste("grep", jaa, "../../data/CoV-AbDab_080224.csv")
    covab_hits <- system(grep_cmd, intern = TRUE)
    
    cat("\U0001F50D Checking:", jaa, "| Matches found:", length(covab_hits), "\n")
    flush.console()
    
    if (length(covab_hits) == 0) next
    
    for (line in covab_hits) {
      covab_id <- str_extract(line, "^[^,;'\\\"]+")
      covab_source <- str_match(line, "(?:[^,]*,){19}([^,]*)")
      ag6_val <- as.character(bcr_row$ag6_umi)
      ag6_val[is.na(ag6_val)] <- "NA"
      
      header <- paste0(">copies=", bcr_row$clone_count,
                       ";umi_count=", bcr_row$umi_count,
                       ";group=", bcr_row$Group,
                       ";ag6_umi=", ag6_val,
                       ";covab=", covab_id)
      fasta_seq <- paste0(header, "\n", jaa)
      fasta_lines <- c(fasta_lines, fasta_seq)
      
      detail_rows[[length(detail_rows) + 1]] <- data.frame(
        junction_aa = jaa,
        umi_count = bcr_row$umi_count,
        group = bcr_row$Group,
        ag6_umi = ag6_val,
        clone_count = bcr_row$clone_count,
        match_names = covab_id,
        match_sources = covab_source,
        stringsAsFactors = FALSE
      )
      
      match_counter <- match_counter + 1
    }
  }
  
  if (length(detail_rows) > 0) {
    write_csv(bind_rows(detail_rows), csv_path)
    cat("\U0001F4C4 Match details written:", csv_path, "\n")
  } else {
    cat("\u26A0\uFE0F No CoV-AbDab hits found for", base, "\n")
  }
  
  cat("\u2705 Done:", base,
      "| Total CoV-AbDab hits:", match_counter, "\n\n")
  flush.console()
}
