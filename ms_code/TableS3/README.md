CDR3 Homology Analysis with Known Dengue Antibodies
This analysis was implemented using the code located in the Figure5CDFG/ directory. The best hit showed 78% sequence identity, with the top 5 alignments summarized here.

We extracted 2,486 unique CDRH3 amino acid sequences from B cells in three Dengue-exposed subjects. These sequences were selected based on dual labeling criteria: DV > 25 (Dengue labeling intensity) HSA < 25 (Human Serum Albumin background) This filtering ensured enrichment for potentially antigen-specific cells. To assess potential Dengue reactivity, we compared these CDR3s against a curated database of 288 known Dengue-reactive CDRH3 sequences from the Natsrita et al. study.

The comparison was performed using blastp with the following command:

blastp -query target_cdr3.fa
-db denguedb_ann
-out aln_results_with_annotation.tsv
-evalue 1e-5
-num_threads 8
