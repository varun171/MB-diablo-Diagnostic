# run DESeq2 analysis
Rscript --vanilla 01-dge-analysis-deseq.R \
--expr_mat "../../data/gene-counts-rsem-expected_count_medullo-collapsed.rds" \
--gtf_file "../../data/gencode.v39.primary_assembly.annotation.gtf.gz" \
--cluster_file "../intNMF/results/intnmf_clusters.tsv" \
--sample_map_file "../intNMF/input/samples_map.tsv" \
--kegg_medicus_file "../../data/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt" \
--results_dir "results/intNMF/deseq" \
--plots_dir "plots/intNMF/deseq"

# run NOISEQ analysis
Rscript --vanilla 01-dge-analysis-noiseq.R \
--expr_mat "../../data/gene-counts-rsem-expected_count_medullo-collapsed.rds" \
--gtf_file "../../data/gencode.v39.primary_assembly.annotation.gtf.gz" \
--cluster_file "../intNMF/results/intnmf_clusters.tsv" \
--sample_map_file "../intNMF/input/samples_map.tsv" \
--kegg_medicus_file "../../data/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt" \
--results_dir "results/intNMF/noiseq" \
--plots_dir "plots/intNMF/noiseq"


