### Author: Komal S. Rathi

### Purpose

This module performs multi-modal clustering using RNA, CNV, SNV, Methylation and Splicing data on Medulloblastoma samples. 

#### Data version

- OT v12 for histologies (molecular subtype etc), SNV and RNA
- OT v12 for methylation data (I think there was some missing info in v11 which has been fixed in v12)
- [Release 20230309](https://cavatica.sbgenomics.com/u/d3b-bixu-ops/monthly-release-data/files/#q?path=20230309_release) for the CNV gainloss file (this is not available in v11).

### Run Analysis

```
# run full analysis
bash run_analysis.sh
```

### Feature selection

Features were selected from OpenPedCan-analysis v12 datasets using the following filters:

1) CNV:

- The `gainloss.txt` file was first filtered to Medulloblatoma samples.
- Adjusted copy number was calculated. 
- Neutral sites were removed.
- The dataset was reduced to features corresponding to "Gain", "Amplification" in Oncogenes and "Loss", "Complete Loss" in TSGs (using the comprehensive cancer gene list from `annoFuse`)
- Finally, features with `standard deviation >= 0.9` were chosen for clustering. (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176278)

2) SNV:

- The Consensus MAF dataset was first filtered to Medulloblatoma samples.
- Features were reduced to non-synonymous variant classifications e.g. `Missense_Mutation, Frame_Shift_Del, In_Frame_Ins, Frame_Shift_Ins, Splice_Site, Nonsense_Mutation, In_Frame_Del, Nonstop_Mutation, Translation_Start_Site`. 
- Finally, features with `mutation rate > 0.02` (i.e. filter out features with low mutation rate) were chosen for clustering. (Reference: https://www.bioconductor.org/packages/release/bioc/vignettes/iClusterPlus/inst/doc/iManual.pdf)

3) mRNA: 

- The expected counts dataset was first filtered to Medulloblatoma samples. 
- Features were reduced to `Top 1000 most variable protein coding genes` followed by `Rank transformation`.

4) Methylation:

- Methylation beta-values matrix was first filtered to Medulloblatoma samples. 
- Features were reduced to `Top 1000 most variable probes`.

5) Splicing:

- Splice matrix was first filtered to Medulloblatoma samples. 
- Features were reduced to `Top 1000 most variable splice variants`.

### Selection of optimal cluster

From this comment: https://github.com/d3b-center/bixu-tracker/issues/1704#issuecomment-1480304368
For each attempted k, the cluster values were mapped to the samples for each mode of data (e.g. CNV, SNV, Methylation, RNA and Splicing), the fpc stats for each mode were computed, sum of the `average.between` and `average.within` across the modes of data at each k were taken, then the difference of those two summed values at each k was taken to select the optimal cluster number (e.g. the k with the largest difference between the two).

Using this method, the samples were classified into `9 clusters`.

### Input data matrices

```
input
├── cnv_data.tsv # cnv data used as input
├── methyl_data.tsv # methylation data used as input
├── norm_counts.tsv # expression data used as input
├── snv_data.tsv # snv data used as input
├── splice_data.tsv # splice data used as input
└── samples_map.tsv # biospecimens + cohort identifiers for samples used for each modality 

```

### Results

```
results
├── intnmf_fit_all.rds # output of nmf.mnnals for all k values
├── intnmf_clusterstats.tsv # cluster stats across all k-values for each modality
├── intnmf_best_fit.rds # output of nmf.mnnals for best fit (selected k)
├── intnmf_clusters.tsv # output clusters with per sample along with corresponding molecular subtype
├── ari_intnmf_vs_subtypes.txt # adjusted rand index between Intnmf clusters and RNA/Methylation derived molecular subtypes
└── chisq_intnmf_vs_subtypes.txt # chisq test of independence between Intnmf clusters and RNA/Methylation derived molecular subtypes
```

### Plots

```
plots
├── opt_k_cpi_ami.pdf # AMI values across k values
├── opt_k_cpi_ari.pdf # ARI values across k values
├── intnmf_consensus_plot.pdf # consensus plot of optimal k
├── intnmf_silhouette_plot.pdf # silhouette plot of optimal k
├── intnmf_clusters_vs_dkfz_v11_methylation_subclass_balloonplot.pdf # balloon plot of Intnmf clusters vs v11_methylation_subclass
├── intnmf_clusters_vs_dkfz_v11_methylation_subclass_corrplot.pdf # corrplot of Intnmf clusters vs v11_methylation_subclass
├── intnmf_clusters_vs_dkfz_v12_methylation_subclass_balloonplot.pdf # balloon plot of Intnmf clusters vs v12_methylation_subclass
├── intnmf_clusters_vs_dkfz_v12_methylation_subclass_corrplot.pdf # corrplot of Intnmf clusters vs v12_methylation_subclass
├── intnmf_clusters_vs_molsubtype_balloonplot.pdf # balloon plot of Intnmf clusters vs RNA-derived molecular subtypes
├── intnmf_clusters_vs_molsubtype_corrplot.pdf # corrplot of Intnmf clusters vs RNA-derived molecular subtypes
├── feature_level_heatmaps.pdf # feature level heatmaps for all modalities
├── sample_level_heatmaps.pdf # sample level heatmaps for all modalities
└── survival
    ├── survival_IntNMF_clusters.pdf # survival curves of IntNMF clusters
    ├── survival_methyl_subtype_v11.pdf # survival curves of v11_methylation_subclass
    ├── survival_methyl_subtype_v12.pdf # survival curves of v12_methylation_subclass
    └── survival_molsubtype.pdf # # survival curves of RNA-derived molecular subtypes
```

### Association between IntNMF clusters and RNA/Methyl-derived molecular subtypes

Adjusted rand indices for all three comparisons:

```
[1] "IntNMF clusters vs RNA-derived molecular_subtypes"
[1] 0.3373173
Chisq p-value < 2.2e-16

[1] "IntNMF clusters vs dkfz_v11_methylation_subclass"
[1] 0.3426125
Chisq p-value < 2.2e-16

[1] "IntNMF clusters vs dkfz_v12_methylation_subclass"
[1] 0.5166781
Chisq p-value < 2.2e-16
```

### Feature level analysis

Expression, Methylation and Splice data: Enrichment of features with non-zero NMF weights per cluster

```
plots/feature_analysis
├── expression_fgsea.pdf # expression fgsea barplots
├── methyl_fgsea.pdf # methylation fgsea barplots
├── splice_barplot.pdf # splice barplot of enriched pathways
├── splice_barplot_graphics.pdf # splice barplot using the Graphics package (same as the other barplot)
├── splice_goplot.pdf # splice goplot
├── splice_gsoa.pdf # splice gsoa
└── splice_wordcloud.pdf # splice wordcloue
```

Mutation and Copy-number: Top 10 feature (by IntNMF weight) per cluster mapped to Group3/4 subgroups comparison to drivers listed in Fig. 8 of PMID: 31076851. 

```
results/feature_analysis
├── cnv_top10_feature_group34_comparison.tsv
└── mutation_top10_feature_group34_comparison.tsv
```
