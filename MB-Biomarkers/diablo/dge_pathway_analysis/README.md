### Differential Gene Expression and Pathway Analysis

#### Author: Komal S. Rathi

#### Description of scripts
***
`01-dge-analysis-deseq.R`:  The function of this script is to perform differential gene expression using `DESeq2` using the `cluster vs rest` approach. 

`01-dge-analysis-noiseq.R`: The function of this script is to perform differential gene expression using `NOISeq` using the `cluster vs rest` approach.

#### Inputs

```
../../data
├── c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt # KEGG MEDICUS gmt file
├── gencode.v39.primary_assembly.annotation.gtf.gz # gencode v39
└── gene-counts-rsem-expected_count_medullo-collapsed.rds

# intNMF inputs but any other clustering method can be used
../intNMF
├── input/samples_map.tsv # sample_id to biospecimen mapping  
└── results/intnmf_clusters.tsv # intNMF clusters
```

#### Outputs

For all the `cluster-vs-rest` comparisons, the DESeq2/NOISeq output with FDR adjusted p-value < 0.05 is saved under `diffexpr_output_per_cluster.tsv` and the pathway enrichment output using `clusterProfiler::GSEA` with FDR adjusted p-value < 0.05 is saved under `_gsea.tsv`.


```
# DESeq2 output
results/intNMF/deseq
├── diffexpr_output_per_cluster.tsv
├── kegg_medicus
│   └── cluster_{n}_vs_rest_gsea.tsv
└── reactome
    └── cluster_{n}_vs_rest_gsea.tsv

# NOISeq output
results/intNMF/noiseq
├── diffexpr_output_per_cluster.tsv
├── kegg_medicus
│   └── cluster_{n}_vs_rest_gsea.tsv
└── reactome
    └── cluster_{n}_vs_rest_gsea.tsv
```

For each `cluster-vs-rest` enrichment, a barplot of top 10 upregulated and top 10 downregulated (i.e. maximum of 20 pathways) identified by `clusterProfiler::GSEA` utilizing `KEGG MEDICUS` and `REACTOME` pathways at `FDR < 0.05` is generated under `*_gsea_barplot.pdf`, dotplot is generated under `*_gsea_dotplot.pdf` and network under `*_gsea_cnet.pdf`. 


When using NOISeq for differential expression analysis, `noiseq_pca.pdf` is generated which has the PCA plot of input matrix before and after NOISeq batch correction. 

```
# here n is cluster number of interest

# DESeq2 output
plots/intNMF/deseq
├── kegg_medicus
│   ├── cluster_{n}_vs_rest_gsea_barplot.pdf
│   ├── cluster_{n}_vs_rest_gsea_cnet.pdf
│   └── cluster_{n}_vs_rest_gsea_dotplot.pdf
└── reactome
    ├── cluster_{n}_vs_rest_gsea_barplot.pdf
    ├── cluster_{n}_vs_rest_gsea_cnet.pdf
    └── cluster_{n}_vs_rest_gsea_dotplot.pdf


# NOISeq output
plots/intNMF/noiseq
├── kegg_medicus
│   ├── cluster_{n}_vs_rest_gsea_barplot.pdf
│   ├── cluster_{n}_vs_rest_gsea_cnet.pdf
│   └── cluster_{n}_vs_rest_gsea_dotplot.pdf
├── noiseq_pca.pdf
└── reactome
    ├── cluster_{n}_vs_rest_gsea_barplot.pdf
    ├── cluster_{n}_vs_rest_gsea_cnet.pdf
    └── cluster_{n}_vs_rest_gsea_dotplot.pdf

```

***

### Run analysis

To run the full analysis, use the bash script as follows:

```
bash run_analysis.sh
```
