# Bulk-RNA-seq-data-analysis_mRNAsLGM

# Bulk RNA-seq pipeline

Standard bulk RNA-seq workflow:
- Trim Galore (adapter + quality trimming)
- Bowtie2 alignment to hg38
- featureCounts gene quantification
- limma-voom differential expression + heatmaps/QC

## Files
- `workflow/Snakefile` : pipeline
- `config/config.yaml` : paths/parameters
- `config/samples.tsv` : sample groups/replicates
- `R/02_limma_voom_DE_heatmap.R` : DE + plots
- `envs/conda.yaml` : environment (optional)

## Note
Raw FASTQ/BAM files are not included in this repository.
