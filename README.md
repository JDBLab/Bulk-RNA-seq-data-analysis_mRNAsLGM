# Bulk RNA-seq pipeline (Trim Galore → Bowtie2 → featureCounts → limma-voom)

This repository contains a **standard bulk RNA-seq analysis workflow** used for paired-end Illumina FASTQ files:

1) **Trim Galore**: adapter trimming + low-quality filtering (includes FastQC)  
2) **Bowtie2**: alignment to **hg38**  
3) **featureCounts (Subread)**: gene-level quantification  
4) **limma-voom**: differential expression + QC + heatmaps/volcano plots  

> Raw FASTQs/BAMs are **not** included in the repository.

---

## Quick start

### 1) Edit paths/parameters
Open and edit:
- `config/config.yaml` (paths to FASTQs, hg38 Bowtie2 index, hg38 GTF, threads, strandedness)

### 2) Run mapping + counting
```bash
bash scripts/01_trim_map_count.sh
```

### 3) Run differential expression + plots
```bash
Rscript R/02_limma_voom_DE_heatmap.R --counts results/04_counts/featureCounts.txt --samples config/samples.tsv --outdir results/06_limma
```

---

## Notes
- The pipeline automatically **merges multi-lane FASTQs** per sample based on the prefix before the first underscore.
  Example: `1CI_S61_L001_R1_001.fastq.gz` → sample ID = `1CI`.
- Update contrasts in the R script if your paper uses comparisons other than `each group vs NT`.
- Update `strand` in `config/config.yaml` if libraries are stranded:
  - `0` unstranded, `1` stranded, `2` reverse-stranded.

---

## Folder structure
- `scripts/` : bash pipeline scripts
- `R/` : differential expression and plotting
- `config/` : sample sheet + config
- `envs/` : optional conda environment file
