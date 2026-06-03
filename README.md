# rnaseq-ism1-lung

Differential expression and gene set enrichment analysis of bulk RNA-seq data
from *Ism1* knockout (*Ism1-/-*) mouse lung tissue at postnatal day 7 (P7).

**MSc Bioinformatics — University of Edinburgh**  
Functional Genomic Technologies, 2025–2026

---

## Biological context

Isthmin 1 (ISM1) is a secreted protein with known anti-apoptotic and
anti-angiogenic functions. *Ism1-/-* mice develop spontaneous progressive lung
emphysema with enhanced immune cell infiltration, implicating ISM1 in pulmonary
inflammation. This analysis characterises the transcriptional changes in
*Ism1-/-* lung tissue at the onset of alveolarization (P7) using bulk RNA-seq
data from GEO study [GSE262122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262122).

---

## Repository contents

| File | Description |
|------|-------------|
| `RNAseq_analysis_fgt_ica.Rmd` | Full DESeq2 pipeline: QC, normalisation, PCA, DEG analysis, annotated results |
| `gsea_fgt_ica.R` | Pre-ranked GSEA using fgsea against MSigDB Hallmark mouse gene sets |

---

## Analysis overview

### 1. Data processing and quality control
- Input: featureCounts output (`mytable_feaures`) from STAR alignment against GRCm39
- Pre-filtering: removed 36,206 genes with total counts ≤ 1
- Normalisation: regularised log (rlog) transformation via DESeq2
- QC: PCA and sample-to-sample Euclidean distance heatmap

### 2. Differential expression (DESeq2)
- Design: `~ Genotype` (WT as reference level)
- Wald test; results ranked by absolute Wald statistic
- Significant DEGs (padj < 0.05): 126 out of 21,602 tested genes
- Top hits include *Cxcl5*, *Col2a1*, *Plac9*, *Tmem254* — consistent with
  published findings on inflammation and lung fibrosis

### 3. Gene Set Enrichment Analysis (fgsea)
- Pre-ranked gene list: sorted by Wald statistic (descending)
- Gene sets: MSigDB Hallmark mouse collection
- Top upregulated pathways: INTERFERON_GAMMA_RESPONSE (NES = 2.36),
  INFLAMMATORY_RESPONSE (NES = 2.35), ALLOGRAFT_REJECTION (NES = 2.28)
- Top downregulated pathway: HEME_METABOLISM
- Results are consistent with ISM1's known anti-inflammatory role

---

## Key results

| Gene | log2FC | Wald stat | padj |
|------|--------|-----------|------|
| Plac9 | -1.76 | -13.52 | 2.6e-37 |
| Cxcl5 | +4.35 | +11.99 | 2.9e-29 |
| Col2a1 | +2.17 | +11.59 | 2.6e-27 |
| Pirt | +1.93 | +10.33 | 1.9e-21 |
| Duxbl1 | -1.52 | -9.32 | 3.5e-17 |

---

## How to run

### Prerequisites

The analysis requires `mytable_feaures` (featureCounts output) and
`resultAnnot.RData` (cached Ensembl annotation). These are not included
in the repository as they contain raw count data from the original study.
Download the source data from GEO: [GSE262122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262122).

Install R dependencies:

```r
# Bioconductor packages
BiocManager::install(c("DESeq2", "Rsubread", "edgeR", "limma",
                       "affy", "QuasR", "pheatmap", "biomaRt", "fgsea"))

# CRAN packages
install.packages(c("ggplot2", "gplots", "stringr", "dplyr",
                   "RColorBrewer", "patchwork", "scatterplot3d"))
```

### Run the analysis

```r
# 1. Open RNAseq_analysis_fgt_ica.Rmd in RStudio and knit to PDF or HTML
# 2. After the .Rmd completes, run the GSEA script:
source("gsea_fgt_ica.R")
```

The `.Rmd` writes two CSV files used as input by the GSEA script:
- `result_treatment_sorted_not_ann.csv` — unannotated DEG results (used by GSEA)
- `result_treatment_sorted_ann.csv` — annotated DEG results

---

## Software and versions

| Tool | Version | Purpose |
|------|---------|---------|
| DESeq2 | 1.44 | Normalisation and differential expression |
| fgsea | 1.30 | Gene set enrichment analysis |
| biomaRt | 2.60 | Ensembl gene annotation |
| STAR | 2.7.3a | Read alignment (upstream, pre-processed) |
| featureCounts | 2.0.8 | Read counting (upstream, pre-processed) |
| R | 4.4 | Analysis environment |

---

## Data source

Shanmugasundaram M. et al. (2024). *Ism1 deficiency in mice exacerbates
bleomycin-induced pulmonary fibrosis with enhanced cellular senescence and
delayed fibrosis resolution.* hLife, 2(7), 342–359.
[doi:10.1016/j.hlife.2024.05.006](https://doi.org/10.1016/j.hlife.2024.05.006)

GEO accession: [GSE262122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262122)

---

## Notes on code attribution

Both scripts are adapted from class templates provided by Dr Simon Tomlinson
(University of Edinburgh, Functional Genomic Technologies, 2026). Sample
annotation, analysis design, statistical interpretation, and GSEA configuration
were authored independently as part of the coursework.

---

## Author

Mario Antonio Rodriguez Diaz  
MSc Bioinformatics, University of Edinburgh  
Chevening Scholar 2025–2026  
[linkedin.com/in/mario-rd](https://linkedin.com/in/mario-rd) · [ORCID: 0009-0008-0104-7421](https://orcid.org/0009-0008-0104-7421)
