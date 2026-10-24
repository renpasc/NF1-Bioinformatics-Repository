# NF1 exon 31 (23a) splicing and its regulators in RNA-seq data

Scripts that estimate the fraction of *NF1* isoform II (exon 31 included) from RNA-seq junction reads, and measure expression of candidate exon 31 splicing regulators (*REST*, *CELF3–6* and *ELAVL2–4*). They provide independent benchmarks for ddPCR measurements of *NF1* isoform I and II fractional abundance.

The analysis was written for a Master of Diagnostic Genomics thesis (Queensland University of Technology).

## Background

*NF1* exon 31 (legacy exon 23a) is a 63 bp alternative exon at chr17:31,252,938–31,253,000 (GRCh38). It adds 21 amino acids to the GAP-related domain of neurofibromin. Transcripts that skip it encode isoform I, which has higher Ras-GAP activity. Transcripts that include it encode isoform II.

CELF and ELAVL RNA-binding proteins promote exon 31 skipping, and REST represses their expression. These scripts ask where and when exon 31 is included, and whether that tracks the expression of these regulators.

## Method

### Exon 31 PSI

All junction-based scripts use the same approach.

1. Collect exon-exon junction counts within *NF1*.
2. Identify exon 31 from its two inclusion junctions (exon 30–31 and exon 31–32) and its skipping junction (exon 30–32), and check it against the expected coordinates.
3. Calculate percent spliced in (PSI), the isoform II fraction, for each sample:
   `PSI = mean(left inclusion, right inclusion) / (mean(left inclusion, right inclusion) + skipping)`
4. Report per-sample estimates for samples with at least `MIN_READS` informative reads (10 by default, 20 in the exploratory scripts).
5. Report a pooled, read-weighted estimate per group with a Wilson 95% confidence interval.

### Regulator expression

Gene-level counts are converted to transcripts per million (TPM) using gene length, with TPM denominators calculated over all genes. Where a gene symbol maps to more than one gene ID, the ID with the most reads is kept. *NF1* is included as a reference gene.

## Scripts

Run the R scripts in RStudio with `source(file.choose(), echo = TRUE)`, or click **Source**. Outputs are written to the R working directory (`getwd()`).

### Public data: blood

| Script | Data | Output |
|---|---|---|
| `nf1_exon23a_gtex_snaptron.R` | GTEx whole blood and EBV-transformed lymphocytes, Snaptron junction counts (Wilks et al., 2018) | Exon 31 PSI per sample and pooled |
| `nf1_regulators_gtex_recount3.R` | GTEx whole blood and EBV-transformed lymphocytes, recount3 gene counts (Wilks et al., 2021) | Regulator TPM per sample and per tissue |
| `nf1_exon23a_pbmc_snaptron.R` | Snaptron SRA compilation, samples described as PBMCs (exploratory) | No healthy PBMC samples were found. Kept to document the search. |
| `nf1_exon23a_monaco_recount3.R` | Monaco et al. (2019) sorted immune cells, GEO GSE107011 (exploratory) | Coverage at exon 31 was too low to use |

### Public data: early neuronal induction

| Script | Data | Output |
|---|---|---|
| `nf1_busskamp_recount3.R` | Busskamp et al. (2014), GEO GSE60548. PGP1 iPSCs induced with NGN1/2, Days 0, 1, 3 and 4, three replicates per day. recount3 junction and gene counts. | Exon 31 PSI and regulator TPM per sample and per day |

The script looks up the SRA project from GEO, checks that it is in recount3 and installs any missing packages. On the first run it writes `busskamp_sample_days.csv` with the day parsed from each sample label. Check it, correct it if needed and re-run.

### Laboratory data: NGN2 neurons

These scripts run on nf-core/rnaseq output (Ewels et al., 2020). No sequencing data or sample names are included in this repository.

| Script | Input | Output |
|---|---|---|
| `nf1_exon31_psi_ngn2.R` | STAR junction tables, `star_salmon/log/*.SJ.out.tab` | Exon 31 PSI per sample, per line and day, and pooled per day, with a plot |
| `nf1_exon31_psi.py` | The same junction tables, or BAM files (Python alternative, needs `pysam` for BAM input) | Exon 31 PSI per sample and per day |
| `nf1_regulators_neurons_nfcore.R` | Salmon gene TPM, `star_salmon/salmon.merged.gene_tpm.tsv` | Regulator TPM per sample, line and day |

The line and day for each sample are parsed from file names of the form `LINE_D14` or `LINE_day14`. To set them explicitly, add a `sample_metadata.csv` with columns `sample`, `line` and `day` to the working directory. Libraries whose names contain `NTC` are skipped.

### Figures

| Script | Input | Output |
|---|---|---|
| `busskamp_regulators_figure.R` | `busskamp_nf1_samples.csv`, plus `nf1_regulators_gtex_summary.csv` and `nf1_regulators_neurons_summary.csv` if available | Heatmap of log₂(TPM + 1) across GTEx, Busskamp and NGN2 neurons, and a time course of each gene (600 dpi PNG) |
| `nf1_regulators_figure.py` | GTEx and NGN2 summary CSVs | Two-panel heatmap and time course (Python version, 600 dpi PNG) |

The figure script searches the working directory, then OneDrive and the home folder, for the optional summary files, and reports which it used. Blocks without data are left out.

## Results from public data

### Exon 31 PSI

| Dataset | Samples | Isoform II (PSI) |
|---|---|---|
| GTEx whole blood | 456 | 0.27 pooled (95% CI 0.25–0.30) |
| GTEx EBV-transformed lymphocytes | 139 | 0.41 pooled (95% CI 0.39–0.42) |
| Busskamp et al. (2014), Day 0 | 3 | 0.59 pooled (95% CI 0.47–0.71) |
| Busskamp et al. (2014), Day 1 | 3 | 0.51 pooled (95% CI 0.40–0.61) |
| Busskamp et al. (2014), Day 3 | 3 | 0.24 pooled (95% CI 0.15–0.36) |
| Busskamp et al. (2014), Day 4 | 3 | 0.19 pooled (95% CI 0.12–0.28) |

Isoform I was the predominant *NF1* transcript in both blood tissues. During NGN1/2 induction, exon 31 inclusion fell by about two-thirds, mainly between Day 1 and Day 3.

### Regulator expression

| Gene | EBV lymphocytes (median TPM) | Whole blood (median TPM) | Busskamp Day 0 (mean TPM) | Busskamp Day 4 (mean TPM) |
|---|---|---|---|---|
| *NF1* | 7.1 | 1.0 | 8.8 | 13 |
| *REST* | 11.3 | 1.8 | 33 | 1.5 |
| *CELF3* | <0.1 | <0.1 | 0.4 | 170 |
| *CELF4* | <0.1 | <0.1 | 0.7 | 15 |
| *CELF5* | <0.1 | <0.1 | 5.8 | 5.7 |
| *CELF6* | 1.1 | 0.9 | <0.1 | 5.3 |
| *ELAVL2* | <0.1 | <0.1 | 12 | 135 |
| *ELAVL3* | 0.1 | <0.1 | 1.5 | 561 |
| *ELAVL4* | <0.1 | <0.1 | 0.9 | 142 |

GTEx data: 196 EBV-transformed lymphocyte and 852 whole blood samples in recount3. *CELF3–5* and *ELAVL2–4* were effectively absent from blood. During NGN1/2 induction, *REST* fell and *CELF3*, *CELF4* and *ELAVL2–4* rose sharply, with *ELAVL3* already elevated by Day 1. Results from laboratory NGN2 neurons are reported in the thesis.

## Limitations

- Junction counts at exon 31 are low in many samples, especially whole blood and the Busskamp time course, so several estimates rely on pooling.
- Pooled confidence intervals treat reads as independent. They do not account for variation between donors or replicates and are likely to understate uncertainty.
- EBV-transformed lymphocytes are activated B cells and whole blood is dominated by granulocytes. Neither is equivalent to PBMCs.
- GTEx and Busskamp data are poly(A)-selected and were quantified with recount3 (GRCh38, GENCODE v26). The NGN2 data are total RNA libraries quantified with nf-core/rnaseq (STAR and Salmon, Ensembl 114). Compare patterns between datasets, not absolute TPM values.
- Gene-level TPM does not show which isoforms of each regulator are expressed.
- Expression changes in regulators and exon 31 inclusion are correlational. They do not establish a causal relationship.

## Requirements

R (4.3 or later) with these packages. The newer scripts install anything missing on first run.

```r
install.packages("BiocManager")
BiocManager::install(c("snapcount", "recount3", "SummarizedExperiment"))
install.packages(c("ggplot2", "ragg", "systemfonts"))
```

Python 3 with `matplotlib` and `numpy` for `nf1_regulators_figure.py`, and `pysam` only if running `nf1_exon31_psi.py` on BAM files.

The public-data scripts need internet access to NCBI GEO and the Snaptron and recount3 servers. Institutional networks may block these.

## Data privacy

This repository contains code only. Sequencing files, sample sheets and analysis outputs are excluded by `.gitignore`. Keep laboratory data out of commits.

## References

Busskamp, V., Lewis, N. E., Guye, P., Ng, A. H. M., Shipman, S. L., Byrne, S. M., Sanjana, N. E., Murn, J., Li, Y., Li, S., Stadler, M., Weiss, R., & Church, G. M. (2014). Rapid neurogenesis through transcriptional activation in human stem cells. *Molecular Systems Biology, 10*(11), 760. https://doi.org/10.15252/msb.20145508

Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines. *Nature Biotechnology, 38*(3), 276–278. https://doi.org/10.1038/s41587-020-0439-x

GTEx Consortium. (2017). Genetic effects on gene expression across human tissues. *Nature, 550*(7675), 204–213. https://doi.org/10.1038/nature24277

Monaco, G., Lee, B., Xu, W., Mustafah, S., Hwang, Y. Y., Carré, C., Burdin, N., Visan, L., Ceccarelli, M., Poidinger, M., Zippelius, A., Pedro de Magalhães, J., & Larbi, A. (2019). RNA-seq signatures normalized by mRNA abundance allow absolute deconvolution of human immune cell types. *Cell Reports, 26*(6), 1627–1640.e7. https://doi.org/10.1016/j.celrep.2019.01.041

Wilks, C., Gaddipati, P., Nellore, A., & Langmead, B. (2018). Snaptron: Querying splicing patterns across tens of thousands of RNA-seq samples. *Bioinformatics, 34*(1), 114–116. https://doi.org/10.1093/bioinformatics/btx547

Wilks, C., Zheng, S. C., Chen, F. Y., Charles, R., Solomon, B., Ling, J. P., Imada, E. L., Zhang, D., Joseph, L., Leek, J. T., Jaffe, A. E., Nellore, A., Collado-Torres, L., Hansen, K. D., & Langmead, B. (2021). recount3: Summaries and queries for large-scale RNA-seq expression and splicing. *Genome Biology, 22*(1), 323. https://doi.org/10.1186/s13059-021-02533-6

## Author

Renee Pascouau
