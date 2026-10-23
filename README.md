# NF1 exon 23a splicing in public RNA-seq data

R scripts that estimate the fraction of NF1 isoform II (exon 23a included) in public RNA-seq data. They provide an independent benchmark for ddPCR measurements of NF1 isoform I and II fractional abundance in peripheral blood mononuclear cells (PBMCs).

A second set of scripts measures expression of candidate exon 23a splicing regulators (REST, CELF3–6 and ELAVL2–4) in GTEx blood and in NGN2-induced neurons.

The analysis was written for a Master of Diagnostic Genomics thesis (Queensland University of Technology).

## Background

NF1 exon 23a is a 63 bp alternative exon. It adds 21 amino acids to the GAP-related domain of neurofibromin. Transcripts that skip it encode isoform I. Transcripts that include it encode isoform II.

Published estimates of the isoform balance in blood differ. Assunto et al. (2019) reported isoform II as the predominant form in blood. ddPCR in PBMCs for this thesis gave roughly 75% isoform I and 25% isoform II. These scripts use public RNA-seq junction counts as an independent check.

## Method

All three scripts follow the same approach.

1. Download all exon-exon junction counts within NF1.
2. Find every 63 bp cassette exon, defined by two inclusion junctions and one skipping junction.
3. Select the candidate whose weakest junction has the most reads. A true alternative exon uses all three junctions. Ranking by total reads instead can select artefacts that borrow one busy constitutive junction.
4. Check the selected exon against the expected exon 23a coordinates, chr17:31,252,938–31,253,000 (GRCh38).
5. Calculate the isoform II fraction for each sample as `mean(left inclusion, right inclusion) / (mean(left inclusion, right inclusion) + skipping)`.
6. Report per-sample estimates for samples with at least 20 informative reads. Also report a pooled, read-weighted estimate per group with a Wilson 95% confidence interval.

## Scripts

### `nf1_exon23a_gtex_snaptron.R` (main result)

Uses the Snaptron GTEx compilation. Compares whole blood with EBV-transformed lymphocytes. Tries the newer `gtexv2` compilation first and falls back to the original `gtex` compilation if snapcount fails to build the result.

Outputs `nf1_exon23a_gtex_samples.csv` and `nf1_exon23a_gtex_summary.csv`.

### `nf1_exon23a_pbmc_snaptron.R` (exploratory)

Searches the Snaptron SRA compilation for samples whose metadata mentions PBMCs. This did not find healthy PBMCs. The samples that passed the filter were K562 cells, in vitro polarised T helper cells, CD34+ cells undergoing erythroid differentiation and monocyte-derived macrophages. Their methods sections mentioned PBMC isolation. The script is kept to document this.

### `nf1_exon23a_monaco_recount3.R` (exploratory)

Uses recount3 junction counts for SRA project SRP125125. This project holds the Monaco et al. (2019) dataset of total PBMCs and 29 sorted immune cell types from healthy donors (GEO GSE107011). Coverage at exon 23a was too low to use. All 127 samples together gave only 79 informative reads.

## Results

| GTEx tissue | Samples | Samples with ≥20 reads | Isoform II fraction |
|---|---|---|---|
| Whole blood | 456 | 5 | 0.27 pooled (95% CI 0.25–0.30) |
| EBV-transformed lymphocytes | 139 | 87 | 0.40 median (IQR 0.33–0.46), 0.41 pooled (95% CI 0.39–0.42) |

Isoform I was the predominant NF1 transcript in both tissues. The whole blood estimate is close to the ddPCR result in PBMCs. Neither tissue showed isoform II predominance.

### Limitations (exon 23a)

- Most whole blood samples have low coverage at exon 23a, so the whole blood estimate relies on pooling.
- Pooled confidence intervals treat reads as independent. They do not account for variation between donors and are likely to understate uncertainty.
- EBV-transformed lymphocytes are activated B cells. They are not equivalent to resting PBMCs.
- Junction counts come from the Rail-RNA aligner used by Snaptron. They are not directly comparable with ddPCR.

## Splicing regulator expression

These scripts compare gene-level expression of candidate exon 23a regulators in blood and in neurons. NF1 is included as a reference gene.

### `nf1_regulators_gtex_recount3.R`

Uses recount3 gene counts for the GTEx BLOOD project. Keeps whole blood and EBV-transformed lymphocytes. Converts base-pair coverage to read counts with `transform_counts()`. Calculates TPM from all genes using gene length, then keeps the genes of interest. If a gene symbol maps to more than one gene ID, keeps the ID with the most reads.

Outputs `nf1_regulators_gtex_samples.csv`, `nf1_regulators_gtex_summary.csv` and, if ggplot2 is installed, `nf1_regulators_gtex.png`.

### `nf1_regulators_neurons_nfcore.R`

Reads `star_salmon/salmon.merged.gene_tpm.tsv` from nf-core/rnaseq output for NGN2 neuronal differentiation (days 7 to 35). Reads the cell line and day from `sample_metadata.csv` (columns `sample`, `line`, `day`) if present. Otherwise it parses the day from the sample name and takes the text before the day as the line.

Use the gene-level file. The transcript-level file (`salmon.merged.transcript_tpm.tsv`) has different columns and will be rejected.

Outputs `nf1_regulators_neurons_samples.csv`, `nf1_regulators_neurons_summary.csv` and, if ggplot2 is installed, `nf1_regulators_neurons.png`.

### `nf1_regulators_figure.py`

Builds the combined figure from the two summary CSVs. Panel A is a heatmap of log2(TPM + 1). It shows GTEx median TPM beside the mean TPM across neuronal lines for each day. Panel B shows the neuronal time course for each gene (mean ± SD per line).

```
pip install matplotlib numpy
python nf1_regulators_figure.py path/to/csv/folder
```

Writes `nf1_regulators_figure.png` at 600 dpi. Uses the Lato font if installed.

![Splicing regulator expression](nf1_regulators_figure.png)

### Regulator results

| Gene | EBV lymphocytes (median TPM) | Whole blood (median TPM) | NGN2 neurons (mean TPM, days 7–35) |
|---|---|---|---|
| NF1 | 7.1 | 1.0 | 67–89 |
| REST | 11.3 | 1.8 | 0.5–4.8 |
| CELF3–5 | <0.1 | <0.1 | 61–210 |
| CELF6 | 1.1 | 0.9 | 4–29 |
| ELAVL2–4 | ≤0.1 | <0.1 | 231–504 |

GTEx samples: 196 EBV-transformed lymphocytes and 852 whole blood. Neuronal data: two iPSC lines (PB001_1 and PB10_15), two replicates per time point.

CELF3–5 and ELAVL2–4 were effectively absent from blood. All six were highly expressed in neurons from day 7. CELF6 rose steadily across differentiation. REST fell after day 7.

### Limitations (regulators)

- GTEx data were quantified with recount3 (GENCODE v26). Neuronal data were quantified with nf-core/rnaseq (Salmon, Ensembl 114). Compare patterns between datasets, not absolute TPM values.
- Gene-level TPM does not show which isoforms of each regulator are expressed.
- Each neuronal time point has two replicates per line.

## Requirements

R with these Bioconductor packages:

```r
install.packages("BiocManager")
BiocManager::install(c("snapcount", "SummarizedExperiment"))   # Snaptron scripts
BiocManager::install("recount3")                                # Monaco and GTEx regulator scripts
install.packages("ggplot2")                                     # optional plots
```

Python 3 with matplotlib and numpy for `nf1_regulators_figure.py`.

The scripts need internet access to the Snaptron and recount3 servers. Institutional networks may block these.

## Usage

Run each script as a whole file so that it stops at the first error.

```r
source("nf1_exon23a_gtex_snaptron.R", echo = TRUE)
```

Output CSVs are written to the R working directory. The Monaco script downloads a large junction file and may need several GB of memory.

## References

Assunto, A., Ferrara, U., De Luca, A., Pivonello, C., Lombardo, L., Piscitelli, A., Tortora, C., Pinna, V., Daniele, P., Pivonello, R., Russo, M. G., Limongelli, G., Colao, A., Tartaglia, M., Strisciuglio, P., & Melis, D. (2019). Isoform-specific NF1 mRNA levels correlate with disease severity in Neurofibromatosis type 1. Orphanet Journal of Rare Diseases, 14(1), 261. https://doi.org/10.1186/s13023-019-1223-1

GTEx Consortium. (2017). Genetic effects on gene expression across human tissues. *Nature, 550*(7675), 204–213. https://doi.org/10.1038/nature24277

Monaco, G., Lee, B., Xu, W., Mustafah, S., Hwang, Y. Y., Carré, C., Burdin, N., Visan, L., Ceccarelli, M., Poidinger, M., Zippelius, A., Pedro De Magalhães, J., & Larbi, A. (2019). RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types. Cell Reports, 26(6), 1627-1640.e7. https://doi.org/10.1016/j.celrep.2019.01.041

Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines. *Nature Biotechnology, 38*(3), 276–278. https://doi.org/10.1038/s41587-020-0439-x

Wilks, C., Zheng, S. C., Chen, F. Y., Charles, R., Solomon, B., Ling, J. P., Imada, E. L., Zhang, D., Joseph, L., Leek, J. T., Jaffe, A. E., Nellore, A., Collado-Torres, L., Hansen, K. D., & Langmead, B. (2021). recount3: Summaries and queries for large-scale RNA-seq expression and splicing. *Genome Biology, 22*(1), 323. https://doi.org/10.1186/s13059-021-02533-6

Wilks, C., Gaddipati, P., Nellore, A., & Langmead, B. (2018). Snaptron: Querying splicing patterns across tens of thousands of RNA-seq samples. *Bioinformatics, 34*(1), 114–116. https://doi.org/10.1093/bioinformatics/btx547

## Author

Renee Pascouau
