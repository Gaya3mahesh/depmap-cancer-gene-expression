# VEGFA and PCNA expression in breast and lung cancer cell lines (DepMap)

A small exploratory project in R using public DepMap data. I work on nanoparticle delivery of doxorubicin and withaferin A in breast and lung cancer cells, and in that work I measure VEGFA (angiogenesis) and PCNA (proliferation). This project asks a simple question: where do the cell lines I actually use sit among other breast and lung cancer cell lines in public RNA-seq data?

I'm a molecular and cell biologist learning bioinformatics, so this is a learning project. It is descriptive, and I have not run statistical tests yet.

## Question

Compared with other breast and lung cancer cell lines, are MDA-MB-231, MCF-7 and NCI-H358 high or low for VEGFA and PCNA at baseline?

## Data

All data is public and comes from the DepMap portal (https://depmap.org/portal), release DepMap Public 25Q2. The data files are not included in this repository. Please download them yourself from the DepMap Data Downloads page and place them in a folder called `data`.

Files used:
- `OmicsExpressionProteinCodingGenesTPMLogp1Stranded.csv` (gene expression, log2 TPM + 1, one row per cell line)
- `Model.csv` (cell line names, lineage and disease annotation)
- `OmicsSomaticMutationsMatrixHotspot.csv` (hotspot mutation calls, used for KRAS status)

Citation for the data: DepMap, Broad (2025). DepMap Public 25Q2. Dataset. depmap.org

## Methods

1. Pulled VEGFA, PCNA, ESR1, PGR and ERBB2 expression for all cell lines and joined it to the cell line annotation.
2. Kept breast and lung lineage lines, and highlighted MDA-MB-231, MCF-7 and NCI-H358.
3. For lung, restricted to non-small cell lung cancer and split lines by KRAS hotspot mutation status.
4. DepMap does not label triple-negative breast cancer, so breast lines were grouped by expression into "TNBC-like" (ESR1 < 2 and ERBB2 < 7), "ER-positive-like" (ESR1 >= 2) and "HER2-high, ER-low" (the rest). These cutoffs are my own rough choices, not a validated classifier. I only checked them against two lines I know: MDA-MB-231 came out TNBC-like and MCF-7 came out ER-positive-like.

## Results

Baseline expression of my cell lines (log2 TPM + 1):

| Cell line | Group | VEGFA | PCNA |
|---|---|---|---|
| MDA-MB-231 | TNBC-like breast | 7.57 | 7.97 |
| MCF-7 | ER-positive-like breast | 5.49 | 7.50 |
| NCI-H358 | KRAS-mutant NSCLC | 8.02 | 8.73 |

What the plots show:
- VEGFA: MDA-MB-231 and NCI-H358 sit towards the high end of their groups. NCI-H358 is well above the bulk of the NSCLC lines. MCF-7 sits near the middle of the ER-positive-like lines.
- PCNA: all three lines sit close to the middle of their groups, so PCNA does not separate them.
- KRAS hotspot mutant and non-mutant NSCLC lines look similar for both genes.
- TNBC-like lines look somewhat higher for VEGFA than ER-positive-like lines, while PCNA looks similar across breast groups. This is a visual impression only.

Figures are in the `figures` folder:
- `expression_breast_lung.png`
- `expression_breast_subtype.png`
- `expression_nsclc_kras.png`

## Limitations

- This is baseline expression in untreated cell lines. It says nothing about how treatment changes these genes.
- Cell lines are not tumours, and the groups are small (tens of lines).
- The TNBC-like grouping is a proxy based on my own expression cutoffs.
- No statistical testing has been done, so nothing here should be read as a significant difference.
- Doxorubicin sensitivity is not included yet.

## How to reproduce

1. Install R and the tidyverse: `install.packages("tidyverse")`
2. Download the three files listed above into a `data` folder.
3. Open `analysis.R` in RStudio, set the working directory to this project folder, and run it.

The script writes the figures to `figures` and a small results table to `results_my_lines.csv`.

## Next steps

- Add doxorubicin sensitivity data and check whether VEGFA or PCNA expression relates to sensitivity across breast and lung lines.
- Add tumour data from TCGA to compare cell lines with real tumours.
- Add statistical tests and a better breast subtype classification.

## About

Built as a learning project, with AI assistance for writing and explaining the R code. I checked the outputs against cell lines whose biology I know.
