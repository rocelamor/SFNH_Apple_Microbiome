# SFNH Apple Microbiome — Analysis Code

Code accompanying the manuscript **"Bacterial diversity loss without community
convergence in a decomposing apple microbiome across a *Caenorhabditis
elegans* abundance gradient"** (Soil-Fruit Natural Habitat [SFNH] mesocosm
study).

Raw 16S rRNA gene and fungal ITS amplicon sequencing data are deposited in
the NCBI Sequence Read Archive under BioProject accession
[PRJNA1535416](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1535416).
This repository contains the processed/derived data tables and R scripts
used to generate every main-text and supplementary figure from those data.

## Repository structure

```
.
├── scripts/    R scripts, one per figure (see table below)
├── data/       Processed data tables read by the scripts
└── output/     Created automatically when a script is run
    ├── Figures/
    └── Tables/
```

## How to run

1. Clone or download this repository.
2. Open it as the working directory (e.g. open `scripts/` as an RStudio
   project at the repo root, or run from a terminal at the repo root).
3. Run any script, e.g.:
   ```r
   source("scripts/AM_Fig2.R")
   ```
   or from the command line:
   ```bash
   Rscript scripts/AM_Fig2.R
   ```
4. Each script reads its inputs from `./data` and writes figures/tables to
   `./output/Figures` and `./output/Tables`, created automatically.

Scripts are independent of each other and can be run in any order.

## Data files (`data/`)

| File | Used by |
|---|---|
| `metadata.csv` | AM_Fig2, AM_Fig3, AM_Fig4, AM_FigS1, AM_FigS2, AM_FigS3, AM_FigS5 |
| `bac_counts_rarefied14000.csv` | AM_Fig3, AM_Fig4, AM_FigS2, AM_FigS3, AM_FigS5 |
| `bac_abundance_class.csv` | AM_Fig3, AM_FigS3 |
| `its_metadata_apple.csv` | AM_FigS4 |
| `its_counts_rarefied14000_apple.csv` | AM_FigS4 |

## Scripts (`scripts/`) and the figures they produce

| Script | Figure | Contents |
|---|---|---|
| `AM_Fig2.R` | Fig. 2 | Local *C. elegans* abundance and developmental-stage composition in decomposing apple |
| `AM_Fig3.R` | Fig. 3 | Bacterial ASV abundance classification; partial-Spearman diversity/richness associations; negative-binomial GLMM dose-response; matched-control richness deficit |
| `AM_Fig4.R` | Fig. 4 | Bray-Curtis distance to matched control; between-cell dispersion; residual within-cell dispersion |
| `AM_FigS1.R` | Fig. S1 | Total nematode abundance per pot; whole-pot developmental-stage composition |
| `AM_FigS2.R` | Fig. S2 | Bacterial community structure: apple vs. soil; apple community vs. local nematode abundance |
| `AM_FigS3.R` | Fig. S3 | Sequential-removal sensitivity analysis of rare ASV richness |
| `AM_FigS4.R` | Fig. S4 | Fungal ITS non-response as a specificity control |
| `AM_FigS5.R` | Fig. S5 | Directional alignment of community displacement; high- vs. low-nematode-abundance convergence |

## Software requirements

- R (tested under R 4.4.3, macOS, aarch64-apple-darwin20)
- R packages: `tidyverse`, `ggplot2`, `patchwork`, `scales`, `vegan`,
  `lme4`, `DHARMa`, `FSA`, `MASS`, `emmeans`, `ggtext`, `grid`

Each script writes its own `sessionInfo()` output to `output/` on
completion, recording the exact package versions used for that run.

Install missing packages with:
```r
install.packages(c(
  "tidyverse", "patchwork", "scales", "vegan",
  "lme4", "DHARMa", "FSA", "MASS", "emmeans", "ggtext"
))
```

## License

MIT License — see `LICENSE`. <!-- confirm/replace if a different license is preferred -->

## Citation

If you use this code, please cite the associated manuscript
(citation to be added upon publication) and the sequence data BioProject
accession PRJNA1535416.
