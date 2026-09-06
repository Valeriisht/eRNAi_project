# Repository files

## Repository structure

| File/folder | Purpose | Status | Separate documentation |
| --- | --- | --- | --- |
| `README.md` | Brief project description and results obtained | up to date | `README.md` |
| `Snakefile` | Snakemake: defines the rules | up to date | `SNAKEMAKE.md` |
| `config/config.yaml` | Set of paths and parameters | needs updating | `SNAKEMAKE.md` |
| `enviromental.yaml` | Conda environment | incomplete | `--` |
| `workflow/rules/` | One Snakemake rule per task | needs further work | `SNAKEMAKE.md` |
| `CCA_analysis/` | R analysis and CCA/statistics plots | needs updating | `R_ANALYSIS.md` |
| `scripts/` |  Python and Bash scripts | partially implemented | `SCRIPTS.md` |
| `results/kraken_res/` | Saved Kraken2/Bracken results | needs re-checking | -- |
| `img/` | Illustrations and plots | -- | `README.md` |
| `tests/` | Code checks | not implemented | `SCRIPTS.md` |

## Service files

| File | What it represents | Note |
| --- | --- | --- |
| `.github/workflows/conda.yml` | GitHub Actions CI | on push, creates the conda environment, runs flake8 and pytest |
| `.gitignore` | list of files Git should not track | should be reviewed and correct |
| `requirements/installations.sh` | list of conda install commands | covered in detail in `SCRIPTS.md` |
| `requirements/minimal.txt` | minimal Python development list | not sufficient for the full pipeline |
| `requirements/tests.txt` | for tests | needs further work |

### Critical pairing rule

The file that expresses the correspondence between the two omics layers is `CCA_analysis/mapping_table.tsv`.

| Field | Meaning |
| --- | --- |
| `meta` | SRA ID of the metagenomic sample |
| `trans` | SRA ID of the host RNA-seq sample |

Before any CCA, the rows of the two matrices must be ordered strictly according to this correspondence. One matrix row = one paired biological unit.

## Classification results

`results/kraken_res/k2/` contains 24 files of the form `SRR*.nt.report`.
These are hierarchical Kraken2 reports: percentage of reads, number of reads, taxonomic rank, NCBI taxid, and taxon name.

`results/kraken_res/bracken/` contains Bracken results for different levels:

| Suffix | Interpretation |
| --- | --- |
| `.S.bracken` | species-level abundance table |
| `.G.bracken` | genus-level abundance table |
| `.P.bracken` | phylum-level abundance table |

The R script `Preprocessing_CCA.R` expects the report files and converts them into an abundance matrix suitable for CCA.

## Plots

| Path | What it shows | Source script |
| --- | --- | --- |
| `img/CCA_graphics.png` | two mutually inverse CCA ordinations | `CCA_analysis/CCA.R` |
| `img/statistic_tests/fisher_cca.png` | 3x2 Fisher comparison of CCA groups | `Statistic_Analysis.R` |
| `img/statistic_tests/binomial_cca.png` | binomial enrichment across groups | `Statistic_Analysis.R` |
| `img/statistic_tests/CCA_out(graphic1).png` | CCA vs. non-CCA Fisher comparison | `Statistic_Analysis.R` |
| `img/statistic_tests/windrose.png` | comparison of distributions across organisms | `Statistic_Analysis.R` |
| `CCA_analysis/cca_image/*.png` | one directionality plot per microorganism, for its transcripts | `Statistic_Analysis.R` |

