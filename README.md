# RNA-mediated interactions: An eRNAi target search algorithm for studying the impact of the metagenome

<img align=right src="https://clipart-library.com/images/BTaKAn6gc.jpg" alt="# Codon-optimization Tool" width="100"/>

- The phenomenon of environmental RNA interference (eRNAi) is based on the transfer of small RNA molecules (sRNA) between organisms to suppress the expression of target genes. Investigating the processes of environmental RNA interference provides novel insights into the dynamics of interactions between living beings and raises the question of whether RNA could be transferred between organisms.

## Aim

- The main goal of the project is the development of an algorithm to identify potential targets in an organism's genome that can be specifically regulated by environmental double-stranded RNAs.

## Contents
- [Pipeline](#Pipeline)
- [Dataset](#Dataset)
- [Installation](#Installation)
- [Usage](#Usage)
- [Testing](#Testing)
- [Results](#Results)
- [Contributing](#Contributing)
- [Literature](#Literature)
- [Citation](#Citation)

## Pipeline

1) Data analysis and selection of suitable target datasets
2) Data pre-processing 
3) Canonical Correspondence Analysis (CCA)
4) Identification of meta-genome and transcriptome overlap regions
5) Statistical validation of results  
4) Pipeline validation and testing 

<img src="https://github.com/Valeriisht/eRNAi_project/blob/dev/imgs/pipeline.png" />

## Dataset

The data were obtained from public sources.

- *Gallus gallus*, **BioProject PRJNA503784** (Tejas M. Shah et al. 2019)
- *Homo Sapiens*,  **NASA GeneLab - OSD-574** (Park J et al. 2024)

## Installation

To get the tool clone the git repository:

```sh
git clnone git@github.com:Valeriisht/eRNAi_project.git
```

Create a conda environment with the necessary packages. 
Activate it:

```sh
conda env create -f environment.yaml
conda activate eRNAi
```

Then you can use your data, such as metagenomic and transcriptomic data.
In order to use your own data, you need to specify the parameters (data references and SRA ID) in the configuration file (```config.yaml```).

## Usage

The scripts for primary data processing can be found in the following folder ```workflow/rules/```

To run the process, select the desired rule all in ```SnakeFile``` to process the data.

Example:

```
snakemake --cores 8
```

## Testing

Three types of tests: unit tests, integration tests, and generation tests cover the project.

To run it, execute the command:

```
pytest tests/ 
```


## Contributing 

[Contributing.md](docs/CONTRIBUTING.md).

## Citation
1) Park, J., Overbey, E. G., Narayanan, S. A., Kim, J., Tierney, B. T., Damle, N., Najjar, D., Ryon, K. A., Proszynski, J., Kleinman, A., Hirschberg, J. W., MacKay, M., Afshin, E. E., Granstein, R., Gurvitch, J., Hudson, B. M., Rininger, A., Mullane, S., Church, S. E., … Mason, C. E. (2024). Spatial multi-omics of human skin reveals KRAS and inflammatory responses to spaceflight. Nature Communications, 15(1), 4773. https://doi.org/10.1038/s41467-024-48625-2

2) Shah, T. M., Patel, J. G., Gohil, T. P., Blake, D. P., & Joshi, C. G. (2019). Host transcriptome and microbiome interaction modulates physiology of full-sibs broilers with divergent feed conversion ratio. Npj Biofilms and Microbiomes, 5(1), 24. https://doi.org/10.1038/s41522-019-0096-3

