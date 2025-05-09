# Enviromental RNA-interference (eRNAi)

<img align=right src="https://clipart-library.com/images/BTaKAn6gc.jpg" alt="# Codon-optimization Tool" width="100"/>

## RNA-mediated interactions: An eRNAi target search algorithm for studying the impact of the metagenome

- The phenomenon of environmental RNA interference (eRNAi) is based on the transfer of small RNA molecules (sRNA) between organisms to suppress the expression of target genes. The design algorithm for finding eRNAi targets is based on multi-omics analysis by integrating metagenomic and transcriptomic data.

### Steps

1) Data analysis and selection of suitable target datasets
2) Data pre-processing 
3) Canonical Correspondence Analysis (CCA)
4) Identification of meta-genome and transcriptome overlap regions
5) Statistical validation of results  
4) Pipeline validation and testing 

<img src="https://github.com/Valeriisht/eRNAi_project/blob/dev/imgs/pipeline.png" />

## Contents
- [Installation](#Installation)
- [Testing](#Testing)
- [Contributing](#Contributing)

## Installation

To get the tool clone the git repository:

```sh
git clnone git@github.com:Valeriisht/eRNAi_project.git
```
Create a conda environment with the necessary packages. 
Activate it:

```sh
conda env create -f environment.yml
conda activate eRNAi
```

## Testing

Three types of tests: unit tests, integration tests, and generation tests cover the project.

To run it, execute the command:

```
pytest tests/ 
```

## Contributing 

[Contributing.md](docs/CONTRIBUTING.md).

