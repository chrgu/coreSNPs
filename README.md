# coreSNPs

MAG: metagenome-assembled genomes

In this repo, we start with assembled and taxonomically annotated contigs from [sunbeam](https://github.com/sunbeam-labs/sunbeam) pipeline:

- assess the draft assemblies quatlify using [checkm](https://github.com/Ecogenomics/CheckM).
- pangenome analysis using [roary](https://sanger-pathogens.github.io/Roary/).

- extract core genes using in-house R script
- calculate SNPs for each core gene using snp-sites
- calcuate hamming distance

## Install `coresnps` Conda environment

  ```bash
  conda env update --name=coresnps --quiet --file env.yml
  ```

## Assess assemblies quality

Checkm requires python2.7 and we also take care of the precalcualted [checkm-database](https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz) in the `checkm_dataset` rule.

  ```bash
  snakemake --configfile config.yml _run_checkm --use-conda --cores 8
  ```
  
## Pangenome analysis

  ```bash
  snakemake --configfile config.yml run_roary --cores 8
  ```
  
