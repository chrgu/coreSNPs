# coreSNPs

In the [pitfalls paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0267-5), we domentrated the benefits of using positive controls in shotgun metagenomics sequencing for high-volume analysis. In particular, we used a marine Vibrio species (Vibrio campbellii) for all of our HiSeq runs at PCMP. Due to the high sequencing coverage and low diversity of the positive control samples, we are able to de novo assemble the Vibrio campbellii genome, and study the distributuion of the SNPs for the same genome being sequenced multiple times.

In this repo, we start with assembled and taxonomically annotated contigs from [sunbeam](https://github.com/sunbeam-labs/sunbeam) pipeline:

- extract V. campbellii contigs using [sbx_contigs](https://github.com/sunbeam-labs/sbx_contigs).
- assess the draft assemblies quatlify using [checkm](https://github.com/Ecogenomics/CheckM).
- pangenome analysis using [roary](https://sanger-pathogens.github.io/Roary/).

- extract core genes using in-house R script
- calculate SNPs for each core gene using snp-sites
- muscle (?) or at least subset muscle
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
  
