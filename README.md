# coreSNPs

In the [pitfalls paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0267-5), we domentrated the benefits of using positive controls in shotgun metagenomics sequencing for high-volume analysis. In particular, we used a marine Vibrio species (Vibrio campbellii) for all of our HiSeq runs. Due to the high sequencing coverage and low diversity of the positive control samples, we are able to de novo assemble the Vibrio campbellii genome, and study the distributuion of the SNPs for the same genome being sequenced multiple times.

In this repo, we start with assembled contigs, and do the following analysis:

- extract V. campbellii contigs using [sbx_contigs](https://github.com/sunbeam-labs/sbx_contigs).
- assess the draft assemblies quatlify using [checkm](https://github.com/Ecogenomics/CheckM).
- pangenome analysis using [roary](https://sanger-pathogens.github.io/Roary/).

- extract core genes using in-house R script
- calculate SNPs for each core gene using snp-sites
- muscle (?) or at least subset muscle
- calcuate hamming distance
