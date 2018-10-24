# coreSNPs

For hybrid assembly using Illumina and Nanopore sequencing, we can easily get complete or near-complete bacterial genomes. When the sequencing depth is deep and/or the diverty is relatively low, we can even recover good quality draft genomes from shotgun metagenomics sequencing data (assessed by [checkm](https://github.com/Ecogenomics/CheckM)). 

In this repo, we showed the following things:
- assess draft genomes quality using checkm
- build pangenome analysis using roary
- extract core genes using in-house R script
- calculate SNPs for each core gene using snp-sites
- muscle (?) or at least subset muscle
- calcuate hamming distance
