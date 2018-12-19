#
# coreSNPs
#
# Author: Chunyu Zhao <zhaoc1@email.chop.edu>
# Created: 2018-10-24
#

import sys
import os
from pathlib import Path
import argparse
from snakemake.utils import listfiles

def read_samples(filepath):
  with open(filepath) as f:
    samples = f.read().splitlines()
    return(samples)

# Load config file
if not config:
  raise SystemExit("No config file specified. Run snakemake with specify with --configfile")

workdir: config['root']

rule list_samples:
  input:
    config['root']+'/'+config['data_fp']
  output:
    config['root']+'/'+config['samplelist_fp']
  params:
    config['filename_fmt']
  run:
    files = list(listfiles( str(Path(input[0])/Path("{sample}." + params[0])) ))
    samples = [t[1]['sample'] for t in files]
    with open(output[0], 'w') as f:
        f.write('\n'.join(samples) + '\n')


SAMPLES = read_samples( config['root'] + '/' + config['samplelist_fp'] )
CONDA_PREFIX = os.getenv("CONDA_PREFIX")

#### make sure checkM has downloaded its precalculated data files
rule checkm_dataset:
  output:
    ancient(config['root']+'/1_checkm/checkm_data')
  conda:
    "checkm.yml"
  params:
    conda_prefix = CONDA_PREFIX
  shell:
    """
    path_db="{params.conda_prefix}/opt/checkm_data"
    mkdir -p $path_db
    
    if [ ! -d $path_db/genome_tree ]; then
        wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz -P {params.conda_prefix}/opt
        tar -xzf {params.conda_prefix}/opt/checkm_data_2015_01_16.tar.gz -C $path_db
    fi
    
    ln -sf $path_db {output}
    """

rule run_checkm:
   input:
     contigs = config['root']+'/'+config['data_fp'],
     checkm_data = ancient(config['root']+'/1_checkm/checkm_data')
   output:
     results = config['root']+'/1_checkm/checkm_results/results.tsv',
     table = config['root']+'/1_checkm/checkm_results/checkm_results.tsv',
     lineage = config['root']+'/1_checkm/checkm_results/lineage.ms'
   params:
     outdir = config['root']+'/1_checkm/checkm_results',
     fmt = config['filename_fmt']
   conda: 
     "checkm.yml"
   threads: 8
   shell:
     """
     ## you can’t specify the path on the command line
     checkm data setRoot {input.checkm_data}
     
     ## analyze your genome using checkM
     checkm lineage_wf {input.contigs} {params.outdir} \
            -f {output.results} -x {params.fmt} -t {threads} --tab_table
     
     ## run checkM again for the better table
     checkm qa {output.lineage} {params.outdir} \
            -f {output.table} -t {threads} --tab_table -o 2
     """

rule _run_checkm:
   input:
     config['root']+'/1_checkm/checkm_results/checkm_results.tsv'


########### now it's time for pangenome analysis
rule run_prokka:
   input:
     config['root']+'/'+config['data_fp']+'/{sample}.fasta'
   output:
     config['root']+'/2_prokka/prokka_{sample}/{sample}.gff'
   params:
     outdir = config['root']+'/2_prokka/prokka_{sample}',
     sample = '{sample}'
   shell:
     """
     prokka --kingdom Bacteria --outdir {params.outdir} --prefix {params.sample} \
            --locustag {params.sample} --force {input}
     """

rule run_roary:
   input:
     expand(config['root']+'/2_prokka/prokka_{sample}/{sample}.gff', sample=SAMPLES)
   output: 
     config['root']+'/3_roary.done'
   params: 
     outdir = config['root']+'/3_roary_results/'
   threads:
     8
   shell:
     """
     roary -p {threads} -f {params.outdir} -z -e -n -v {input} && touch {output}
     """

