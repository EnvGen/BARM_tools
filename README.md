# BARM Tools
A set of tools that are useful to reuse the Baltic Sea Reference Metagenome.

## Requirements

```
  bowtie2
  
```

## Installation

```
  conda create -n barm_tools_env python=3
  source activate barm_tools_env
  conda install snakemake
  git clone --recursive https://github.com/EnvGen/BARM_tools.git
  cd BARM_tools

```

## Data
Trimmed reads filename should match finished_reads/{r,d]}na/<sample_group>/*_R{1,2}.fq.gz

BARM assembly should be placed in assembly/BARM.fa.gz

## Quantify a set of reads against the barm reference assembly/gene set

snakemake --dryrun quantify_all
