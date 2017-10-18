# BARM Tools
A set of tools that are useful to reuse the Baltic Sea Reference Metagenome.

## Requirements
These are the required external softwares in order to perform mapping, processing of bam files and counting mapped reads per gene. 

```
  bowtie2
  picard tools
  samtools
  htseq-count
```

## Installation
We recommend the usage of a conda virtual environment, see https://conda.io/miniconda.html 

```
  conda create -n barm_tools_env python=3
  source activate barm_tools_env
  conda install snakemake pandas
  git clone --recursive https://github.com/EnvGen/BARM_tools.git
  cd BARM_tools
```

Modify config.json (e.g. the the "picard_jars" variable which should correspond to the directory where your picard tools jars are located).

## Data
The data needs to be stored in a very precise manner since Snakemake recognizes files based on their exact path.

Trimmed reads should be put in:

    finished_reads/{r,d]}na/<sample_group>/*_R{1,2}.fq.gz

The BARM Assembly (available on ENA) as the file:

    assembly/BARM.fa.gz

The gene locations are fetched from the gff which should be located in the file:

    annotation/prodigal/all_annotated_sequences/BARM/proteins.gff2

Annotations in:

    annotation/EggNOG/summary_annotation/prodigal/all.EggNOG.standardized.tsv
   
BARM assembly should be placed in assembly/BARM.fa.gz

## Quantify a set of reads against the barm reference assembly/gene set

snakemake --dryrun quantify_all
