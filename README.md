# BARM Tools
A set of tools that are useful to reuse the Baltic Sea Reference Metagenome. This tutorial will show how to set up and perform mapping of your own fastq reads files against the Baltic Sea Reference Metagenome to achieve annotation counts for different types of annotations.

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

This will set up a virtual environment, install snakemake and pandas within the environment and further clone this repository in order to start the analysis. This tutorial assumes you're using this directory as the working directory.

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

Finally, the annotations per gene should go in:

    annotation/<annotation_type>/summary_annotation/prodigal/all.<annotation_type>.standardized.tsv

For example:
    
    annotation/EggNOG/summary_annotation/prodigal/all.EggNOG.standardized.tsv


## Quantify a set of reads against the barm reference assembly/gene set
Now you can let snakemake figure out what steps it needs to run in order to achieve the quantified annotations and genes.

If you have downloaded all annotations, you can ask snakemake to figure out how to create all annotation tables:

    snakemake --dryrun quantify_all

Otherwise you can ask for any of the files individually, e.g.:

    snakemake --dryrun summary_tables/<sample_group>/EggNOG.tpm.annotated.tsv.gz

Snakemake will run these steps for you if you remove the ```--dryrun```. Do check further options offered by snakemake by running

    snakemake --help


## Recommendations
It is recommended that you're samples are from the Baltic Sea or from an environment rather similar to this in order to achieve as high mapping rate as possible. A low mapping rate will otherwise most likely lead to low annotation rates.
