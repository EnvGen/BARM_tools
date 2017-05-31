__author__ = "Johannes Alneberg"
__license__ = "MIT"


import os
import glob

megahit_coassembly_contigs = "assembly/BARM.fa.gz"


config["bowtie2_quant_rules"]["samples"] = {}
config["bowtie2_quant_rules"]["sample_groups"] = sample_groups

# Add reads for mapping
for read_file in glob.glob("finished_reads/*na/*/*.fq.gz"):
    read_basename = os.path.basename(read_file)
    read_name = read_basename.replace(".fq.gz", "")
    sample_name = read_name.replace("_R1", "").replace("_R2", "")
    
    sample_type = read_file.split('/')[1]
    if sample_type not in ['dna','rna']:
        continue
    elif sample_type == 'rna':
        # RNA-samples have slightly different handling
        config["bowtie2_quant_rules"]["rna_samples"].append(sample_name)

    # Add sample to sample group
    sample_group = read_file.split('/')[2]
    if sample_group not in sample_groups:
        sample_groups[sample_group] = []

    if sample not in sample_groups[sample_group]:
        sample_groups[sample_group].append(sample)

    if sample_name in config["bowtie2_quant_rules"]["samples"]:
        config["bowtie2_quant_rules"]["units"][sample_name].append(read_file)
        config["bowtie2_quant_rules"]["units"][sample_name].sort()
    else:
        config["bowtie2_quant_rules"]["units"][sample_name] = [read_file]
        config["bowtie2_quant_rules"]["samples"][sample_name] = [sample_name]

config["bowtie2_quant_rules"]["references"]["megahit_coassembly_contigs"] = megahit_coassembly_contigs
config["bowtie2_quant_rules"]["reference_for_ref_set"]["megahit_coassembly"] = "megahit_coassembly_contigs"

config["bowtie2_rules"]["references"] = config["bowtie2_quant_rules"]["references"]
config["bowtie2_rules"]["units"] = config["bowtie2_quant_rules"]["units"]
config["bowtie2_rules"]["samples"] = config["bowtie2_quant_rules"]["samples"]
config["bowtie2_rules"]["mapping_params"] = config["bowtie2_quant_rules"]["mapping_params"]

config["taxonomic_annotation"]["sample_sets"] = config['bowtie2_quant_rules']["split_ref_sets"]

WORKFLOW_DIR = "snakemake-workflows/"

include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/mapping/bowtie2.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/annotation/prodigal.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/mapping/samtools.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/blast/rpsblast.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/quantification/bowtie2.rules")

# The index is large, need to use .bt2l indices
ruleorder: bowtie2_map_large > bowtie2_map

annotations = ['EggNOG', 'dbCAN', 'EC', 'pfam']
normalizations = ['raw_counts', 'tpm']


rule link_quantification:
    input: "quantification/bowtie2_genes/local/megahit_coassembly/annotations/{sample_group}/{annotation}.{normalization}.annotated.tsv.gz"
    output: "quantification/{sample_group}/{annotation}.{normalization}.annotated.tsv.gz"


rule quantify_all:
    input: expand("quantification/{sample_group}/{annotation}.{normalization}.annotated.tsv.gz", 
            sample_group = sample_groups,
            annotation = annotations,
            normalization = normalizations)
            
