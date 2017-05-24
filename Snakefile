__author__ = "Johannes Alneberg"
__license__ = "MIT"


import os
import sys
import shutil
import glob
from subprocess import check_output

megahit_coassembly_contigs = "assembly/megahit_coassembly/default/final.contigs.fa"


config["bowtie2_quant_rules"]["samples"] = {}

for sample_t in config["megahit_rules"]["samples"].items()):
    sample, units = sample_t
    config["bowtie2_quant_rules"]["units"][sample] = units
    config["bowtie2_quant_rules"]["samples"][sample] = [sample]

config["bowtie2_quant_rules"]["sample_groups"] = sample_groups

# Add reads for mapping that was not part of assembly
config["bowtie2_quant_rules"]["rna_samples"] = []
for read_file in glob.glob("only_for_mapping/finished_reads/rna/*.fq.gz"):
    read_basename = os.path.basename(read_file)
    read_name = read_basename.replace(".fq.gz", "")
    sample_name = read_name.replace("_R1", "").replace("_R2", "")
    
    # RNA-samples have slightly different handling
    config["bowtie2_quant_rules"]["rna_samples"].append(sample_name)
    if sample_name in config["bowtie2_quant_rules"]["samples"]:
        config["bowtie2_quant_rules"]["units"][sample_name].append(read_file)
        config["bowtie2_quant_rules"]["units"][sample_name].sort()
    else:
        config["bowtie2_quant_rules"]["units"][sample_name] = [read_file]
        config["bowtie2_quant_rules"]["samples"][sample_name] = [sample_name]

    config["fastqc_rules"]["reads"]["finished_reads_only_for_mapping_rna_" + read_name] = read_file


for read_file in glob.glob("only_for_mapping/finished_reads/dna/*.fq.gz"):
    read_basename = os.path.basename(read_file)
    read_name = read_basename.replace(".fq.gz", "")
    sample_name = read_name.replace("_R1", "").replace("_R2", "")
    
    if sample_name in config["bowtie2_quant_rules"]["samples"]:
        config["bowtie2_quant_rules"]["units"][sample_name].append(read_file)
        config["bowtie2_quant_rules"]["units"][sample_name].sort()
    else:
        config["bowtie2_quant_rules"]["units"][sample_name] = [read_file]
        config["bowtie2_quant_rules"]["samples"][sample_name] = [sample_name]

config["bowtie2_quant_rules"]["references"] = {"megahit_coassembly_genes": megahit_coassembly_genes}
config["bowtie2_quant_rules"]["references"]["megahit_coassembly_contigs"] = megahit_coassembly_contigs
config["bowtie2_quant_rules"]["reference_for_ref_set"]["megahit_coassembly"] = "megahit_coassembly_contigs"

config["bowtie2_rules"]["references"] = config["bowtie2_quant_rules"]["references"]
config["bowtie2_rules"]["units"] = config["bowtie2_quant_rules"]["units"]
config["bowtie2_rules"]["samples"] = config["bowtie2_quant_rules"]["samples"]
config["bowtie2_rules"]["mapping_params"] = config["bowtie2_quant_rules"]["mapping_params"]

# Add reads that will go through internal standards protocol before added to the pipeline
for read_file in glob.glob("only_for_mapping/with_internal_standards/*/*.fq.gz"):
    read_basename = os.path.basename(read_file)
    read_name = read_basename.replace(".fq.gz", "")
    sample_name = read_name.replace("_R1", "").replace("_R2", "")
    
    # Adding to bowtie2_rules instead of bowtie2_quant_rules
    # So that the mapping against the reference genome can 
    # work without adding the samples to the pipeline 
    if sample_name in config["bowtie2_quant_rules"]["samples"]:
        config["internal_standards"]["samples"].append(sample_name)
        config["bowtie2_rules"]["units"][sample_name].append(read_file)
        config["bowtie2_rules"]["units"][sample_name].sort()
    else:
        config["bowtie2_rules"]["units"][sample_name] = [read_file]
        config["bowtie2_rules"]["samples"][sample_name] = [sample_name]

# Load information on which internal standard have been used
with open("internal_standards.json") as i_s:
    internal_standards = json.load(i_s)

config["internal_standards"]["sample_to_reference"] = {}

for sample, reference in internal_standards.items():
    ref_name = reference.split('/')[-2]
    if ref_name not in config["bowtie2_rules"]["references"]:
        config["bowtie2_rules"]["references"][ref_name] = reference

    config["internal_standards"]["sample_to_reference"][sample] = ref_name

config["taxonomic_annotation"]["sample_sets"] = config['bowtie2_quant_rules']["split_ref_sets"]


WORKFLOW_DIR = "../../snakemake-workflows/"

include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/mapping/bowtie2.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/annotation/prodigal.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/annotation/eggnog.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/annotation/ec.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/annotation/dbcan.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/annotation/pfam.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/mapping/samtools.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/blast/rpsblast.rules")
#include: os.path.join(WORKFLOW_DIR, "rules/quantification/rpkm.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/trimming/cutadapt.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/quality_control/fastqc.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/duplicate_removal/fastuniq.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/assembly/megahit.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/annotation/prokka.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/quantification/bowtie2.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/taxonomic_annotation/megan.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/quality_control/internal_standards.rules")


# The index is large, need to use .bt2l indices
ruleorder: bowtie2_map_large > bowtie2_map

## Enable this when internal_standards mapping is done
#ruleorder: bowtie2_map > bowtie2_map_large

rule preprocess_all:
    input:
        htmls=expand("fastqc/{reads}/{reads}_fastqc.html", reads=config["fastqc_rules"]["reads"]),
        zips=expand("fastqc/{reads}/{reads}_fastqc.zip", reads=config["fastqc_rules"]["reads"])

rule quantify_all:
    input:
        expand("quantification/{assembly}/orf/{samples}/{samples}.rpkm",
            samples = ["P1414_101", "P1414_102", "P1414_103"],
            assembly = "assembly_v1"),


# Testing rules and configs
rule cutadapt_all_test:
    """Trim all reads with all supplied trimming parameters"""
    input:
        trimmed_reads=expand("cutadapt/adapt_cutting/{trim_params}/{reads}_{ext}.fq.gz",
        reads={"120628": ['samples/raw/120628_R1.fq.gz', 'samples/raw/120628_R2.fq.gz']},
        trim_params=config["cutadapt_rules"]["trim_params"],
        ext=["R1","R2"])

test_reads_orig = {
            "P1994_101_R1": "samples/raw/P1994_101_R1.fq.gz",
            "P1994_101_R2": "samples/raw/P1994_101_R2.fq.gz",
            "P1994_102_R1": "samples/raw/P1994_102_R1.fq.gz",
            "P1994_102_R2": "samples/raw/P1994_102_R2.fq.gz",
            "P1994_110_R1": "samples/raw/P1994_110_R1.fq.gz",
            "P1994_110_R2": "samples/raw/P1994_110_R2.fq.gz"
        }

test_reads = {}
for read_name, read_file in test_reads_orig.items():
    test_reads[read_name] = read_file
    read_basename = os.path.basename(read_file)
    for trim_params_name, trim_params_dict in config["cutadapt_rules"]["trim_params"].items():
        test_reads["cutadapt_" + trim_params_name + "_" + read_name] = \ 
            "cutadapt/adapt_cutting/{trim_params}/{read}".format(
                trim_params=trim_params_name,
                read = read_basename
            )

        test_reads["fastuniq_"+trim_params_name+"_"+read_name] = \
            "fastuniq/{trim_params}/{read}".format(
                trim_params=trim_params_name,
                read = read_basename
            ) 

rule fastqc_all_test:
    input:
        htmls=expand("fastqc/{reads}/{reads}_fastqc.html", reads=test_reads),
        zips=expand("fastqc/{reads}/{reads}_fastqc.zip", reads=test_reads)

