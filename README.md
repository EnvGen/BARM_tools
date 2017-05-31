Trimmed reads filename should match finished_reads/[r,d]na/<sample_group>/*.fq.gz

BARM assembly should be placed in assembly/BARM.fa.gz

## Quantify a set of reads against the barm reference assembly/gene set

snakemake --dryrun quantify_all
