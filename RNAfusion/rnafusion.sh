#!/bin/bash

# conda activate nfcore

module load java/17.0.4.1
module load singularity/3.4.1
module load nextflow/21.10.6

# nextflow run /home/j_wang@intra.igr.fr/rnafusion/main.nf --allow-setuid --starindex --build_references --all --input samplesheet.csv --outdir ./ --genome GRCh38 -profile singularity -resume

nextflow run /home/j_wang@intra.igr.fr/rnafusion/main.nf --allow-setuid --starindex --input samplesheet.csv --outdir /mnt/beegfs/scratch/j_wang/Projects/01_BCC/20230110_RNAfusion --genome GRCh38 -profile singularity -resume
