#!/usr/bin/env bash

#SBATCH -c 2
#SBATCH --mem=5G

module load Python
conda activate /scratch/blazqul/envs/DGE_analysis

cd /scratch/blazqul/alvaro/CGGA

outputDir="salmon_out"
metric=lengthScaledTPM
designFile="salmon_out/design.tab"
gtf='/scratch/blazqul/indexes/h_index/gencode.v41.primary_assembly.annotation.gtf'
outname="BIOD"

Rscript GetSalmonTPMs.R \
    $metric \
    $outputDir \
    $designFile \
    $gtf \
    $outname