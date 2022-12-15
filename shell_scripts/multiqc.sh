#!/bin/bash -l
#SBATCH -J qc   #jobname
#SBATCH --partition=medium 
#SBATCH --cpus-per-task=2
#SBATCH --mem=4GB

cd /home/pthorpe/scratch/mike_r_tanya_rnaseq_sep2022/
conda activate multiqc

multiqc --interactive /home/pthorpe/scratch/mike_r_tanya_rnaseq_sep2022/



