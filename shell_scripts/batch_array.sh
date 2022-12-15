#!/bin/bash -l
#SBATCH -J tri_Batch  #jobname
#SBATCH --partition=medium 
#SBATCH --cpus-per-task=4
#SBATCH --mem=10GB
#SBATCH --array=1-54

cd /home/pthorpe/scratch/mike_r_tanya_rnaseq_sep2022/raw/
conda activate shell_example

./trimmo_$SLURM_ARRAY_TASK_ID.sh
