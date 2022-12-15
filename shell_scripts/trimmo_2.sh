#!/bin/bash -l
#SBATCH -J trimmo_2   #jobname
#SBATCH --partition=medium 
#SBATCH --cpus-per-task=4
#SBATCH --mem=10GB

cd /home/pthorpe/scratch/mike_r_tanya_rnaseq_sep2022/raw/
conda activate shell_example

trimmomatic PE -threads 4 -phred33 CSHA1_1.fq.gz CSHA1_2.fq.gz ../trimmed/CSHA1_R1_paired.fq.gz ../trimmed/CSHA1_R1_unpaired.fq.gz ../trimmed/CSHA1_R2_paired.fq.gz  ../trimmed/CSHA1_R2_unpaired.fq.gz ILLUMINACLIP:/home/pthorpe/scratch/HS_inhertied_snps/raw_data/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:69
