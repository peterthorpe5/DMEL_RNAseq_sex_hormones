

conda activate R

perl /mnt/shared/scratch/pthorpe/apps/trinityrnaseq-v2.14.0/util/abundance_estimates_to_matrix.pl --est_method salmon  --gene_trans_map none --name_sample_by_basedir --out_prefix DMEL.genes ./*/quant.genes.sf


perl /mnt/shared/scratch/pthorpe/apps/trinityrnaseq-v2.14.0/util/abundance_estimates_to_matrix.pl --est_method salmon  --gene_trans_map none --name_sample_by_basedir --out_prefix DMEL.all ./*/quant.sf