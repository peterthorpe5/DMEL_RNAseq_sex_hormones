import os
handle = open("names_wanted.txt","w")
fastqc_script = open("fastqc_script.sh","w")

fastqc_script.write("fastqc ")


count = 0
for filename in os.listdir(".") :
    if not filename.endswith("1.fq.gz") : continue
    count = count + 1
    out = "trimmo_%d.sh" % count
    trimmo = open(out,"w")
    name = filename
    prefix = filename.split("_1.f")[0]
    exp = prefix
    
    handle.write("Sample_%s/\n/%s\n" % (prefix, name))
    fastqc_script.write("%s_paired.fq.gz " % (prefix))
    trimmo.write("#!/bin/bash -l\n#SBATCH -J trimmo_%d   #jobname\n#SBATCH --partition=medium \n#SBATCH --cpus-per-task=4\n#SBATCH --mem=10GB\n\n" % count)
    trimmo.write("cd /home/pthorpe/scratch/mike_r_tanya_rnaseq_sep2022/raw/\n")
    trimmo.write("conda activate shell_example\n\ntrimmomatic PE -threads 4 -phred33 %s_1.fq.gz %s_2.fq.gz ../trimmed/%s_R1_paired.fq.gz ../trimmed/%s_R1_unpaired.fq.gz ../trimmed/%s_R2_paired.fq.gz  ../trimmed/%s_R2_unpaired.fq.gz ILLUMINACLIP:/home/pthorpe/scratch/HS_inhertied_snps/raw_data/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:69\n" \
                 % (prefix, prefix, exp, exp, exp, exp))
    #print("java -jar /storage/home/users/pjt6/fly_selective_sweeps/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 %s_1.fq.gz %s_2.fq.gz ../trimmed/%s_R1_paired.fq.gz ../trimmed/%s_R1_unpaired.fq.gz ../trimmed/%s_R2_paired.fq.gz ../trimmed/%s_R2_unpaired.fq.gz ILLUMINACLIP:/storage/home/users/pjt6/fly_selective_sweeps/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:69\n" \
    #             % (prefix, prefix, exp, exp, exp, exp))
    trimmo.close()
    print("sbatch trimmo_%d.sh\n" % count)
fastqc_script.close()
handle.close()

