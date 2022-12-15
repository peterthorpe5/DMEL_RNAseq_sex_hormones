import os
handle = open("names_wanted.txt","w")
fastqc_script = open("fastqc_script.sh","w")

fastqc_script.write("fastqc ")


count = 0
for filename in os.listdir(".") :
    if not filename.endswith("R1_paired.fq.gz") : continue
    count = count + 1
    out = "salmon_%d.sh" % count
    trimmo = open(out,"w")
    name = filename
    prefix = filename.split("_R")[0]
    exp = prefix
    
    handle.write("Sample_%s/\n/%s\n" % (prefix, name))
    fastqc_script.write("%s_paired.fq.gz " % (prefix))
    trimmo.write("#!/bin/bash -l\n#SBATCH -J salmon_%d   #jobname\n#SBATCH --partition=medium \n#SBATCH --cpus-per-task=4\n#SBATCH --mem=10GB\n\n" % count)
    trimmo.write("cd /home/pthorpe/scratch/mike_r_tanya_rnaseq_sep2022/trimmed\n")
    trimmo.write("conda activate salmon\n\n")
    
    trimmo.write("/home/pthorpe/scratch/mike_r_tanya_rnaseq_sep2022/genome/salmon-1.9.0_linux_x86_64/bin/salmon quant --libType A --geneMap /home/pthorpe/scratch/mike_r_tanya_rnaseq_sep2022/genome/dmel-all-r6.47.gtf --threads 4 -o %s -1 %s_R1_paired.fq.gz -2 %s_R2_paired.fq.gz --index /home/pthorpe/scratch/mike_r_tanya_rnaseq_sep2022/genome/DMEL_all_gene " % (exp, exp, exp) )
    #print("java -jar /storage/home/users/pjt6/fly_selective_sweeps/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 %s_1.fq.gz %s_2.fq.gz ../trimmed/%s_R1_paired.fq.gz ../trimmed/%s_R1_unpaired.fq.gz ../trimmed/%s_R2_paired.fq.gz ../trimmed/%s_R2_unpaired.fq.gz ILLUMINACLIP:/storage/home/users/pjt6/fly_selective_sweeps/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:69\n" \
    #             % (prefix, prefix, exp, exp, exp, exp))
    trimmo.close()
    print("sbatch salmon_%d.sh\n" % count)
fastqc_script.close()
handle.close()

