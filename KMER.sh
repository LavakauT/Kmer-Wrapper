#!/bin/bash
#SBATCH -J KMER_v3
#SBATCH -o KMER_v3.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 03-00:00

module load Miniconda3/4.9.2
# or source your source file [your own conda.sh(~/miniconda3/etc/profile.d/conda.sh)] **if you don't have one, please inatall miniconda**

start_time=$(date +%s)
echo ">>>START RUNNING KMER<<<"

source activate python_Kmer
cd /RAID1/working/R425/lavakau/Kmer-Wrapper-main
python 0_intergrated_codes_modified.py -main_dir /RAID1/working/R425/lavakau/Kmer-Wrapper-main -gff3_filename Araport11_201606.Evolinc.gff3  -genome_filename TAIR10_Chr.all_cleaned.fasta  -TP_filename TPMeristem.txt -TN_filename TNMeristem.txt  -set_num 10 -group_name Meristem

echo ">>>FINISHED KMER RUN<<<"
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "RUNNING TIME:${elapsed_time}"

# Note!!
# If the program can't find kmer, please check whether feature names are the same between your database and gff_file
