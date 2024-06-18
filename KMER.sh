#!/bin/bash
#SBATCH -J KMER_v3
#SBATCH -o KMER_v3.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 03-00:00

module load Miniconda3/4.9.2

start_time=$(date +%s)
echo ">>>START RUNNING KMER<<<"

source activate python_Kmer
cd /RAID1/working/R425/lavakau/Kmer-Wrapper-main
python 0_intergrated_codes_modified.py -main_dir /RAID1/working/R425/lavakau/Kmer-Wrapper-main -gff3_filename Mus_musculus.GRCm39.112.chr.gff3  -genome_filename Mus_musculus.GRCm39.dna.toplevel.fa  -TP_filename TP.txt -TN_filename TN.txt  -set_num 10 -group_name Mus

echo ">>>FINISHED KMER RUN<<<"
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "RUNNING TIME:${elapsed_time}"
