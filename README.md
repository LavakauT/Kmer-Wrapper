# KmersDiscovery Wrapper
### A wrapper to identify k-mers and perform machine learning on the identified k-mers.
<img src="https://github.com/LavakauT/Kmer-Wrapper/assets/132649549/9bb29229-36a4-4698-8050-837487a63682" width="80%"> 

### Update:
> 19 Jun 2024: The version that supports multiple sets simultaneously has been released.
> 
> 18 Jun 2024: The current version supports only one set of TP genes and one set of TN genes. A version that supports multiple sets simultaneously will be released soon.

## Authors
- Wrapper developer: Mr. Pan Kai-Cheng (Dr. Cheng Chia-Yi)
- K-mer/Machine Learning source code: Shiulab (https://github.com/ShiuLab/MotifDiscovery; https://github.com/bmmoore43/ML-Pipeline)
- Optimize and develop applications: Mr. Lavakau Thalimaraw (https://github.com/LavakauT/KmersDiscovery_MEME; https://github.com/LavakauT/ML-pipeline)
- Modification concept: Dr. Liu Ming-Jung, Dr. Wu Ting-Ying, and Dr. Cheng Chia-Yi

## Simple guide for user
### FOR FIRST-TIME USE, you need to do some additional setup:
```
   1. change your direction to 'code_intergration' file (cd ~/code_intergration)
   2. source [your own conda.sh(~/miniconda3/etc/profile.d/conda.sh)]     **if you don't have one, please install miniconda**
   3. conda env create -n python_Kmer -f ./environment/python_Kmer.yml
   4. conda env create -n R_Kmer -f ./environment/R_Kmer.yml
   5. replace your ~/miniconda3/envs/R_Kmer/lib/R/library folder by ./enviroment/library
```

### FOR NORMAL USE, follow these steps:
```
   0. (put your .gff3 & genome file in ~/code_intergration/gffandgenome) AND (put your TP/TN data in ~/code_intergration/data/Kmer/TPTN)
   # You can follow README in gffandgenome folder to download Arabidopsis thaliana genome fasta and gff3 file for testing.
   # Name must start with "TP_(subgroup)_(number)" or "TN_(subgroup)_(number)" (mutiple files are avaliable). 
   # Which have same (subgroup)_(number) will be in the same subgroup; there's no interaction between different subgroups.


   1. change your direction to 'code_intergration' file (cd ~/code_intergration)
   2. source [your own conda.sh(~/miniconda3/etc/profile.d/conda.sh)]     **if you don't have one, please install miniconda**
   3. conda activate python_Kmer
   4. python ~/code_intergration/Intergrated_pipeline_BATCH.py -main_dir [YOUR MAIN FOLDER ABSOLUTE PATH] -conda_dir [YOUR conda.sh file ABSOLUTE PATH] -gff3_filename [YOUR .GFF3 FILE NAME] -features [FEATURES FROM GFF]-genome_filename [YOUR GENOME FILE NAME] -set_num [SETS OF TEST/TRAIN] -group_name [YOUR GROUP FOLDER NAME]
      Please only operate "Intergrated_pipeline_BATCH.py" but not "Intergrated_pipeline.py", even though you have only one set of TP/TN**
      (Please see the next part for more explanation of these arguments.)
```

*For "0_intergrated_codes.py", there are eight arguments you have to declare:*
   - *-main_dir*, Your main folder path(usually: ~/code_intergration)
   - *-conda_dir*, Your conda.sh path(usually: ~/miniconda3/etc/profile.d/conda.sh)
   - *-gff3_filename*, Your .gff3 file name(NOT PATH, only filename)
   - *-genome_filename*, Your genome file name(NOT PATH, only filename)
   - *-features*, Features from gff3 file that you want(For instance: exon, gene,mRNA). The format must be separated by "," and include NO SPACE!!!!
   -  *-set_num*, How many test/train sets you want to split (1~40)
   -  *-group_name*, Your new create group folder name
   - Deprecation, directly run in RF, SVM, and LogReg >>> *-ML_method*, ML Algorithm to run (RF, SVM, SVMpoly, SVMrbf, GB, LogReg)
   - Deprecation >> *-TP_filename*, Your TP data file name(NOT PATH)
   - Deprecation >> *-TN_filename*, Your TN data file name(NOT PATH)
   
   (The program will create a new folder with this name to store the output files, organizing them into a group. 
   This helps to ensure that outputs from different runs or different groups of data do not overlap.)

### HPC bash script
>Once you finish all the required settings, you can submit KMER.sh as a testing run with Sample Data (TP_Meristem_1 & TN_Meristem_1) and Sample Genome (TAIR10_Chr.all_cleaned.fasta & Araport11_201606.Evolinc.gff3) in your HPC system.

### Web Tool: PREDICT
>Try our web tool [PREDICT](https://predict.southerngenomics.org/kmers/kmers.php) for rapid k-mer identification and TFBM association.
>
><img src="https://github.com/LavakauT/Kmer-Wrapper/assets/132649549/b731ec02-d6b1-4853-8c29-811d4c190dff" width="60%">
