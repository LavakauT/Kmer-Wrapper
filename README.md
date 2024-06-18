# KmersDiscovery Wrapper

A wrapper to identify k-mers and perform machine learning on the identified k-mers.
> Note: The current version supports only one set of TP genes and one set of TN genes. A version that supports multiple sets simultaneously will be released soon.

## Authors
- Wrapper developer: Mr. Pan Kai-Cheng (Dr. Cheng Chia-Yi)
- K-mer/Machine Learning source code: Shiulab (https://github.com/ShiuLab/MotifDiscovery; https://github.com/bmmoore43/ML-Pipeline)
- Optimize and develop applications: Mr. Lavakau Thalimaraw (https://github.com/LavakauT/KmersDiscovery_MEME; https://github.com/LavakauT/ML-pipeline)
- Modification concept: Dr. Liu Ming-Jung, Dr. Wu Ting-Ying, and Dr. Cheng Chia-Yi

## Simple guide for user
### FOR FIRST-TIME USE, you need to do some additional setup:
```
   1. change your direction to 'code_intergration' file (cd ~/code_intergration)
   2. source [your own conda.sh(~/miniconda3/etc/profile.d/conda.sh)]     **if you don't have one, please inatall miniconda**
   3. conda env create -n python_Kmer -f ./environment/python_Kmer.yml
   4. conda env create -n R_Kmer -f ./environment/R_Kmer.yml
   5. replace your ~/miniconda3/envs/R_Kmer/lib/R/library folder by ./enviroment/library
```

### FOR NORMAL USE, follow these steps:
```
   0. (put your .gff3 & genome file in ~/code_intergration/gffandgenome) AND (put your TP/TN data in ~/code_intergration/data/Kmer/TPTN)
   # You can follow README in gffandgenome folder to download Arabidopsis thaliana genome fasta and gff3 file for testing.

   1. change your direction to 'code_intergration' file (cd ~/code_intergration)
   2. source [your own conda.sh(~/miniconda3/etc/profile.d/conda.sh)]     **if you don't have one, please inatall miniconda**
   3. conda activate python_Kmer
   4. python ~/code_intergration/0_intergrated_codes_modified.py -main_dir [YOUR MAIN FOLDER ABSOLUTE PATH] -gff3_filename [YOUR .GFF3 FILE NAME] -genome_filename [YOUR GENOME FILE NAME] -TP_filename [YOUR TP DATA FILE NAME] -TN_filename [YOUR TN DATA FILE NAME] -set_num [SETS OF TEST/TRAIN] -ML_method [ML_ALG] -group_name [YOUR GROUP FOLDER NAME]
      (More explanation of these arguments, please see next part.)
```

*For "0_intergrated_codes_modified.py", there are eight arguments you have to declare:*
   1. *-main_dir*, Your main folder path(default: code_intergration)

   2. *-gff3_filename*, Your .gff3 file name(NOT PATH, only filename)
   3. *-genome_filename*, Your genome file name(NOT PATH, only filename)
   4. *-TP_filename*, Your TP data file name(NOT PATH)
   5. *-TN_filename*, Your TN data file name(NOT PATH)

   6. *-set_num*, How many test/train set you want to split (1~40)
   7. Deprecation, directly run in RF, SVM, and LogReg >>> *-ML_method*, ML Algorithm to run (RF, SVM, SVMpoly, SVMrbf, GB, LogReg)

   8. *-group_name*, Your new create group folder name
   (The program will create a new folder with this name to store the output files, organizing them into a group. 
   This helps to ensure that outputs from different runs or different groups of data do not overlap.)
