import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
	'-main_dir',
	help='Direction of the main folder (code_intergration)',
  required=True)
parser.add_argument(
	'-gff3_filename',
	help='Your .gff3 file name(NOT PATH)',
  required=True)
parser.add_argument(
	'-genome_filename',
	help='Your genome file name(NOT PATH)',
  required=True)
parser.add_argument(
	'-TP_filename',
	help='Your TP data file name(NOT PATH)',
  required=True)
parser.add_argument(
	'-TN_filename',
	help='Your TN data file name(NOT PATH)',
  required=True)

parser.add_argument(
	'-set_num',
	help='How many test/train set you want to split --> 1~40',
  default=1)
parser.add_argument(
	'-ML_method',
	help='ML Algorithm to run --> RF, SVM, SVMpoly, SVMrbf, GB, LogReg',
  required=False)
parser.add_argument(
	'-group_name',
	help='Customize your output folder name',
  required=True)

args = parser.parse_args()

args.TP_filename = args.TP_filename.split(".")[0]
args.TN_filename = args.TN_filename.split(".")[0]

dir_this_funtion_at = args.main_dir
import sys
sys.path.append(dir_this_funtion_at + "/package/KmersDiscovery")
sys.path.append(dir_this_funtion_at + "/package")
sys.path.append(dir_this_funtion_at + "/package/ML_pipeline")
import os
import subprocess


os.chdir(dir_this_funtion_at)
print("Current working directory:" + os.getcwd())
new_group_folder = os.getcwd() + "/" + str(args.group_name)


def main(gff3_filename, genome_filename,TP_filename, TN_filename, set_num, ML_method, new_group_folder):
    ##  create ./Data/ML/Train/...
    command = f"python ./package/create_TN_filename_folder.py -TN_filename {args.TN_filename}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)

    ##  (all) gff3 ---> 5utr.coord
    print("--------------------.coord file generating--------------------")
    gff3_path = "./gffandgenome/" + gff3_filename
    command = f"python ./package/FastaManager_modified_Ray.py -f gff_prom_to_coord_5utr -gff {gff3_path}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    # FM.fasta_manager().gff_prom_to_coord_5utr(gff=gff3_path)
    print("----------------.coord file generation COMPLETE----------------\n")
    
    ##  (all) 5utr.coord ---> 5utr.coord.fa
    print("-------------------.coord.fa file generating-------------------")
    coords_path = gff3_path + "_prom-5utr.coord"
    genome_path = "./gffandgenome/" + genome_filename
    command = f"python ./package/FastaManager_modified_Ray.py -f get_stretch4 -coords {coords_path} -fasta {genome_path}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    # FM.fasta_manager().get_stretch4(fasta=genome_path, coords=coords_path, seqid=0)
    print("---------------.coord.fa file generation COMPLETE---------------\n")

    ##  (target) split test_set & train_set
    print("-----------------test_set & train_set splitting-----------------")
    command = f"python ./package/testset_split.py -set_num {set_num} -TP_filename {args.TP_filename} -TN_filename {args.TN_filename}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    # tss.split_test_train_set()
    print("-------------test_set & train_set splitting COMPLETE-------------\n")

    ##  (target) test/train_set.txt ---> txt.fasta
    print("-----------------targer genes getting seqence-------------------")
    for folder_path in ["./Data/Kmer/Test", "./Data/Kmer/Train", f"./Data/Kmer/Train/{args.TN_filename}"]:
        all_files = os.listdir(folder_path)
        txt_files = [file for file in all_files if file.endswith(".txt")]
        for file_name in txt_files:
            file_path = os.path.join(folder_path, file_name)
            command = f"python ./package/FastaManager_modified_Ray.py -f getseq2 -fasta {coords_path}.fa -name {file_path}"
            subprocess.run("bash -c '" + command + "'", shell=True, check=True)
            # FM.fasta_manager().getseq2(fasta=coords_path+".fa", name=file_path)
    print("----------------targer genes get seqence COMPLETE----------------\n")

    ##  (target) making dataframe
    print("---------------------------making dataframe-----------------------")
    command = "for inp in ./Data/Kmer/Train/" + args.TN_filename + "/*.txt.fa; do neg=$inp; pos=./Data/Kmer/Train/${inp##*/}; pos=${pos/" + args.TN_filename + "/" + args.TP_filename + "}; python ./package/pCRE_Finding_FET.py -pos $pos -neg $neg -k ./package/6mer.txt FDR Y -save \"$inp.pcre\"; done"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("-----------------------making dataframe COMPLETE------------------\n")

    ##  (target) pos_neg_convertion
    print("--------------------------pos_neg_convertion----------------------")
    command = f"source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh && conda deactivate && conda activate R_Kmer && Rscript ./package/pos_neg_conversion.R -TN_filename {args.TN_filename}"
    subprocess.run("bash -c '" + command + "'", shell=True)
    print("----------------------pos_neg_convertion COMPLETE------------------\n")

    ##  (target) pcc filtering
    print("------------------------------pcc_filtering--------------------------")
    command = f"source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh && conda deactivate && conda activate R_Kmer && Rscript ./package/pcc_filtering.R -set_num {set_num}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("--------------------------pcc_filtering COMPLETE----------------------\n")

    ##  (target) ML_preprocess
    print("-----------------------------ML_preprocessing-------------------------")
    command = "for inp in ./Data/ML/Train/" + args.TN_filename +"/*.txt.fa.pcre_df_p0.01.txt; do python ./package/ML_preprocess.py -df $inp; done"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("-------------------------ML_preprocessing COMPLETE---------------------\n")

    ##  (target) generate ML test set
    print("---------------------------ML test_set generating-----------------------")
    command = f"python ./package/test_Kmercount.py -set_num {set_num} -TP_filename {args.TP_filename} -TN_filename {args.TN_filename}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("-----------------------ML test_set generation COMPLETE-------------------\n")

    ## (target)  pos_neg_convertion (ML_test)
    print("-------------------------pos_neg_convertion(ML_test)---------------------")
    command = "source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh && conda deactivate && conda activate R_Kmer && Rscript ./package/pos_neg_conversion_Testset.R"
    subprocess.run("bash -c '" + command + "'", shell=True)
    print("---------------------pos_neg_convertion(ML_test) COMPLETE-----------------\n")

    ##  (target) pre_RF
    print("--------------------------------doing pre_RF----------------------------")
    command = f"python ./package/pre_RF.py -set_num {set_num} -TN_filename {args.TN_filename}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("-------------------------------pre_RF COMPLETE---------------------------\n")
    
    ##  (target) ML
    print("-----------------------------------doing ML-----------------------------")
    command = f"for A in SVM LogReg RF; do for G in full random; do for inp in ./Data/ML/Train/{args.TN_filename}/RF/*df_p0.01_mod.txt; do python ./package/ML_classification_modified.py -df $inp -test $inp.test -cl_train pos,neg -alg $A -x_norm T -n 10 -cv_num 10 -plots T -n_jobs 30 -gs_type $G -gs_n 10 -tag x_norm_T_gs_type_$G; python ./package/ML_classification_modified.py -df $inp -test $inp.test -cl_train pos,neg -alg $A -x_norm F -n 10 -cv_num 10 -plots T -n_jobs 30 -gs_type $G -gs_n 10 -tag x_norm_F_gs_type_$G; done; done; done"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("---------------------------------ML COMPLETE-----------------------------\n")

    ##  files sorting
    print("----------------------------------sorting files----------------------------")
    command = f"python ./package/files_sorting.py -new_group_folder {new_group_folder}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("------------------------------files sorting COMPLETE------------------------")



main(gff3_filename=args.gff3_filename, genome_filename=args.genome_filename,TP_filename=args.TP_filename, TN_filename=args.TN_filename, set_num=int(args.set_num), new_group_folder=args.group_name)



    
