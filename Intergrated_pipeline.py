import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
	'-main_dir',
	help='Direction of the main folder (code_intergration)',
  required=True)
parser.add_argument(
	'-conda_dir',
	help='Direction of your conda.sh',
  required=True)
parser.add_argument(
	'-gff3_filename',
	help='Your .gff3 file name(NOT PATH)',
  required=True)
parser.add_argument(
	'-features',
	help='Features from gff3 file that you want',
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
new_group_folder = os.getcwd() + "/Results/" + str(args.group_name)


def main(conda_dir, gff3_filename, features, genome_filename,TP_filename, TN_filename, set_num, new_group_name):
    sub_group_name = TP_filename + "__X__"  + TN_filename
    #  create ./Data/ML/Train/...
    command = f"python ./package/Create_subgroup_folder.py -sub_group_name {sub_group_name} -new_group_name {new_group_name}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)

    ##  (all) gff3 ---> 5utr.coord
    print("--------------------.coord file generating--------------------")
    gff3_path = f"{new_group_folder}/gffandgenome/" + gff3_filename
    command = f"python ./package/FastaManager_modified_Ray.py -f gff_prom_to_coord_5utr -gff {gff3_path} -features {features}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("----------------.coord file generation COMPLETE----------------\n")
    
    ##  (all) 5utr.coord ---> 5utr.coord.fa
    print("-------------------.coord.fa file generating-------------------")
    coords_path = gff3_path + "_prom-5utr.coord"
    genome_path = f"{new_group_folder}/gffandgenome/" + genome_filename
    command = f"python ./package/FastaManager_modified_Ray.py -f get_stretch4 -coords {coords_path} -fasta {genome_path}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("---------------.coord.fa file generation COMPLETE---------------\n")

    ##  (target) split test_set & train_set
    print("-----------------test_set & train_set splitting-----------------")
    command = f"python ./package/testset_split.py -set_num {set_num} -TP_filename {TP_filename} -TN_filename {TN_filename} -new_group_folder {new_group_folder} -sub_group_name {sub_group_name}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("-------------test_set & train_set splitting COMPLETE-------------\n")

    ##  (target) test/train_set.txt ---> txt.fasta
    print("-----------------target genes getting seqence-------------------")
    for folder_path in [f"{new_group_folder}/Data/Kmer/Test", f"{new_group_folder}/Data/Kmer/Train", f"{new_group_folder}/Data/Kmer/Train/{sub_group_name}"]:
        all_files = os.listdir(folder_path)
        txt_files = [file for file in all_files if file.endswith(".txt")]
        for file_name in txt_files:
            file_path = os.path.join(folder_path, file_name)
            command = f"python ./package/FastaManager_modified_Ray.py -f getseq2 -fasta {coords_path}.fa -name {file_path}"
            subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("----------------target genes get seqence COMPLETE----------------\n")

    ##  (target) making dataframe
    print("---------------------------making dataframe-----------------------")
    command = "for inp in " + new_group_folder + "/Data/Kmer/Train/" + sub_group_name + "/*.txt.fa; do neg=$inp; pos=" + new_group_folder + "/Data/Kmer/Train/${inp##*/}; pos=${pos/" + TN_filename + "/" + TP_filename + "}; python ./package/pCRE_Finding_FET.py -pos $pos -neg $neg -k ./package/6mer.txt FDR Y -save \"$inp.pcre\"; done"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("-----------------------making dataframe COMPLETE------------------\n")

    ##  (target) pos_neg_convertion
    print("--------------------------pos_neg_convertion----------------------")
    command = f"source {conda_dir} && conda deactivate && conda activate R_Kmer && Rscript ./package/pos_neg_conversion.R -TN_filename {TN_filename} -new_group_folder {new_group_folder} -sub_group_name {sub_group_name}"
    subprocess.run("bash -c '" + command + "'", shell=True)
    print("----------------------pos_neg_convertion COMPLETE------------------\n")

    ##  (target) pcc filtering
    print("------------------------------pcc_filtering--------------------------")
    command = f"source {conda_dir} && conda deactivate && conda activate R_Kmer && Rscript ./package/pcc_filtering.R -set_num {set_num} -new_group_folder {new_group_folder}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("--------------------------pcc_filtering COMPLETE----------------------\n")

    ##  (target) ML_preprocess
    print("-----------------------------ML_preprocessing-------------------------")
    command = f"for inp in {new_group_folder}/Data/ML/Train/{sub_group_name}/*.txt.fa.pcre_df_p0.01.txt; do python ./package/ML_preprocess.py -df $inp; done"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("-------------------------ML_preprocessing COMPLETE---------------------\n")

    ##  (target) generate ML test set
    print("---------------------------ML test_set generating-----------------------")
    command = f"python ./package/test_Kmercount.py -set_num {set_num} -TP_filename {TP_filename} -TN_filename {TN_filename} -new_group_folder {new_group_folder} -sub_group_name {sub_group_name}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("-----------------------ML test_set generation COMPLETE-------------------\n")

    ## (target)  pos_neg_convertion (ML_test)
    print("-------------------------pos_neg_convertion(ML_test)---------------------")
    command = f"source {conda_dir} && conda deactivate && conda activate R_Kmer && Rscript ./package/pos_neg_conversion_Testset.R -new_group_folder {new_group_folder}"
    subprocess.run("bash -c '" + command + "'", shell=True)
    print("---------------------pos_neg_convertion(ML_test) COMPLETE-----------------\n")

    ##  (target) pre_ML
    print("--------------------------------doing pre_ML----------------------------")
    command = f"python ./package/pre_ML.py -set_num {set_num} -TN_filename {TN_filename} -new_group_folder {new_group_folder} -sub_group_name {sub_group_name}"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("-------------------------------pre_ML COMPLETE---------------------------\n")
    
    ##  (target) ML
    print("-----------------------------------doing ML-----------------------------")
    command = f"for A in SVM LogReg RF; do for G in full random; do for inp in {new_group_folder}/Data/ML/Train/{sub_group_name}/ML/*df_p0.01_mod.txt; do python ./package/ML_classification_modified.py -df $inp -test $inp.test -cl_train pos,neg -alg $A -x_norm T -n 10 -cv_num 10 -plots T -n_jobs 30 -gs_type $G -gs_n 10 -tag x_norm_T_gs_type_$G; python ./package/ML_classification_modified.py -df $inp -test $inp.test -cl_train pos,neg -alg $A -x_norm F -n 10 -cv_num 10 -plots T -n_jobs 30 -gs_type $G -gs_n 10 -tag x_norm_F_gs_type_$G; done; done; done"
    subprocess.run("bash -c '" + command + "'", shell=True, check=True)
    print("---------------------------------ML COMPLETE-----------------------------\n")



main(conda_dir=args.conda_dir, gff3_filename=args.gff3_filename, features=args.features, genome_filename=args.genome_filename,TP_filename=args.TP_filename, TN_filename=args.TN_filename, set_num=int(args.set_num), new_group_name=args.group_name)