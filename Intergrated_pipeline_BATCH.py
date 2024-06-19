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
	'-set_num',
	help='How many test/train set you want to split --> 1~40',
  default=1)
parser.add_argument(
	'-group_name',
	help='Customize your output folder name',
  required=True)

args = parser.parse_args()

import subprocess
import os
import time

dir_this_funtion_at = args.main_dir
os.chdir(dir_this_funtion_at)
print("Current working directory:" + os.getcwd())
new_group_folder = os.getcwd() + "/Results/" + str(args.group_name)

def get_suffix(filename):
    return filename.split("_")[1:2]

command = f"python ./package/Create_group_folder.py  -new_group_name {args.group_name}"
subprocess.run("bash -c '" + command + "'", shell=True, check=True)
TPTNs = os.listdir(f"./Results/{args.group_name}/Data/Kmer/TPTN")
TPs = []
TNs = []

for file in TPTNs:
    if file.startswith("TP_"):
        TPs.append(file)
    elif file.startswith("TN_"):
        TNs.append(file)
    else:
        print("Please put files which start with 'TP_' or 'TN_' !!!")

for TP_filename in TPs:
    for TN_filename in TNs:
        TP_suffix = get_suffix(TP_filename)
        TN_suffix = get_suffix(TN_filename)
        if TP_suffix == TN_suffix:
          print(f"Now processing {TP_filename} and {TN_filename}")
          time.sleep(0.5)
          command = f"python {args.main_dir}/Intergrated_pipeline.py -main_dir {args.main_dir} -conda_dir {args.conda_dir} -gff3_filename {args.gff3_filename} -features {args.features} -genome_filename {args.genome_filename} -TP_filename {TP_filename} -TN_filename {TN_filename} -set_num {args.set_num} -group_name {args.group_name}"
          subprocess.run("bash -c '" + command + "'", shell=True)