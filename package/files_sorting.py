import os
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
	'-new_group_folder',
	help='New group folder neme',
    required=True)

args = parser.parse_args()

source_directories = ["./Data", "./gffandgenome"]
destination_directory = f"./Results/{args.new_group_folder}"
for source_directory in source_directories:
    print(source_directory)
    try:
        dir_name = os.path.basename(source_directory.rstrip('/'))
        dest_dir = os.path.join(destination_directory, dir_name)
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

        for item in os.listdir(source_directory):
            s = os.path.join(source_directory, item)
            d = os.path.join(dest_dir, item)
            if os.path.isdir(s):
                shutil.copytree(s, d)
            else:
                shutil.copy2(s, d)
        print(f"Directory {source_directory} copied to {dest_dir}")
    except Exception as e:
        print(f"Error: {e}")

if (os.path.isdir(f"./Results/{args.new_group_folder}/Data") == True) and (os.path.isdir(f"./Results/{args.new_group_folder}/gffandgenome") == True) :
    shutil.rmtree("./Data")
    shutil.rmtree("./gffandgenome")

build_pathes = ["./Data/Kmer/Test", "./Data/Kmer/TPTN", "./Data/Kmer/Train", "./Data/ML/Test", "./gffandgenome"]
for build_path in build_pathes :
    os.makedirs(build_path)
