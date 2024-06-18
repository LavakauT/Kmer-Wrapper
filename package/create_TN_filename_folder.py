import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
	'-TN_filename',
	help='Your TN data file name(NOT PATH)',
  required=True)
args = parser.parse_args()

build_pathes = [f"./Data/ML/Train/{args.TN_filename}/RF", f"./Data/Kmer/Train/{args.TN_filename}"]

for build_path in build_pathes:
	os.makedirs(build_path)