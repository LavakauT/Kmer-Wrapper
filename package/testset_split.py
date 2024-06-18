import pandas as pd
from sklearn.model_selection import train_test_split
import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument(
	'-set_num',
	help='How many files you need to do',
  required=True)
parser.add_argument(
	'-TP_filename',
	help='Your TP data file name(NOT PATH)',
  required=True)
parser.add_argument(
	'-TN_filename',
	help='Your TN data file name(NOT PATH)',
  required=True)

args = parser.parse_args()
work_dir = "./Data"

def split_test_train_set(work_dir=work_dir, set_num=int(args.set_num), test_size=0.2):
    TP= pd.read_csv(work_dir+f"/Kmer/TPTN/{args.TP_filename}.txt",sep = "\t",names=['x'])
    TN= pd.read_csv(work_dir+f"/Kmer/TPTN/{args.TN_filename}.txt",sep = "\t",names=['x'])
    TP=TP['x'].tolist()
    TN=TN['x'].tolist()

    min_samples = min(len(TP), len(TN))
    if len(TP) > min_samples:
        TP = random.sample(TP, min_samples)
        TN = TN
    elif len(TN) > min_samples:
        TN = random.sample(TN, min_samples)
        TP = TP
      
    for i in range(set_num):
        TP_Test = []
        TN_Test = []
        TP_Train = []
        TN_Train = []
        TP_Train,TP_Test,TN_Train,TN_Test = train_test_split(TP,TN,random_state=i,test_size=test_size,shuffle=True)
        with open(work_dir+f"/Kmer/Train/{args.TP_filename}_Train_"+str(i+1)+".txt","w") as TP_Train_w:
            for gene in TP_Train:
                TP_Train_w.write(str(gene)+"\n")
        with open(work_dir+f"/Kmer/Train/{args.TN_filename}/{args.TN_filename}_Train_"+str(i+1)+".txt","w") as TN_Train_w:
            for gene in TN_Train:
                TN_Train_w.write(str(gene)+"\n")
        with open(work_dir+f"/Kmer/Test/{args.TP_filename}_Test_"+str(i+1)+".txt","w") as TP_Test_w:
            for gene in TP_Test:
                TP_Test_w.write(str(gene)+"\n")
        with open(work_dir+f"/Kmer/Test/{args.TN_filename}_Test_"+str(i+1)+".txt","w") as TN_Test_w:
            for gene in TN_Test:
                TN_Test_w.write(str(gene)+"\n")

split_test_train_set()