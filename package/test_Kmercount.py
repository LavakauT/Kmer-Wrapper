import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.stats import fisher_exact
from os import chdir
import argparse

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

wkdir = './Data'
chdir(wkdir)


for y in range(1, int(args.set_num)+1):
  df_dir = f"ML/Train/{args.TN_filename}/neg_{args.TN_filename}_"+str(y)+"_distinct_pcc_enriched_kmer.txt"
  POS = f"Kmer/Test/{args.TP_filename}_Test_"+str(y)+".txt.fa"
  NEG = f"Kmer/Test/{args.TN_filename}_Test_"+str(y)+".txt.fa"
  df = pd.read_csv(df_dir,sep = "\t")
  motifs = df['motif'].tolist()

  p = open(POS,'r')
  n = open(NEG,'r')
  num_pos = 0
  num_neg = 0
  genes = ['Skip_this_line']
  numpy_header = ['Class']
  for i in motifs:
      numpy_header.append(i)  

  dataframe = np.zeros([1,len(motifs)+1]) 
  positive_present = {}.fromkeys(motifs, 0) 
  negative_present = {}.fromkeys(motifs, 0)

  for seq_record in SeqIO.parse(p, 'fasta'):
      num_pos += 1
      header = seq_record.id
      genes.append(header)
      seq = str(seq_record.seq)
      gene_array =np.array([1])       # Array of P/A (1/0) for each gene - starts with '1' For Positive Class
      for ki in motifs:
        if " " in ki:                         #Checks to see if motif is a pair - pairs are separated by a space
          k1 = Seq(ki.split(" ")[0])
          k2 = Seq(ki.split(" ")[1])
          if str(k1) in seq or str(k1.reverse_complement()) in seq and str(k2) in seq or str(k2.reverse_complement()) in seq: 
            gene_array = np.append(gene_array, 1)
            positive_present[ki] = positive_present[ki]+1
          else:
            gene_array = np.append(gene_array, 0)

        else:                                #If no separation by a space, assumes you're looking at singletons.
          kmer = Seq(ki)
          if str(kmer) in seq or str(kmer.reverse_complement()) in seq:
            gene_array = np.append(gene_array, 1)
            positive_present[ki] = positive_present[ki]+1
          else:
            gene_array = np.append(gene_array, 0)
      dataframe = np.vstack((dataframe,gene_array))

  for seq_record in SeqIO.parse(n, 'fasta'):
      num_neg += 1
      header = seq_record.id
      genes.append(header)
      seq = str(seq_record.seq)
      gene_array =np.array([0])       # Array of P/A (1/0) for each gene - starts with '1' For Positive Class
      for ki in motifs:
        if " " in ki:                         #Checks to see if motif is a pair - pairs are separated by a space
          k1 = Seq(ki.split(" ")[0])
          k2 = Seq(ki.split(" ")[1])
          if str(k1) in seq or str(k1.reverse_complement()) in seq and str(k2) in seq or str(k2.reverse_complement()) in seq: 
            gene_array = np.append(gene_array, 1)
            positive_present[ki] = positive_present[ki]+1
          else:
            gene_array = np.append(gene_array, 0)

        else:                                #If no separation by a space, assumes you're looking at singletons.
          kmer = Seq(ki)
          if str(kmer) in seq or str(kmer.reverse_complement()) in seq:
            gene_array = np.append(gene_array, 1)
            positive_present[ki] = positive_present[ki]+1
          else:
            gene_array = np.append(gene_array, 0)
      dataframe = np.vstack((dataframe,gene_array))

  DF= pd.DataFrame(dataframe, index=genes, columns=numpy_header, dtype=int)  # , dtype=int # Turn numpy into pandas DF
  DF= DF.drop("Skip_this_line",0)
  print(DF)
  DF.to_csv("ML/Test/data_Test_"+str(y)+".txt", sep='\t')