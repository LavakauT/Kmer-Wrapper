import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
	'-set_num',
	help='How many files you need to do',
  required=True)
parser.add_argument(
	'-TN_filename',
	help='Your TN data file name(NOT PATH)',
  required=True)
args = parser.parse_args()

work_dir = "./Data"
print(os.getcwd())
for y in range(1,int(args.set_num)+1):
    #Train DF
    Test= pd.read_csv(work_dir+"/ML/Test/neg_data_Test_"+str(y)+".txt",sep = "\t")
    Train= pd.read_csv(work_dir+f"/ML/Train/{args.TN_filename}/neg_{args.TN_filename}_Train_"+str(y)+".fa.pcre_df_p0.01_mod.txt",sep = "\t")
    combine = pd.concat([Train,Test],join='inner')
    combine.to_csv(work_dir+f"/ML/Train/{args.TN_filename}/RF/neg_{args.TN_filename}_"+str(y)+".fa.pcre_df_p0.01_mod.txt",sep='\t',index=False)

    #Test list
    testlist = Test['X'].tolist()
    with open(work_dir+f"/ML/Train/{args.TN_filename}/RF/neg_{args.TN_filename}_"+str(y)+".fa.pcre_df_p0.01_mod.txt.test","w") as test_w:
        for gene in testlist:
            test_w.write(str(gene)+'\n')

