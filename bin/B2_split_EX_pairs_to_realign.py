#!/usr/bin/env python

import argparse  
import pandas as pd
import re

#read arguments
parser = argparse.ArgumentParser(description="Split the pool of exon pairs to be realigned into batches of aln_num pairs")
parser.add_argument("--input_file", "-i", required=True)
parser.add_argument("--aln_num", "-n", required=True)

args = parser.parse_args()
my_input = args.input_file
my_aln_num = int(args.aln_num)

#### Main
#read input
ex_to_split_df = pd.read_table(my_input, header=0, index_col=False, sep="\t")
output_suffix = re.sub(".txt", "_", re.sub("split", "realign", str(my_input))) #It changes ${sp1}-${sp2}_EXs_to_split_part_*.txt to ${sp1}-${sp2}_EXs_to_realign_part_*_ 

if ex_to_split_df.shape[0] == 0:
  ex_to_split_df.to_csv(output_suffix+"1.txt", sep="\t", header=True, index=False, na_rep="NA") #modified 05/03/21
else:
  #split in df of aln_num entries and save to file
  list_df = [ex_to_split_df[i:i+my_aln_num] for i in range(0,ex_to_split_df.shape[0],my_aln_num)]
  num = 0
  for df in list_df:
    num = num + 1
    df.to_csv(output_suffix+str(num)+".txt", sep="\t", header=True, index=False, na_rep="NA") #modified 05/03/21
