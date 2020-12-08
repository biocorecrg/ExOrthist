#!/usr/bin/env python

import argparse  
import pandas as pd

#read arguments
parser = argparse.ArgumentParser(description="Recluster the gene clusters by species pairs based on orthopairs")
parser.add_argument("--output_folder", "-o", required=True)
parser.add_argument("--input_file", "-i", required=True)
parser.add_argument("--aln_num", "-n", required=True)

args = parser.parse_args()
my_output_folder = str(args.output_folder)
my_input = args.input_file
my_aln_num = int(args.aln_num)

#### Main
#read input
ex_to_split_df = pd.read_table(my_input, header=0, index_col=False, sep="\t")
if ex_to_split_df.shape[0] == 0:
  ex_to_split_df.to_csv(my_output_folder+"/EXs_to_realign_part_1.txt", sep="\t", header=True, index=False, na_rep="NA")
else:
  #split in df of aln_num entries and save to file
  list_df = [ex_to_split_df[i:i+my_aln_num] for i in range(0,ex_to_split_df.shape[0],my_aln_num)]
  num = 0
  for df in list_df:
    num = num + 1
    df.to_csv(my_output_folder+"/EXs_to_realign_part_"+str(num)+".txt", sep="\t", header=True, index=False, na_rep="NA")
