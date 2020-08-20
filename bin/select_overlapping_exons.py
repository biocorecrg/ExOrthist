#!/usr/bin/env python

import pandas as pd
import argparse  

#read arguments
parser = argparse.ArgumentParser(description="Select only one exon from a set of overlapping ones")
parser.add_argument("--input_file", "-i", required=True)
parser.add_argument("--output_file", "-o", required=True)

args = parser.parse_args()
my_input = args.input_file
my_output = args.output_file

#read  input table
my_df = pd.read_table(str(my_input), header=None, names=["geneID", "ExOverlapID", "ExonID", "frequency", "length"])
my_df["frequency"] = pd.to_numeric(my_df["frequency"], errors="coerce") #transform frequency in the right data type
my_df = my_df.fillna(0)
 
final_df = pd.DataFrame(columns=["geneID", "ExOverlapID", "ExonID", "frequency", "length"])
#group by overlap ID
my_df = my_df.groupby("ExOverlapID")
for name, group in my_df:
  all_freq_list = list(group["frequency"])
  max_freq = max(all_freq_list)
  selected_elements_df = group.loc[group.frequency==max_freq]
                   
  if selected_elements_df.shape[0] > 1:
    selected_elements_df = selected_elements_df.loc[selected_elements_df.length==max(list(selected_elements_df["length"]))]
 
  final_df = final_df.append(selected_elements_df, ignore_index=True) # will add the selected element to the final dataframe.
 
#save  final df to file
final_df.to_csv(str(my_output), sep="\t", index=False, index_label=False)