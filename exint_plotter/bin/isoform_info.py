#!/usr/bin/env python

import pandas as pd
import argparse  

#read arguments
parser = argparse.ArgumentParser(description="For each exon, add upstream and downstream phase")
parser.add_argument("--input_file", "-i", required=True)
parser.add_argument("--output_file1", "-o1", required=True)
parser.add_argument("--output_file2", "-o2", required=True)

args = parser.parse_args()
my_input = args.input_file
my_output_1 = args.output_file1
my_output_2 =  args.output_file2

#main
import pandas as pd
import re

my_gtf = pd.read_table(my_input, sep="\t", header=None, index_col=False)
#subset GTF only for CDS exons. Otherwise the last last exon will probably be non-coding.
my_gtf = my_gtf.loc[my_gtf.iloc[:,2]=="CDS"]

my_gtf_subset = my_gtf.iloc[:,my_gtf.shape[1]-1]

######### Exon number by isoform
#subset only the lines containing an exon number
my_gtf_subset = my_gtf_subset[my_gtf_subset.str.contains("exon_number")]
#building the output
my_raw_gene_id = [part for element in list(my_gtf_subset) for part in element.split(";") if "gene_id" in part]
my_gene_id = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_gene_id]
my_raw_transcript_id = [part for element in list(my_gtf_subset) for part in element.split(";") if "transcript_id" in part]
my_transcript_id = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_transcript_id]
my_raw_exon_num = [part for element in list(my_gtf_subset) for part in element.split(";") if "exon_number" in part]
my_exon_num = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_exon_num]
ex_num_df = pd.concat([pd.Series(my_gene_id), pd.Series(my_transcript_id), pd.Series(my_exon_num)], axis=1)
ex_num_df.columns = ["GeneID", "TranscriptID", "ExNum"]
ex_num_df = ex_num_df.drop_duplicates() #remove duplicated rows (exons with different splices sites)
ex_num_df["ExNum"] = [int(element) for element in list(ex_num_df["ExNum"])] #change data type of ExNum to integer
grouped_df = ex_num_df.groupby("TranscriptID")
final_df = pd.DataFrame()
for name, group in grouped_df:
  my_group = group.loc[group.ExNum == max(list(group["ExNum"]))]
  final_df = pd.concat([final_df, my_group])
final_df.to_csv(my_output_1, sep="\t", header=False, index=False) #save first output file

######## Exon number in isoform
my_chr = list(my_gtf.iloc[:,0])
my_start = list(my_gtf.iloc[:,3])
my_stop = list(my_gtf.iloc[:,4])
my_last = list(my_gtf.iloc[:,my_gtf.shape[1]-1]) #last field of the GTF
my_coords = [element+":"+str(element1)+"-"+str(element2) for element, element1, element2 in zip(my_chr, my_start, my_stop)]
#Generate another dataframe
data_tuples = list(zip(my_coords, my_last))
intermediate_df = pd.DataFrame(data_tuples, columns=["Coords", "Info"])
filtered_df = intermediate_df[intermediate_df["Info"].str.contains("exon_number")]
#In theory, the gene_id, transcript_id, exon_num lists are in the same order, and I can use them directly
filtered_df["GeneID"] = my_gene_id
filtered_df["TranscriptID"] = my_transcript_id
filtered_df["ExNum"] = my_exon_num
filtered_df = filtered_df.drop(columns=["Info"])
final_df = filtered_df[["Coords", "ExNum", "GeneID", "TranscriptID"]] #reorder the columns
final_df = final_df.drop_duplicates()
final_df.to_csv(my_output_2, sep="\t", header=False, index=False) #save first output file
