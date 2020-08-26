#!/usr/bin/env python

import argparse
import pandas as pd
import collections

#Read arguments
parser = argparse.ArgumentParser(description="Subset input files to only genes in a given gene cluster")
parser.add_argument("--index", "-i", required=True)
parser.add_argument("--tot", "-t", required=True)
parser.add_argument("--gene_clusters", "-g", required=True)
parser.add_argument("--output", "-out", required=True)

args = parser.parse_args()
my_index_file = args.index
my_tot_file = args.tot
my_gene_clusters_file = args.gene_clusters
my_output_file = args.output

### Main
my_ex_index_df = pd.read_table(my_index_file, sep="\t", header=None, names=["Coords", "Pos", "GeneID", "TranscriptID"])
my_ex_tot_df = pd.read_table(my_tot_file, sep="\t", header=None, names=["GeneID", "TranscriptID", "TotEx"])

my_ex_index_df["Pos"] = [int(element) for element in list(my_ex_index_df["Pos"])] #make sure that the pos in an integer
my_ex_index_df_grouped = my_ex_index_df.groupby("TranscriptID")
my_first_ex_list = []
for name, group in my_ex_index_df_grouped:
  my_first_ex_list.append(list(group.loc[group.Pos==min(list(group["Pos"]))]["Coords"])[0]) #select the first element of the list

my_last_ex_dict = pd.Series(my_ex_tot_df.TotEx.values, index=my_ex_tot_df.TranscriptID).to_dict()
my_ex_index_df["TotEx"] = my_ex_index_df["TranscriptID"].map(my_last_ex_dict)
my_last_ex_list = list(my_ex_index_df.loc[my_ex_index_df["Pos"]==my_ex_index_df["TotEx"]]["Coords"])
my_internal_ex_list = [element for element in list(my_ex_index_df["Coords"]) if element not in my_first_ex_list+my_last_ex_list]
my_first_last_ex_raw = list(set(my_first_ex_list))+list(set(my_last_ex_list))
my_first_last_ex_list = [element for element, count in collections.Counter(my_first_last_ex_raw).items() if count > 1]

#create dictionary for exon status. This works because update will replace the last value.
my_first_ex_dict = {element : "first" for element in my_first_ex_list}
my_last_ex_dict = {element : "last" for element in my_last_ex_list}
my_firstlast_ex_dict = {element : "first;last" for element in my_first_last_ex_list}
my_internal_ex_dict = {element : "Internal" for element in my_internal_ex_list}
my_all_ex_dict = {**my_first_ex_dict, **my_last_ex_dict, **my_firstlast_ex_dict, **my_internal_ex_dict}

#create dictionary for gene clusters
my_clusters = pd.read_table(my_gene_clusters_file, sep="\t", header=None)
my_gene_clusters_dict = pd.Series(list(my_clusters.iloc[:,0]), index=list(my_clusters.iloc[:,2])).to_dict()

#Join results in final table
my_ex_index_df["Status"] = my_ex_index_df["Coords"].map(my_all_ex_dict)
my_ex_index_df["ClusterID"] = my_ex_index_df["GeneID"].map(my_gene_clusters_dict)
my_final_df = my_ex_index_df[["GeneID", "ClusterID", "Coords", "Status"]].drop_duplicates()
my_final_df.to_csv(my_output_file, sep="\t", header=False, index=False, na_rep="NA")
