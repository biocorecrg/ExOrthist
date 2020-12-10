#!/usr/bin/env python

import pandas as pd
import argparse  

#several inputs to set here
parser = argparse.ArgumentParser(description="For each exon, add upstream and downstream phase")
parser.add_argument("--best_hits", "-b", required=True)
parser.add_argument("--chosenID_query", "-q", required=True)
parser.add_argument("--chosenID_target", "-t", required=True)
parser.add_argument("--output_file", "-out", required=True)

args = parser.parse_args()
my_best_hits = args.best_hits
my_chosenID_query = args.chosenID_query
my_chosenID_target = args.chosenID_target
my_output = args.output_file

######## main
import pandas as pd
import re

my_df =  pd.read_csv(str(my_best_hits), sep="\t", index_col=False, header=0, names=['CID', 'Feature', 'Exon_prot_loc', 'GeneID_query', 'Exon_number_query', 'ExonID_query', 'GeneID_target', 'Exon_number_target', 'ExonID_target', 'Score_C1', 'Score_I1', 'Score_A', 'Score_I2', 'Score_C2', 'Total_exon_score', 'Species_query', 'Species_target'])
#select and reorder the columns
my_subsetted_df = my_df.loc[:,['CID', 'Species_query', 'GeneID_query', 'ExonID_query', 'Species_target', 'GeneID_target', 'ExonID_target', 'Score_C1', 'Score_I1', 'Score_A', 'Score_I2', 'Score_C2', 'Total_exon_score']]
#remove the protein isoform
my_subsetted_df["GeneID_query"] = pd.Series([re.sub(".*\|", "", element) for element in list(my_subsetted_df.GeneID_query)])
my_subsetted_df["GeneID_target"] = pd.Series([re.sub(".*\|", "", element) for element in list(my_subsetted_df.GeneID_target)])
#remove the strand information from the exon coords
my_subsetted_df["ExonID_query"] = pd.Series([re.sub(":-", "", re.sub(":\+", "", str(element))) for element in list(my_subsetted_df.ExonID_query)])
my_subsetted_df["ExonID_target"] =  pd.Series([re.sub(":-", "", re.sub(":\+", "", str(element))) for element in list(my_subsetted_df.ExonID_target)])

#read dictionary ExonID - ChosenID for each overlapping group
my_query_chosen_df = pd.read_table(my_chosenID_query, sep="\t", header=0, index_col=False) #query species
query_ExonID_ChosenID_dict = pd.Series(my_query_chosen_df.ChosenID.values, index=my_query_chosen_df.ExonID).to_dict()
my_target_chosen_df = pd.read_table(my_chosenID_target, sep="\t", header=0, index_col=False) #target species
target_ExonID_ChosenID_dict = pd.Series(my_target_chosen_df.ChosenID.values, index=my_target_chosen_df.ExonID).to_dict()
#operate translation
my_subsetted_df["OverlappingID_query"] = my_subsetted_df["ExonID_query"].map(query_ExonID_ChosenID_dict)
my_subsetted_df["OverlappingID_target"] = my_subsetted_df["ExonID_target"].map(target_ExonID_ChosenID_dict)

#reorder columns
my_subsetted_df = my_subsetted_df[['CID', 'Species_query', 'GeneID_query', 'ExonID_query', 'OverlappingID_query', 'Species_target', 'GeneID_target', 'ExonID_target', 'OverlappingID_target', 'Score_C1', 'Score_I1', 'Score_A', 'Score_I2', 'Score_C2', 'Total_exon_score']]
#write to file
my_subsetted_df.to_csv(str(my_output), sep="\t", index=False, header=True, na_rep="NA")
