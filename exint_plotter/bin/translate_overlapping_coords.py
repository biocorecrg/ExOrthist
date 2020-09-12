#!/usr/bin/env python

import pandas as pd
import argparse  

#several inputs to set here
parser = argparse.ArgumentParser(description="For each exon, add upstream and downstream phase")
parser.add_argument("--best_hits", "-b", required=True)
parser.add_argument("--overlapping_query", "-q", required=True)
parser.add_argument("--overlapping_target", "-t", required=True)
parser.add_argument("--selected_overlapping_query", "-sq", required=True)
parser.add_argument("--selected_overlapping_target", "-st", required=True)
parser.add_argument("--output_file", "-out", required=True)

args = parser.parse_args()
my_best_hits = args.best_hits
my_overlapping_query = args.overlapping_query
my_overlapping_target = args.overlapping_target
my_selected_overlapping_query = args.selected_overlapping_query
my_selected_overlapping_target = args.selected_overlapping_target
my_output = args.output_file

######## main
import pandas as pd
import re

my_df =  pd.read_csv(str(my_best_hits), sep="\t", index_col=False, header=0, names=['CID', 'Feature', 'Exon_prot_loc', 'GeneID_query', 'Exon_number_query', 'ExonID_query', 'GeneID_target', 'Exon_number_target', 'ExonID_target', 'Score_C1', 'Score_I1', 'Score_A', 'Score_I2', 'Score_C2', 'Total_exon_score', 'Species_query', 'Species_target'])

my_subsetted_df = my_df.loc[:,['CID', 'Species_query', 'GeneID_query', 'ExonID_query', 'Species_target', 'GeneID_target', 'ExonID_target', 'Score_C1', 'Score_I1', 'Score_A', 'Score_I2', 'Score_C2', 'Total_exon_score']] #selecting and reordering the columns

#modify the input, since I am at this:
#remove the protein isoform
my_subsetted_df["GeneID_query"] = pd.Series([re.sub(".*\|", "", element) for element in list(my_subsetted_df.GeneID_query)])
my_subsetted_df["GeneID_target"] = pd.Series([re.sub(".*\|", "", element) for element in list(my_subsetted_df.GeneID_target)])

#remove the strand information from the exon coords
my_subsetted_df["ExonID_query"] = pd.Series([re.sub(":-", "", re.sub(":\+", "", str(element))) for element in list(my_subsetted_df.ExonID_query)])
my_subsetted_df["ExonID_target"] =  pd.Series([re.sub(":-", "", re.sub(":\+", "", str(element))) for element in list(my_subsetted_df.ExonID_target)])

#get dictionary for query and target species with key=ExonID, value=OverlappingID (meaning choords)
query_overlapping = pd.read_csv(str(my_overlapping_query), sep="\t", index_col=False, header=None, names=["GeneID", "OverlappingID", "ExonID", "Frequency", "Length"])
target_overlapping = pd.read_csv(str(my_overlapping_target), sep="\t", index_col=False, header=None, names=["GeneID", "OverlappingID", "ExonID", "Frequency", "Length"])
query_overlapping_dict = pd.Series(query_overlapping.OverlappingID.values, index=query_overlapping.ExonID).to_dict()
target_overlapping_dict = pd.Series(target_overlapping.OverlappingID.values, index=target_overlapping.ExonID).to_dict()

#add column with overlappingID for both query and target
my_subsetted_df["OverlappingID_query"] = my_subsetted_df["ExonID_query"].map(query_overlapping_dict) #map is much faster than replace
my_subsetted_df["OverlappingID_target"] = my_subsetted_df["ExonID_target"].map(target_overlapping_dict)
#get dictionary for query and target species with key=OverlappingID, value=ExonID
#it is the ExonID of the exon representing each overlapping group.
query_chosen_overlapping = pd.read_csv(str(my_selected_overlapping_query), sep="\t", index_col=False, header=0, names=["GeneID", "OverlappingID", "ExonID", "Frequency", "Length"])
target_chosen_overlapping = pd.read_csv(str(my_selected_overlapping_target), sep="\t", index_col=False, header=0, names=["GeneID", "OverlappingID", "ExonID", "Frequency", "Length"])
query_chosen_dict = pd.Series(query_chosen_overlapping.ExonID.values, index=query_chosen_overlapping.OverlappingID).to_dict()
target_chosen_dict = pd.Series(target_chosen_overlapping.ExonID.values, index=target_chosen_overlapping.OverlappingID).to_dict()

#replace the OverlappingIDs with the relative ExonIDs
my_subsetted_df["OverlappingID_query"] = my_subsetted_df["OverlappingID_query"].map(query_chosen_dict)
my_subsetted_df["OverlappingID_target"] = my_subsetted_df["OverlappingID_target"].map(target_chosen_dict)

#reorder columns
my_subsetted_df = my_subsetted_df[['CID', 'Species_query', 'GeneID_query', 'ExonID_query', 'OverlappingID_query', 'Species_target', 'GeneID_target', 'ExonID_target', 'OverlappingID_target', 'Score_C1', 'Score_I1', 'Score_A', 'Score_I2', 'Score_C2', 'Total_exon_score']]

#write to file
my_subsetted_df.to_csv(str(my_output), sep="\t", index=False, header=True, na_rep="NA")
