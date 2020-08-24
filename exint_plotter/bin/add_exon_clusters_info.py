#!/usr/bin/env python

import pandas as pd
import argparse  

#read arguments
parser = argparse.ArgumentParser(description="Add info regarding evolutionary conservation to each exon")
parser.add_argument("--input_file", "-i", required=True)
parser.add_argument("--output_file", "-o", required=True)
parser.add_argument("--exon_clusters", "-c", required=True)
parser.add_argument("--species", "-s", required=True)

args = parser.parse_args()
my_input = args.input_file
my_output = args.output_file
my_clusters = args.exon_clusters
my_species = args.species

#### Main
import pandas as pd
import re

#read input
my_df = pd.read_table(str(my_input), header=None, names=["GeneID", "ExonID", "Strand", "UpPhase", "DownPhase"])
exon_clusters_df = pd.read_table(str(my_clusters), header=0, index_col=False, sep="\t")
exon_clusters_df_species = exon_clusters_df.loc[exon_clusters_df["Sps"] == my_species]

#create dictionary with key=exon_coords, value=exon_cluster_id
#it also considers only start or only stop coords
coords_cluster_dict = pd.Series(exon_clusters_df_species.ClusterID.values, index=exon_clusters_df_species.Exon_coords).to_dict()
for my_id in list(coords_cluster_dict.keys()):
  start_coord = re.sub('-.*', '', my_id)
  stop_coord = re.sub(':.*-', ':', my_id)
  coords_cluster_dict[start_coord] = coords_cluster_dict[my_id]
  coords_cluster_dict[stop_coord] = coords_cluster_dict[my_id]

#get exon_cluster_id
my_final_dict = {}
for element in list(my_df["ExonID"]): #for each exon_coord
  if element in list(coords_cluster_dict.keys()): #for each ex_coords in the exon_cluster       
    my_final_dict[element] = coords_cluster_dict[element] #save clusterID as value for ex_coord key
  else:
    element_start = re.sub('-.*', '', element) #select start and stop
    if element_start in list(coords_cluster_dict.keys()):
      my_final_dict[element] = coords_cluster_dict[element_start]
    else:
      element_stop = re.sub(':.*-', ':', element)
      if element_stop in list(coords_cluster_dict.keys()):
        my_final_dict[element] = coords_cluster_dict[element_stop]
      else:
        my_final_dict[element] = element


my_df["ClusterID"] = [my_final_dict[my_id] for my_id in list(my_df["ExonID"])]
my_df["Cluster_status"] = ["not_in_clusters"  if x not in list(coords_cluster_dict.values()) else "Conserved" for x in list(my_df["ClusterID"])]

#write pandas df to file
my_df.to_csv(str(my_output), sep="\t", header=True, index=False)
