#!/usr/bin/env python

import argparse  
import pandas as pd
import re

#read arguments
parser = argparse.ArgumentParser(description="Recluster the gene clusters by species pairs based on orthopairs")
parser.add_argument("--orthopairs", "-op", required=True)
parser.add_argument("--orthogroups", "-og", required=True)
parser.add_argument("--species1", "-s1", required=True)
parser.add_argument("--species2", "-s2", required=True)
parser.add_argument("--output_file", "-out", required=True)

args = parser.parse_args()
my_orthopairs = args.orthopairs
my_orthogroup = args.orthogroups
species1 = args.species1
species2 = args.species2
my_output = args.output_file

#### Main
#read input
orthopairs_df = pd.read_table(my_orthopairs, sep="\t", header=0, names=["GeneID1", "GeneID2"])
orthopairs_df["GeneID1"] = [re.sub(".*\\|", "", element) for element in list(orthopairs_df["GeneID1"])]
orthopairs_df["GeneID2"] = [re.sub(".*\\|", "", element) for element in list(orthopairs_df["GeneID2"])]
orthogroups_df = pd.read_table(my_orthogroup, sep="\t", header=0, names=["ClusterID", "Species", "GeneID", ""], index_col=None)
#add Species to orthopairs
geneID_species_dict = pd.Series(orthogroups_df.Species.values, index=orthogroups_df.GeneID).to_dict()
orthopairs_df["Species1"] = orthopairs_df["GeneID1"].map(geneID_species_dict)
orthopairs_df["Species2"] = orthopairs_df["GeneID2"].map(geneID_species_dict)
#add ClusterID to orthopairs
geneID_clusterID_dict = pd.Series(orthogroups_df.ClusterID.values, index=orthogroups_df.GeneID).to_dict()
orthopairs_df["ClusterID"] = orthopairs_df["GeneID1"].map(geneID_clusterID_dict) #GeneID1 and GeneID2 would do the same job here.
#filter only for species1 and species2 orthologs.
species_pair_df = orthopairs_df.loc[orthopairs_df.Species1.isin([species1, species2])]
species_pair_df = species_pair_df.loc[orthopairs_df.Species2.isin([species1, species2])]
#group by gene cluster
species_pair_grouped_df = species_pair_df.groupby("ClusterID")
final_df = pd.DataFrame()
for name, group in species_pair_grouped_df:
  recluster_ID=1
  cluster_genes_list = list(set(list(group["GeneID1"])+list(group["GeneID2"])))
  my_cluster_counter = 1
  while my_cluster_counter == 1:
    focus_gene = cluster_genes_list[0]
    subgroup = group.loc[(group.GeneID1==focus_gene) | (group.GeneID2==focus_gene)]
    subgroup_genes = list(set(list(subgroup["GeneID1"])+list(subgroup["GeneID2"])))
    my_subcluster_counter = 1
    while my_subcluster_counter == 1:
      initial_len = len(subgroup_genes)
      subgroup1 = group.loc[(group.GeneID1.isin(subgroup_genes) | group.GeneID2.isin(subgroup_genes))]  
      subgroup1_genes = list(set(list(subgroup1["GeneID1"])+list(subgroup1["GeneID2"])))
      subgroup_genes = list(set(subgroup_genes+subgroup1_genes))
      if len(subgroup_genes) == initial_len:
        my_subcluster_counter = 2
    #generate subclusterID
    if len(str(recluster_ID)) == 1:
      subcluster_ID = list(group["ClusterID"])[0]+".0"+str(recluster_ID)
    else:
      subcluster_ID = list(group["ClusterID"])[0]+"."+str(recluster_ID)
    #add sub-cluster to final dataframe
    subgroup_df = pd.DataFrame({"ClusterID" : subcluster_ID, "GeneID" : subgroup_genes})
    final_df = pd.concat([final_df, subgroup_df])
    #update original list
    cluster_genes_list = [element for element in cluster_genes_list if element not in subgroup_genes]
    #update recluster_ID
    recluster_ID = recluster_ID+1
    if len(cluster_genes_list) == 0:
      my_cluster_counter = 2

#add species
final_df["Species"] = final_df["GeneID"].map(geneID_species_dict)
#reorder columns
final_df = final_df[["ClusterID","Species","GeneID"]]
#save final df to file
final_df.to_csv(my_output, sep="\t", index=False, header=False, na_rep="NA")
