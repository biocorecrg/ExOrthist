#!/usr/bin/env python

import argparse  
import pandas as pd
import re

#read arguments
parser = argparse.ArgumentParser(description="Subset the exon clusters by species pairs based on the pairwise reclustered gene orthogroups")
parser.add_argument("--exon_pairs", "-ep", required=True)
parser.add_argument("--reclustered_genes", "-rg", required=True)
parser.add_argument("--exon_clusters", "-ec", required=True)
parser.add_argument("--species1", "-s1", required=True)
parser.add_argument("--species2", "-s2", required=True)
parser.add_argument("--output_file", "-out", required=True)

args = parser.parse_args()
my_exon_pairs = args.exon_pairs
my_reclustered_genes = args.reclustered_genes
my_exon_clusters = args.exon_clusters
species1 = args.species1
species2 = args.species2
my_output = args.output_file

###### Main
#read input
reclustered_genes_df = pd.read_table(my_reclustered_genes, sep="\t", header=None, names=["ClusterID", "Species", "GeneID"])
exon_pairs_df = pd.read_table(my_exon_pairs, sep="\t", header=None, names=["GeneID1", "ExonID1", "GeneID2", "ExonID2", "Species1", "Species2"])
exon_clusters_df = pd.read_table(my_exon_clusters, sep="\t", header=0) #header=[ExCluster_ID, GeneID, Coordinate, Species, Membership_score]
#filter only for species pair of interest
exon_pairs_sub_df = exon_pairs_df.loc[exon_pairs_df.Species1.isin([species1, species2])]
exon_pairs_sub_df = exon_pairs_sub_df.loc[exon_pairs_df.Species2.isin([species1, species2])] 

#Get geneID-reclusteredID dictionary.
geneID_reclusteredID_dict = pd.Series(reclustered_genes_df.ClusterID.values, index=reclustered_genes_df.GeneID).to_dict()
#translate geneID with relative reclusteredID.
exon_pairs_sub_df["ClusterID1"] = exon_pairs_sub_df["GeneID1"].map(geneID_reclusteredID_dict)
exon_pairs_sub_df["ClusterID2"] = exon_pairs_sub_df["GeneID2"].map(geneID_reclusteredID_dict)
#filter only for exon hits within genes in the gene subcluster.
filtered_pairs_df = exon_pairs_sub_df.loc[exon_pairs_sub_df.ClusterID1 == exon_pairs_sub_df.ClusterID2]

#get a ExonID-ExClusterID dictionary.
exonID_exclusterID_dict = pd.Series(exon_clusters_df.ExCluster_ID.values, index=exon_clusters_df.Coordinate).to_dict()
#add ExCluster ID (just one, it MUST be identical for the two entries).
filtered_pairs_df["ExClusterID"] = filtered_pairs_df["ExonID1"].map(exonID_exclusterID_dict)
#generate new exon reclustering ID.
ExReclusteringID = [str(element)+re.sub(".*\.", ".", str(element1)) for (element, element1) in list(zip(list(filtered_pairs_df["ClusterID1"]), list(filtered_pairs_df["ExClusterID"])))]
filtered_pairs_df["ExReclusterID"] = ExReclusteringID

#get long format final_df
final_df_species1 = filtered_pairs_df[["ExReclusterID", "GeneID1", "ExonID1", "Species1"]]
final_df_species2 = filtered_pairs_df[["ExReclusterID", "GeneID2", "ExonID2", "Species2"]]
#rename columns
final_df_species1 = final_df_species1.rename(columns={"ExReclusterID" : "ExCluster_ID", "GeneID1" : "GeneID", "ExonID1" : "Coordinate", "Species1" : "Species"})
final_df_species2 = final_df_species2.rename(columns={"ExReclusterID" : "ExCluster_ID", "GeneID2" : "GeneID", "ExonID2" : "Coordinate", "Species2" : "Species"})
#join dataframes
final_df = pd.concat([final_df_species1, final_df_species2])
#get only unique entries
final_df = final_df.drop_duplicates()
#save final dataframe to output
final_df.to_csv(my_output, sep="\t", index=False, header=True, na_rep="NA")
