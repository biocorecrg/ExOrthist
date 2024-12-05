#!/usr/bin/env python
import argparse
import pandas as pd
import re  

#Read arguments
parser = argparse.ArgumentParser(description="Subset input files to only genes in a given gene cluster")
parser.add_argument("--annotation", "-a", required=True)
parser.add_argument("--overlap", "-o", required=True)
parser.add_argument("--ref_prots", "-r", required=True)
parser.add_argument("--gene_clusters", "-c", required=True)
parser.add_argument("--gene_clusterID", "-g", required=True)
parser.add_argument("--species", "-s", required=True)

args = parser.parse_args()
my_annot_file = args.annotation
my_overlap_file = args.overlap
my_ref_prots_file = args.ref_prots
my_gene_clusters = args.gene_clusters
my_gene_clusterID = args.gene_clusterID
my_species = args.species

#### Main
my_cluster_df = pd.read_csv(my_gene_clusters, sep="\t", header=None, index_col=False)
selected_genes = list(my_cluster_df[my_cluster_df.iloc[:,0]==str(my_gene_clusterID)].iloc[:,2])
#Subset GTF
my_gtf_df = pd.read_table(my_annot_file, header=None, index_col=False)
my_raw_data = [element for element in list(my_gtf_df.iloc[:,my_gtf_df.shape[1]-1]) if "protein_id" in element and "exon_number" in element]

my_raw_gene_id = [part for element in list(my_gtf_df.iloc[:,my_gtf_df.shape[1]-1]) for part in element.split(";") if "gene_id" in part]
my_gene_id = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_gene_id]
my_gtf_df["GeneID"] = my_gene_id
my_gtf_sub_df = my_gtf_df.loc[my_gtf_df.GeneID.isin(selected_genes)]
my_gtf_sub_df = my_gtf_sub_df.drop(columns=["GeneID"])
my_gtf_sub_df.to_csv(my_species+"_subsetted_annot.gtf", header=False, index=False, sep="\t", na_rep="NA")

#Subset Overlap file
my_overlap_df = pd.read_table(my_overlap_file, header=None, index_col=False, sep="\t", names=["OverlapID", "GeneID", "ExCoords", "Strand"])
#added after modification of main.nf output, 17/11/2020
my_overlap_df = my_overlap_df.drop(columns=["Strand"])
my_overlap_sub_df = my_overlap_df.loc[my_overlap_df["GeneID"].isin(selected_genes)]
my_overlap_sub_df.to_csv(my_species+"_subsetted_overlap_info.txt", header=False, index=False, sep="\t", na_rep="NA")

#Get exint_info (format: ProteinID|GeneID (ref protein), ex num.
my_refprot_df = pd.read_table(my_ref_prots_file, header=None, index_col=False, sep="\t", names=["GeneID", "RefProt"])
my_refprot_sub_df = my_refprot_df.loc[my_refprot_df.GeneID.isin(selected_genes)]
my_refprot_sub_df = my_refprot_sub_df.drop(columns=["GeneID"])
my_refprot_sub_df.to_csv(my_species+"_subsetted_ref_proteins.txt", sep="\t", header=True, index=False, na_rep="NA")
