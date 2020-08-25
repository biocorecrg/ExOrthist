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
my_overlap_df = pd.read_table(my_overlap_file, header=None, index_col=False, sep="\t")
my_overlap_sub_df = my_overlap_df[my_overlap_df.iloc[:,1].isin(selected_genes)]
my_overlap_sub_df.to_csv(my_species+"_subsetted_overlap_info.txt", header=False, index=False, sep="\t", na_rep="NA")

#Get exint_info (format: ProteinID|GeneID (ref protein), ex num.
my_refprot_df = pd.read_table(my_ref_prots_file, header=None, index_col=False, sep="\t", names=["GeneID", "RefProt"])
my_refprot_sub_df = my_refprot_df.loc[my_refprot_df.GeneID.isin(selected_genes)]
my_refprot_sub_df = my_refprot_sub_df.drop(columns=["GeneID"])
my_refprot_sub_df.to_csv(my_species+"_subsetted_ref_proteins.txt", sep="\t", header=True, index=False, na_rep="NA")

#selected_proteins = list(my_refprot_sub_df["ProteinID"])
#create a dictionary ProteinID - exon_num
#my_protein_id = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in my_raw_data for part in element.split(";") if "protein_id" in part]
#my_ex_num = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in my_raw_data for part in element.split(";") if "exon_number" in part]
#protein_id_ex_num_df = pd.DataFrame({"ProteinID" : my_protein_id, "ExNum" : [int(element) for element in my_ex_num]})
#protein_id_ex_num_df = protein_id_ex_num_df.loc[protein_id_ex_num_df.ProteinID.isin(selected_proteins)] #filter only for the ref proteins relative to this cluster file
#protein_id_ex_num_df_grouped = protein_id_ex_num_df.groupby("ProteinID")
#
#my_exintnum_df = pd.DataFrame()
#for name, group in protein_id_ex_num_df_grouped:
#  my_entry = group.loc[group.ExNum==max(list(group["ExNum"]))] 
#  my_exintnum_df = pd.concat([my_exintnum_df, my_entry])
#
#my_protID_exnum_dict = pd.Series(my_exintnum_df.ExNum.values, index=my_exintnum_df.ProteinID).to_dict()
#my_refprot_sub_df["ExNum"] = my_refprot_sub_df["ProteinID"].map(my_protID_exnum_dict)
#my_refprot_sub_df["RefProt"] = [element+"|"+element1 for element, element1 in zip(list(my_refprot_sub_df["ProteinID"]), list(my_refprot_sub_df["GeneID"]))]
#my_exint_final_df = my_refprot_sub_df.drop(columns=["GeneID", "ProteinID"]) #drop extra columns
#my_exint_final_df = my_exint_final_df[["RefProt", "ExNum"]] #reorder columns
#my_exint_final_df.to_csv(my_species+"_subsetted_ref_exon_int_number.txt", sep="\t", header=True, index=False, na_rep="NA")
