#!/usr/bin/env python

import argparse
import pandas as pd
import re

#Read arguments
parser = argparse.ArgumentParser(description="Subset input files to only genes in a given gene cluster")
parser.add_argument("--annotation", "-a", required=True)
parser.add_argument("--ref_prots", "-r", required=True)
parser.add_argument("--output", "-out", required=True)

args = parser.parse_args()
my_annot_file = args.annotation
my_ref_prots_file = args.ref_prots
my_output_file = args.output

my_gtf = pd.read_table(my_annot_file, sep="\t", header=None)
my_ex_int_num_df = pd.read_table(my_ref_prots_file, sep="\t", header=0)

#select only protein coding
my_gtf = my_gtf[my_gtf.iloc[:,my_gtf.shape[1]-1].str.contains("protein_id")]
my_gtf_subset = my_gtf.iloc[:,my_gtf.shape[1]-1]
my_gtf_subset = my_gtf_subset[my_gtf_subset.str.contains("protein_id")]
my_raw_prot_id = [part for element in list(my_gtf_subset) for part in element.split(";") if "protein_id" in part]
my_gtf["proteinID"] = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_prot_id]

ref_proteins_list = [re.sub("\\|.*", "", element) for element in list(my_ex_int_num_df["RefProt"])]
my_int_gtf = my_gtf.loc[my_gtf["proteinID"].isin(ref_proteins_list)] 
my_coords = [str(element)+":"+str(element1)+"-"+str(element2) for element, element1, element2 in zip(list(my_int_gtf.iloc[:,0]), list(my_int_gtf.iloc[:,3]), list(my_int_gtf.iloc[:,4]))]
my_final_gtf = pd.concat([pd.Series(my_coords), pd.Series(list(my_int_gtf.iloc[:,7]))], axis=1) #get a dataframe with exonID, exonPhase
my_final_gtf.to_csv(my_output_file, sep="\t", header=False, index=False, na_rep="NA")
