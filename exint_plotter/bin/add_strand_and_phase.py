#!/usr/bin/env python

import argparse
import pandas as pd
import re
import collections

#Read arguments
parser = argparse.ArgumentParser(description="Subset input files to only genes in a given gene cluster")
parser.add_argument("--annotation", "-a", required=True)
parser.add_argument("--phases", "-p", required=True)
parser.add_argument("--exons", "-e", required=True)
parser.add_argument("--output", "-out", required=True)

args = parser.parse_args()
my_annot_file = args.annotation
my_phases_file = args.phases
my_exons_file = args.exons
my_output_file = args.output

my_selected_exons = pd.read_table(my_exons_file, sep="\t", header=0)
my_gtf = pd.read_table(my_annot_file, sep="\t", header=None)
my_exon_phases_df = pd.read_table(my_phases_file, sep="\t", header=None, names=["Coords", "Phase"]) 

my_gtf_subset = my_gtf.iloc[:,my_gtf.shape[1]-1] #select the last columns of the GTF
my_raw_gene_id = [part for element in list(my_gtf_subset) for part in element.split(";") if "gene_id" in part] #select geneID
my_gtf["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_gene_id]

#Remove genes with exons annotated on different strands
geneID_strand_df = my_gtf.iloc[:,[6,my_gtf.shape[1]-1]].drop_duplicates() #If a gene has exons annotated on both strands, the geneID will be duplicated. 
selected_geneIDs = [item for item, count in collections.Counter(list(geneID_strand_df["geneID"])).items() if count == 1]
my_gtf = my_gtf.loc[my_gtf["geneID"].isin(selected_geneIDs)]
my_gtf["Coords"] = [str(element)+":"+str(element1)+"-"+str(element2) for element, element1, element2 in zip(list(my_gtf.iloc[:,0]), list(my_gtf.iloc[:,3]), list(my_gtf.iloc[:,4]))]
my_gtf = my_gtf.rename(columns={6 : "Strand"})
#my_gtf_filt = my_gtf.loc[:,["Coords","CompleteCoords","Strand"]].drop_duplicates() #select only coords and strand
my_gtf_filt = my_gtf.loc[:,["Coords","Strand"]].drop_duplicates() #select only coords and strand

#Create a dictionary with key=Coords, value=strand
my_coords_strand_dict = pd.Series(my_gtf_filt.Strand.values, index=my_gtf_filt.Coords).to_dict()
my_selected_exons["Strand"] = my_selected_exons["ExonID"].map(my_coords_strand_dict)
#Create a dictionary with key=Coords, value=phase
my_coords_phase_dict = pd.Series(my_exon_phases_df.Phase.values, index=my_exon_phases_df.Coords).to_dict()
my_selected_exons["Phase"] = my_selected_exons["ExonID"].map(my_coords_phase_dict).fillna(0)
my_selected_exons["Phase"] = [int(element) for element in list(my_selected_exons["Phase"])] #transform to integer
#Create a dictionary with key=Coords, value=CompleteCoords
my_final_df = my_selected_exons[["geneID", "ExonID", "Phase", "Strand"]]
my_final_df.to_csv(my_output_file, sep="\t", header=False, index=False, na_rep="NA")
