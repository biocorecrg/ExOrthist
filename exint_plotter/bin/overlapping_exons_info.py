#!/usr/bin/env python

import argparse
import pandas as pd
import re
import collections

#Read arguments
parser = argparse.ArgumentParser(description="Subset input files to only genes in a given gene cluster")
parser.add_argument("--annotation", "-a", required=True)
parser.add_argument("--overlap", "-o", required=True)
parser.add_argument("--output_file", "-out", required=True)

args = parser.parse_args()
my_annot_file = args.annotation
my_overlap_file = args.overlap
my_output_file = args.output_file

##### Main
my_overlap_df = pd.read_table(my_overlap_file, sep="\t", header=None, names=["OverlapID", "GeneID", "ExCoords"])
my_gtf = pd.read_table(my_annot_file, sep="\t", header=None)
my_gtf_subset = my_gtf.iloc[:,my_gtf.shape[1]-1]

my_raw_gene_id = [part for element in list(my_gtf_subset) for part in element.split(";") if "gene_id" in part]
my_gtf["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_gene_id]

#Remove genes with exons annotated on different strands
geneID_strand_df = my_gtf.iloc[:,[6,my_gtf.shape[1]-1]].drop_duplicates() #If a gene has exons annotated on both strands, the geneID will be duplicated. 
selected_geneIDs = [item for item, count in collections.Counter(list(geneID_strand_df["geneID"])).items() if count == 1]
my_gtf = my_gtf.loc[my_gtf["geneID"].isin(selected_geneIDs)]
my_gtf["Coords"] = [str(element)+"-"+str(element1) for element, element1 in zip(list(my_gtf.iloc[:,3]), list(my_gtf.iloc[:,4]))]
my_gtf["CompleteCoords"] = [str(element)+":"+str(element1)+"-"+str(element2) for element, element1, element2 in zip(list(my_gtf.iloc[:,0]), list(my_gtf.iloc[:,3]), list(my_gtf.iloc[:,4]))]

#Add frequency and exon length
my_gtf_exons = my_gtf[my_gtf.iloc[:,2]=="CDS"] #Create a dictionary with key=coords, value=freq
my_exon_freq_dict = {key : value for key, value in collections.Counter(list(my_gtf_exons["Coords"])).items()}
my_overlap_df["Freq"] = my_overlap_df["ExCoords"].map(my_exon_freq_dict).fillna(0) #add frequency
my_overlap_df["Length"] = [int(re.sub(".*-", "",element))-int(re.sub("-.*", "", element)) for element in list(my_overlap_df["ExCoords"])] #add exon lenght
#Add complete coords
complete_coords_dict = pd.Series(my_gtf.CompleteCoords.values, index=my_gtf.Coords).to_dict()
my_overlap_df["CompleteCoords"] = my_overlap_df["ExCoords"].map(complete_coords_dict)
#Put a filter on the Freq: I think for now it is necessary because we don't have the exons from the FakeTranscripts (thus, there are exons from the clusters which have frequency 0).
my_overlap_df = my_overlap_df.loc[my_overlap_df.Freq > 0]
#Order df
my_overlap_df = my_overlap_df[["GeneID","OverlapID","CompleteCoords","Freq","Length"]]
#Write to file
my_overlap_df.to_csv(my_output_file, sep="\t", header=False, index=False, na_rep="NA")
