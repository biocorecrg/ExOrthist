#!/usr/bin/env python

import argparse
import pandas as pd
import re
import collections

#Read arguments
parser = argparse.ArgumentParser(description="Add annotation status for all the exons in the considered genes")
parser.add_argument("--annotation", "-a", required=True)
parser.add_argument("--exons", "-e", required=True)
parser.add_argument("--output", "-out", required=True)

args = parser.parse_args()
my_annot_file = args.annotation
my_exons_file = args.exons
my_output_file = args.output

my_selected_exons = pd.read_table(my_exons_file, sep="\t", header=None, names=["GeneID", "ExonID", "Strand", "UpPhase", "DownPhase"])
my_gtf = pd.read_table(my_annot_file, sep="\t", header=None)
my_gtf = my_gtf[my_gtf.iloc[:,my_gtf.shape[1]-1].str.contains("transcript_id")]
my_gtf["Coords"] = [str(element)+":"+str(element1)+"-"+str(element2) for element, element1, element2 in zip(list(my_gtf.iloc[:,0]), list(my_gtf.iloc[:,3]), list(my_gtf.iloc[:,4]))] #Add exons coordinates

my_gtf_subset = my_gtf.iloc[:,my_gtf.shape[1]-2] #select the column with transcript info, which is now second to last.
my_raw_transcriptID = [part for element in list(my_gtf_subset) for part in element.split(";") if "transcript_id" in part] #select transcriptID
my_gtf["transcriptID"] = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_transcriptID]
my_coords_transcript_df = my_gtf.loc[:,["Coords", "transcriptID"]]

#select all the exons in fake transcripts
fake_transcript_exons = list(set(list(my_coords_transcript_df.loc[my_coords_transcript_df.transcriptID.str.contains("fB")]["Coords"])))
#select all the exons in annotated transcripts
annotated_exons = list(set(list(my_coords_transcript_df.loc[my_coords_transcript_df.transcriptID.str.contains("fB")==False]["Coords"])))
#select the not-annotated exons (which are present exclusively in the fake transcripts)
not_annotated_exons = [element for element in fake_transcript_exons if element not in annotated_exons]
#create a dictionary with ExonID - annotation status
exonID_annot_status_dict = {}
for exon in annotated_exons:
  exonID_annot_status_dict[exon] = "annotated"
for exon in not_annotated_exons:
  exonID_annot_status_dict[exon] = "not_annotated"
#add annotation status to selected exons.
my_selected_exons["AnnotStatus"] = my_selected_exons["ExonID"].map(exonID_annot_status_dict)
#save to output file
my_selected_exons.to_csv(my_output_file, sep="\t", header=False, index=False, na_rep="NA")
