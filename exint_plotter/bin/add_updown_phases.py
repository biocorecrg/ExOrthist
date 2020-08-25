#!/usr/bin/env python

import pandas as pd
import argparse  

#read arguments
parser = argparse.ArgumentParser(description="For each exon, add upstream and downstream phase")
parser.add_argument("--input_file", "-i", required=True)
parser.add_argument("--output_file", "-out", required=True)

args = parser.parse_args()
my_input = args.input_file
my_output_file = args.output_file

my_output = open(str(my_output_file), "w")

##### start processing
my_df = pd.read_table(str(my_input), header=None, names=["GeneID", "ExonID", "Phase", "Strand"])
exon_phases_dict = pd.Series(my_df.Phase.values, index=my_df.ExonID).to_dict() #key=exon, value=phase
exon_strand_dict = pd.Series(my_df.Strand.values, index=my_df.ExonID).to_dict() #key=exon, value=strand
exon_gene_dict = pd.Series(my_df.GeneID.values, index=my_df.ExonID).to_dict() #key=exon, value=geneID

gene_exons_dict = {} #initialize empty dictionary {GeneID : [list of ordered exons]}
my_grouped_df = my_df.groupby("GeneID") #for each gene
for name, group in my_grouped_df:
  gene_exons_dict[name] = sorted(list(group["ExonID"])) #dict value = sorted list of exons for each gene (key)

for my_exon in list(my_df.loc[:,"ExonID"]): #for each exon
  my_gene = exon_gene_dict[my_exon]     #select the gene
  my_exon_index = gene_exons_dict[my_gene].index(my_exon) #select the position of the exon in the gene
  my_strand = exon_strand_dict[my_exon] #select the strand

  if my_strand == "-": #in case of genes on the negative strand
    if my_exon_index > 0:
      my_exon_previous_ex = gene_exons_dict[my_gene][my_exon_index-1] #isolate the preceding exon
      up_phase = exon_phases_dict[my_exon_previous_ex] #isolate the phase of the preceding exon (phase of the upstream intron of query exon)
      down_phase = exon_phases_dict[my_exon]
    else:
      up_phase = "last"
      down_phase = exon_phases_dict[my_exon]

  elif my_strand == "+": #in case of genes on the positive strand
    if len(gene_exons_dict[my_gene]) > my_exon_index+1:
      my_exon_following_ex = gene_exons_dict[my_gene][my_exon_index+1] #isolate the following exon 
      up_phase = exon_phases_dict[my_exon] #isolate query exon phase (phase of the upstream intron)
      down_phase = exon_phases_dict[my_exon_following_ex] #isolate the phase of the following exon (phase of the downstream intron of query exon)
    else:
      up_phase = exon_phases_dict[my_exon]
      down_phase = "last"

  my_output.write(my_gene+"\t"+my_exon+"\t"+my_strand+"\t"+str(up_phase)+"\t"+str(down_phase)+"\n")
