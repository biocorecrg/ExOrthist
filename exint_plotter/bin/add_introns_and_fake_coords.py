#!/usr/bin/env python

import pandas as pd
import argparse  

#read arguments
parser = argparse.ArgumentParser(description="For each exon, add upstream and downstream phase")
parser.add_argument("--input_file", "-i", required=True)
parser.add_argument("--output_file", "-o", required=True)

args = parser.parse_args()
my_input = args.input_file
my_output = args.output_file

#main
my_df = pd.read_table(str(my_input), sep="\t")
my_grouped_df = my_df.groupby("GeneID") #group exons by gene

#define the functions
def reorder_exons(group):
  if list(group["Strand"])[0] == "+": #check first element of strand list. They should be all the same anyways.
    my_group = group.sort_values(by=["ExonID"]) #order in ascending order
  elif list(group["Strand"])[0] == "-":
    my_group = group.sort_values(by=["ExonID"], ascending=False)
  return(my_group)

def add_fake_coords(group):
  group["FakeStart"] = list(range(0,15*group.shape[0],15))
  group["FakeStop"] = list(range(15,15*(group.shape[0]+1),15))
  return(group)

#I  don't know why but the apply returns a pandas df, not a GroupBy object.
#I need the drop because it was somehow creating a multi level index in the output df.
my_ordered_grouped_df =  my_grouped_df.apply(reorder_exons).reset_index(drop=True)
my_ordered_grouped_df["State"] =  pd.Series(["Exon"]*my_ordered_grouped_df.shape[0])

my_introns_df = my_ordered_grouped_df.loc[(my_ordered_grouped_df.UpPhase!="last") & (my_ordered_grouped_df.DownPhase!="last")]
#NB: this shitty use of index.values is necessary because it was operating a translation
my_introns_df["State"] = pd.Series(["Intron"]*list(my_introns_df.index.values)[-1])
my_introns_df["ExonID"] = pd.Series(["Intron"+str(my_index) for my_index in list(range(0, list(my_introns_df.index.values)[-1]))])
my_introns_df["UpPhase"] = pd.Series(["NA"]*list(my_introns_df.index.values)[-1])
my_introns_df["DownPhase"] = pd.Series(["NA"]*list(my_introns_df.index.values)[-1])
my_introns_df["Cluster_status"] = pd.Series(["Conserved"]*list(my_introns_df.index.values)[-1])
my_introns_df["ClusterID"] = my_introns_df["ExonID"]

#colpo di genio
my_introns_df = my_introns_df.set_index(pd.Index([value+0.5 for value in list(my_introns_df.index.values)])) #get floating point indexes for the introns df
final_df = my_ordered_grouped_df.append(my_introns_df, ignore_index=False).sort_index() #join the dataframes and order them by index

#adding the fake coordinates here already
my_grouped_df = final_df.groupby("GeneID") #group exons by gene
my_FakeCoords_df = my_grouped_df.apply(add_fake_coords)

my_FakeCoords_df.to_csv(str(my_output), sep="\t", index=True, header=True) #write to file
