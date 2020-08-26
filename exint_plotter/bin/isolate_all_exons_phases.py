#!/usr/bin/env python

import argparse
import pandas as pd
import re
import collections

#Read arguments
parser = argparse.ArgumentParser(description="Subset input files to only genes in a given gene cluster")
parser.add_argument("--annotation", "-a", required=True)
parser.add_argument("--ref_phases", "-r", required=True)
parser.add_argument("--output", "-out", required=True)

args = parser.parse_args()
my_annot_file = args.annotation
my_ref_phases_file = args.ref_phases
my_output_file = args.output

#Main
my_gtf = pd.read_table(my_annot_file, sep="\t", header=None)
my_phases = pd.read_table(my_ref_phases_file, sep="\t", header=None, names=["Coords", "Phase"])

my_int_gtf = my_gtf.loc[my_gtf.iloc[:,2]=="CDS"]
my_coords = [str(element)+":"+str(element1)+"-"+str(element2) for element, element1, element2 in zip(list(my_int_gtf.iloc[:,0]), list(my_int_gtf.iloc[:,3]), list(my_int_gtf.iloc[:,4]))]
my_int_gtf["Coords"] = my_coords
my_int_gtf = my_int_gtf.rename(columns={7:"Phase"})
my_filt_gtf = my_int_gtf.loc[:,["Coords", "Phase"]].drop_duplicates()
my_unique_coords = [key for key, value in collections.Counter(list(my_filt_gtf["Coords"])).items() if value == 1] 
my_duplicated_coords = [key for key, value in collections.Counter(list(my_filt_gtf["Coords"])).items() if value > 1] 
my_duplicated_phases = my_phases.loc[my_phases["Coords"].isin(my_duplicated_coords)]
my_unique_phases = my_filt_gtf.loc[my_filt_gtf["Coords"].isin(my_unique_coords)]
my_final_phases = pd.concat([my_duplicated_phases, my_unique_phases]).sort_values(by=["Coords"])
my_final_phases.to_csv(my_output_file, sep="\t", header=False, index=False, na_rep="NA")
