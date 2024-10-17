#!/usr/bin/env python

import argparse
import pandas as pd
import re
import sys
import collections

#Read arguments
parser = argparse.ArgumentParser(description="Generate input for exint plotter")
parser.add_argument("--annotation", "-a", required=True)
parser.add_argument("--overlap", "-o", required=True)
parser.add_argument("--gene_clusters", "-c", required=True)
parser.add_argument("--ex_clusters", "-e", required=True)
parser.add_argument("--ref_prots", "-r", required=True)
parser.add_argument("--species", "-s", required=True)
parser.add_argument("--output_file", "-out", required=True)

args = parser.parse_args()
my_annot_file = args.annotation
my_overlap_file = args.overlap
my_exon_clusters = args.ex_clusters
my_gene_clusters_file = args.gene_clusters
my_ref_prots_file = args.ref_prots
my_species = args.species
my_output_file = args.output_file

################ OVERLAPPING EXONS INFO #####################
my_overlap_df = pd.read_table(my_overlap_file, sep="\t", header=None, names=["ExOverlapID", "GeneID", "ExCoords","Strand"])
my_gtf = pd.read_table(my_annot_file, sep="\t", header=None)
#Exit if gtf does not have the expected number of fields.
if my_gtf.shape[1] != 9:
  sys.exit("GTF does not have the expected number of fields")
#rename GTF entries.
my_gtf = my_gtf.rename(columns={0:"Chr", 1:"Source", 2:"Feature", 3:"Start", 4:"Stop", 5:"Score", 6:"Strand", 7:"Phase", 8:"Info"})
#subset GTF to only entries of CDS exons (which automatically have the ProteinID)
#my_gtf = my_gtf.loc[my_gtf["Info"].str.contains("protein_id")]
my_gtf = my_gtf.loc[my_gtf["Feature"]=="CDS"]
#add extra into (GeneID, ProteinID, ExNum) as GTF Separate columns
my_gtf_subset = my_gtf["Info"]
my_raw_gene_id = [part for element in list(my_gtf_subset) for part in element.split(";") if "gene_id" in part]
my_gtf["GeneID"] = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_gene_id]
# Select the first subfield containing the protein ID. Useful in case of weird GTF structure
#protein_id_subfield = list(my_gtf_subset)[0].split(";").index([element for element in list(my_gtf_subset)[0].split(";") if "protein_id" in element][0])
#my_raw_prot_id = [element.split(";")[protein_id_subfield] for element in list(my_gtf_subset)]
my_raw_prot_id = [part for element in list(my_gtf_subset) for part in element.split(";") if "protein_id" in part]
my_gtf["ProteinID"] = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_prot_id]
#The transcriptID will be used to derive the annotation status
my_raw_transcriptID = [part for element in list(my_gtf_subset) for part in element.split(";") if "transcript_id" in part] #select transcriptID
my_gtf["TranscriptID"] = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_transcriptID]
#The exon number is useful to define the relative position.
my_raw_exon_num = [part for element in list(my_gtf_subset) for part in element.split(";") if "exon_number" in part]
my_exon_num = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_exon_num]
my_gtf["ExNum"] = my_exon_num

#Remove genes with exons annotated on different strands (very weird cases)
geneID_strand_df = my_gtf.loc[:,["Strand","GeneID"]].drop_duplicates() #If a gene has exons annotated on both strands, the geneID will be duplicated.
selected_geneIDs = [item for item, count in collections.Counter(list(geneID_strand_df["GeneID"])).items() if count == 1]
my_gtf = my_gtf.loc[my_gtf["GeneID"].isin(selected_geneIDs)]
#Add coordinates to GTF
#my_gtf["Coords"] = [str(element)+"-"+str(element1) for element, element1 in zip(list(my_gtf["Start"]), list(my_gtf["Stop"]))]
my_gtf["Coords"] = [str(element)+":"+str(element1)+"-"+str(element2) for element, element1, element2 in zip(list(my_gtf["Chr"]), list(my_gtf["Start"]), list(my_gtf["Stop"]))]
#Add chr to start-stop coords in the overlapping group GTF
geneID_chr_dict = pd.Series(my_gtf.Chr.values, index=my_gtf.GeneID).to_dict() #the duplicated keys are automatically overwritten
my_overlap_df["ExonID"] = [str(element)+":"+str(element1) for element, element1 in zip(list(my_overlap_df["GeneID"].map(geneID_chr_dict)), list(my_overlap_df["ExCoords"]))]
#Add frequency and exon length
my_gtf_exons = my_gtf.loc[my_gtf.Feature=="CDS"]
my_exon_freq_dict = {key : value for key, value in collections.Counter(list(my_gtf_exons["Coords"])).items()} #Create a dictionary with key=coords, value=freq
my_overlap_df["Freq"] = my_overlap_df["ExonID"].map(my_exon_freq_dict).fillna(0) #add frequency
my_overlap_df["Length"] = [int(re.sub(".*-", "",element))-int(re.sub(".*:", "", re.sub("-.*", "", element))) for element in list(my_overlap_df["ExonID"])] #add exon lenght

#Put a filter on the Freq: I think for now it is necessary because we don't have the exons from the FakeTranscripts (thus, there are exons from the clusters which have frequency 0).
my_overlap_df = my_overlap_df.loc[my_overlap_df.Freq > 0]
my_overlap_df = my_overlap_df[["GeneID", "ExOverlapID", "ExonID", "Freq", "Length"]] #Order df

################## SELECT OVERLAPPING EXONS #####################
my_overlap_df = my_overlap_df.fillna(0)
#select exons from exon clusters
#header: ExCluster_ID, GeneID, Coordinate, Species, Membership_score
exon_clusters_df = pd.read_table(str(my_exon_clusters), header=0, sep="\t")
exons_in_clusters = [re.sub(":-", "", re.sub(":\+", "", element)) for element in list(exon_clusters_df["Coordinate"])]

my_selected_overlap_df = pd.DataFrame(columns=["GeneID", "ExOverlapID", "ExonID", "Freq", "Length"])
#group by overlap ID
my_grouped_overlap_df = my_overlap_df.groupby("ExOverlapID")
for name, group in my_grouped_overlap_df:
  all_exs = list(group["ExonID"])
  my_ex = [ex for ex in all_exs if ex in exons_in_clusters] #select the exon in the exon clusters for each overalpping group (there should be only one).
  if len(my_ex) == 1:
    selected_elements_df = group.loc[group.ExonID==my_ex[0]]
  else: #if none of the exons in the overlapping group make it to the exon clusters, select the most frequent form.
    all_freq_list = list(group["Freq"])
    max_freq = max(all_freq_list)
    selected_elements_df = group.loc[group.Freq==max_freq]                  
  if selected_elements_df.shape[0] > 1: #if there are some forms with equal frequency, select the longest.
    selected_elements_df = selected_elements_df.loc[selected_elements_df.Length==max(list(selected_elements_df["Length"]))]
  #header: ["GeneID", "ExOverlapID", "ExonID", "Freq", "Length"]
  my_selected_overlap_df = my_selected_overlap_df.append(selected_elements_df, ignore_index=True) #add the selected element to the final dataframe.

#Print out Coords - Overlapping chosen coords file. This will be used to translate the scores from the best-hits
my_overlapID_chosenID_df = my_overlap_df.loc[:,["ExonID", "ExOverlapID"]]
#Create an ExOverlapID - ChosenID dictionary
overlapID_chosenID_dict = pd.Series(my_selected_overlap_df.ExonID.values, index=my_selected_overlap_df.ExOverlapID).to_dict()
my_overlapID_chosenID_df["ExOverlapID"] = my_overlapID_chosenID_df["ExOverlapID"].map(overlapID_chosenID_dict)
my_overlapID_chosenID_df = my_overlapID_chosenID_df.rename(columns={"ExOverlapID" : "ChosenID"})
my_overlapID_chosenID_df.to_csv(my_species+"_overlapID_chosenID.txt", sep="\t", header=True, index=False, na_rep="NA") #save to file
 
################## ISOLATE REF PROTEIN EXONS PHASES #####################
my_ex_int_num_df = pd.read_table(my_ref_prots_file, sep="\t", header=None, names=["GeneID", "RefProt"])
ref_proteins_list = list(my_ex_int_num_df["RefProt"])
#ref_proteins_list = [re.sub("\\|.*", "", element) for element in list(my_ex_int_num_df["RefProt"])]
my_ref_gtf = my_gtf.loc[my_gtf["ProteinID"].isin(ref_proteins_list)] 
my_ref_phases_df = pd.concat([my_ref_gtf["Coords"], pd.Series(list(my_ref_gtf.iloc[:,7]))], axis=1) #get a dataframe with exonID, RefExonPhase

################## ISOLATE ALL EXONS PHASES #####################
my_all_phases_df = my_gtf.loc[:,["Coords", "Phase"]].drop_duplicates()
my_unique_coords = [key for key, value in collections.Counter(list(my_all_phases_df["Coords"])).items() if value == 1] #exons in the same phase across all isoforms 
my_duplicated_coords = [key for key, value in collections.Counter(list(my_all_phases_df["Coords"])).items() if value > 1] #exons in different phases across isoforms
my_duplicated_phases = my_ref_phases_df.loc[my_ref_phases_df["Coords"].isin(my_duplicated_coords)] #select the reference phase for the exons annotated with differnet phases.
my_unique_phases = my_all_phases_df.loc[my_all_phases_df["Coords"].isin(my_unique_coords)]
my_final_phases = pd.concat([my_duplicated_phases, my_unique_phases]).sort_values(by=["Coords"])
my_final_phases.to_csv(my_output_file, sep="\t", header=False, index=False, na_rep="NA")

################## ADD STRAND AND PHASES ################
my_strand_df = my_gtf.loc[:,["Coords","Strand"]].drop_duplicates() #select only coords and strand

#Create a dictionary with key=Coords, value=strand
my_coords_strand_dict = pd.Series(my_strand_df.Strand.values, index=my_strand_df.Coords).to_dict()
my_selected_overlap_df["Strand"] = my_selected_overlap_df["ExonID"].map(my_coords_strand_dict)
#Create a dictionary with key=Coords, value=phase
my_coords_phase_dict = pd.Series(my_final_phases.Phase.values, index=my_final_phases.Coords).to_dict()
my_selected_overlap_df["Phase"] = my_selected_overlap_df["ExonID"].map(my_coords_phase_dict).fillna(0)
my_selected_overlap_df["Phase"] = [int(element) for element in list(my_selected_overlap_df["Phase"])] #transform to integer
#Create a dictionary with key=Coords, value=Coords
my_phase_strand_df = my_selected_overlap_df[["GeneID", "ExonID", "Phase", "Strand"]]

################## ADD UPDOWN PHASES ###################
exon_phases_dict = pd.Series(my_phase_strand_df.Phase.values, index=my_phase_strand_df.ExonID).to_dict() #key=exon, value=phase
exon_strand_dict = pd.Series(my_phase_strand_df.Strand.values, index=my_phase_strand_df.ExonID).to_dict() #key=exon, value=strand
exon_gene_dict = pd.Series(my_phase_strand_df.GeneID.values, index=my_phase_strand_df.ExonID).to_dict() #key=exon, value=geneID

gene_exons_dict = {} #initialize empty dictionary {GeneID : [list of ordered exons]}
UpDown_phases_df = pd.DataFrame(columns=["GeneID", "ExonID", "Strand", "UpPhase", "DownPhase"])
my_grouped_df = my_phase_strand_df.groupby("GeneID") #for each gene
for name, group in my_grouped_df:
  gene_exons_dict[name] = sorted(list(group["ExonID"])) #dict value = sorted list of exons for each gene (key)

for my_exon in list(my_phase_strand_df.loc[:,"ExonID"]): #for each exon
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
  element_df = pd.DataFrame([[my_gene, my_exon, my_strand, str(up_phase), str(down_phase)]], columns=["GeneID", "ExonID", "Strand", "UpPhase", "DownPhase"])
  UpDown_phases_df = pd.concat([UpDown_phases_df, element_df])
  #my_output.write(my_gene+"\t"+my_exon+"\t"+my_strand+"\t"+str(up_phase)+"\t"+str(down_phase)+"\n")

################## ADD ANNOTATION STATUS ####################
#The annotation status is derived from the name of the transcripID
my_coords_transcript_df = my_gtf.loc[:,["Coords", "TranscriptID"]]

#select all the exons in fake transcripts
fake_transcript_exons = list(set(list(my_coords_transcript_df.loc[my_coords_transcript_df.TranscriptID.str.contains("fB")]["Coords"])))
#select all the exons in annotated transcripts
annotated_exons = list(set(list(my_coords_transcript_df.loc[my_coords_transcript_df.TranscriptID.str.contains("fB")==False]["Coords"])))
#select the not-annotated exons (which are present exclusively in the fake transcripts)
not_annotated_exons = [element for element in fake_transcript_exons if element not in annotated_exons]
#create a dictionary with ExonID - annotation status
exonID_annot_status_dict = {}
for exon in annotated_exons:
  exonID_annot_status_dict[exon] = "annotated"
for exon in not_annotated_exons:
  exonID_annot_status_dict[exon] = "not_annotated"
#add annotation status to selected exons.
UpDown_phases_df["AnnotStatus"] = UpDown_phases_df["ExonID"].map(exonID_annot_status_dict)

################## ADD EXON CLUSTER INFO ###################
exon_clusters_df_species = exon_clusters_df.loc[exon_clusters_df["Species"] == my_species]
#create dictionary with key=exon_coords, value=exon_cluster_id
#it also considers only start or only stop coords
coords_cluster_dict = pd.Series(exon_clusters_df_species.ExCluster_ID.values, index=exon_clusters_df_species.Coordinate).to_dict()
for my_id in list(coords_cluster_dict.keys()):
  start_coord = re.sub('-.*', '', my_id)
  stop_coord = re.sub(':.*-', ':', my_id)
  coords_cluster_dict[start_coord] = coords_cluster_dict[my_id]
  coords_cluster_dict[stop_coord] = coords_cluster_dict[my_id]

#get exon_cluster_id
my_final_dict = {}
for element in list(UpDown_phases_df["ExonID"]): #for each exon_coord
  if element in list(coords_cluster_dict.keys()): #for each ex_coords in the exon_cluster       
    my_final_dict[element] = coords_cluster_dict[element] #save clusterID as value for ex_coord key
  else:
    element_start = re.sub('-.*', '', element) #select start and stop
    if element_start in list(coords_cluster_dict.keys()):
      my_final_dict[element] = coords_cluster_dict[element_start]
    else:
      element_stop = re.sub(':.*-', ':', element)
      if element_stop in list(coords_cluster_dict.keys()):
        my_final_dict[element] = coords_cluster_dict[element_stop]
      else:
        my_final_dict[element] = element

UpDown_phases_df["EX_clusterID"] = [my_final_dict[my_id] for my_id in list(UpDown_phases_df["ExonID"])]
UpDown_phases_df["Cluster_status"] = ["not_in_clusters"  if x not in list(coords_cluster_dict.values()) else "Conserved" for x in list(UpDown_phases_df["EX_clusterID"])]

################## ADD FAKE COORDS ##################
exon_info_df = UpDown_phases_df
my_grouped_df =  exon_info_df.groupby("GeneID") #group exons by gene

#define the functions
def reorder_exons(group):
  if list(group["Strand"])[0] == "+": #check first element of strand list. They should be all the same anyways.
    my_group = group.sort_values(by=["ExonID"]) #order in ascending order
  elif list(group["Strand"])[0] == "-":
    my_group = group.sort_values(by=["ExonID"], ascending=False)
  return(my_group)

def add_fake_coords(group):
  #group["FakeStart"] = list(range(0,15*group.shape[0],15))
  #group["FakeStop"] = list(range(15,15*(group.shape[0]+1),15))
  group["FakeStart"] = list(range(0,30*group.shape[0],30))
  group["FakeStop"] = list(range(15,30*(group.shape[0]),30))
  return(group)

#I  don't know why but the apply returns a pandas df, not a GroupBy object.
#I need the drop because it was somehow creating a multi level index in the output df.
my_ordered_grouped_df =  my_grouped_df.apply(reorder_exons).reset_index(drop=True)
my_ordered_grouped_df["State"] =  pd.Series(["Exon"]*my_ordered_grouped_df.shape[0])

#adding the fake coordinates here already
my_grouped_df = my_ordered_grouped_df.groupby("GeneID")
my_FakeCoords_df = my_grouped_df.apply(add_fake_coords)

################## ISOFORM INFO ###############
#Exon number in isoform (exon index) and by isoform
ExNum_in_isoform_df = my_gtf[["Coords", "ExNum", "GeneID", "ProteinID"]].drop_duplicates() #remove duplicated rows (exons with different splices sites)
ExNum_in_isoform_df["ExNum"] = [int(element) for element in list(ExNum_in_isoform_df["ExNum"])] #change data type of ExNum to integer
grouped_df = ExNum_in_isoform_df.groupby("ProteinID")
ExNum_by_isoform_df = pd.DataFrame(columns=["GeneID", "ProteinID", "ExNum"])
for name, group in grouped_df:
  my_group = group.loc[group.ExNum == max(list(group["ExNum"])),["GeneID", "ProteinID", "ExNum"]].drop_duplicates()
  ExNum_by_isoform_df = pd.concat([ExNum_by_isoform_df, my_group])

#Save ExNum_in_isoform_df to file
ExNum_in_isoform_df.to_csv(my_species+"_ExNum_in_isoform", sep="\t", header=False, index=False, na_rep="NA")

################## EXON POSITION INFO #############
my_ex_index_df_grouped = ExNum_in_isoform_df.groupby("ProteinID")
my_first_ex_list = []
for name, group in my_ex_index_df_grouped:
  my_first_ex_list.append(list(group.loc[group.ExNum==min(list(group["ExNum"]))]["Coords"])[0]) #select the first element of the list

my_last_ex_dict = pd.Series(ExNum_by_isoform_df.ExNum.values, index=ExNum_by_isoform_df.ProteinID).to_dict()
ExNum_in_isoform_df["TotEx"] = ExNum_in_isoform_df["ProteinID"].map(my_last_ex_dict) #Add total number of exons for each entries
my_last_ex_list = list(ExNum_in_isoform_df.loc[ExNum_in_isoform_df["ExNum"]==ExNum_in_isoform_df["TotEx"]]["Coords"])
my_internal_ex_list = [element for element in list(ExNum_in_isoform_df["Coords"]) if element not in my_first_ex_list+my_last_ex_list]
my_first_last_ex_raw = list(set(my_first_ex_list))+list(set(my_last_ex_list))
my_first_last_ex_list = [element for element, count in collections.Counter(my_first_last_ex_raw).items() if count > 1]

#create dictionary for exon status. This works because update will replace the last value.
my_first_ex_dict = {element : "first" for element in my_first_ex_list}
my_last_ex_dict = {element : "last" for element in my_last_ex_list}
my_firstlast_ex_dict = {element : "first;last" for element in my_first_last_ex_list}
my_internal_ex_dict = {element : "Internal" for element in my_internal_ex_list}
my_all_ex_dict = {**my_first_ex_dict, **my_last_ex_dict, **my_firstlast_ex_dict, **my_internal_ex_dict}

#create dictionary for gene clusters key=geneID : value=clusterID
my_clusters = pd.read_table(my_gene_clusters_file, sep="\t", header=None)
my_gene_clusters_dict = pd.Series(list(my_clusters.iloc[:,0]), index=list(my_clusters.iloc[:,2])).to_dict()

#Join results in final table
my_FakeCoords_df["ExPosition"] = my_FakeCoords_df["ExonID"].map(my_all_ex_dict)
my_FakeCoords_df["ClusterID"] = my_FakeCoords_df["GeneID"].map(my_gene_clusters_dict)
my_FakeCoords_df["State"] = "Exon"
#Reorder columns
my_final_df = my_FakeCoords_df[["ClusterID", "GeneID", "ExonID", "Strand", "UpPhase", "DownPhase", "AnnotStatus", "Cluster_status", "EX_clusterID", "State", "ExPosition", "FakeStart", "FakeStop"]]
#write to file
my_final_df.to_csv(my_output_file, sep="\t", header=True, index=True, na_rep="NA") 
