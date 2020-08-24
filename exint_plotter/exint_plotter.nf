#!/usr/bin/env nextflow

/*
 * Exint plotter for Exhortist pipeline
 *
 */

/*
================================
Fede's attempt at exint plotter
================================
 */

/*
 * Input parameters:
 */

log.info  """
Executing with the following parameters:

annotations (GTF files)			: ${params.annotations}
exon overlap info			: ${params.overlap}
exon introns nums			: ${params.exintnum}
exon clusters (pipeline output)		: ${params.exonclusters}
exon best scores (pipeline output)	: ${params.bestscores}
gene clusters (same as pipeline input)	: ${params.geneclusters}
geneID					: ${params.geneID}
"""

/*
 * Create channels for input data
 */


//This channel will contain a list of the GTF files, in theory each with a key
//The key corresponds to the value assumed by the wildcard in the annotation variable (which is defined in the params.config)
//annotations  = "$baseDir/data/GTF/*_annot.gtf"
//The set operator assign the channel to the variable specified as a closure parameter
Channel 
    .fromFilePairs(params.annotations, size: 1)
    .ifEmpty{error "Cannot find any annotation matching: ${params.annotations}"}
    .set{annotations}
//I need more channels because I use the same input in more processes
//Create channel for overlapping exon information
//The key is the species, same as for the annotations channel
Channel
    .fromFilePairs( params.overlap, size: 1)
    .ifEmpty{error "Cannot find any overlap info: ${params.overlap}"}
    .set{overlap_info}

//Create channel for files with info regarding ex/int numbers
Channel
    .fromFilePairs(params.exintnum, size: 1)
    .ifEmpty{error "Cannot find any overlap info: ${params.overlap}"}
    .set{exint_num_info_all}

//create a joint channel where each key is paired with the corresponding files
annotations.join(overlap_info).join(exint_num_info_all).into{all_input_info_raw}

my_geneID = "${params.geneID}"
gene_clusters = file(params.geneclusters)
process isolate_cluster_id {
	tag {"${geneID}"}
	input:
	val(my_geneID)
	file(gene_clusters)
	output:
	stdout into (gene_clusterID, gene_clusterID1, gene_clusterID2, gene_clusterID3)
	script:
	"""
	#!/usr/bin/env python

	import pandas as pd
	
	my_df = pd.read_csv("${gene_clusters}", sep="\t", header=None, index_col=False)
	clusterID = list(my_df[my_df.iloc[:,2]=="${my_geneID}"].iloc[:,0])[0]
	print(clusterID, end='')
	"""
}

/*
 * Filter only for the species actually having exons in the exon clusters.
 */
//The pipeline crushes otherwise. But there is a distinction: if a gene gets plotted with no exons, it means that its exons make it to the exon clusters. They simply do not have any ortholog exons to the query species.
my_exon_clusters = file(params.exonclusters)
process select_species_with_orthologs {
	tag {"${clusterID}"}
	input:
	file(my_exon_clusters)
	val(clusterID) from gene_clusterID3
	output:
	stdout into selected_species
	script:
	"""
	#!/usr/bin/env python
	import pandas as pd
	import re

	my_df = pd.read_table("${my_exon_clusters}", sep="\t", header=0, index_col=False)
	my_df["GeneClusterID"] = [re.sub("\\..*", "", element) for element in my_df["ClusterID"]]
	selected_species = list(set(list(my_df[my_df.GeneClusterID=="${clusterID}"]["Sps"])))
	print(selected_species)
	for element in selected_species:
	  print(str(element), end=',')
	"""
}

selected_species.map{it -> it.split(",")}.flatMap().set{my_try}
all_input_info_raw.join(my_try).set{all_input_info}

/*
 * Filter input files for genes belonging to the gene cluster of interest
 */

gene_clusters = file(params.geneclusters)
process subset_input_files {
	tag {"${gene_clusterID}"}
	publishDir "${params.output}/${params.geneID}/inputs", mode: 'symlink'
	input:
	val(gene_clusterID) from gene_clusterID
	file(gene_clusters)
	set species, file(annotations), file(overlap_info), file(exint_num_info_all) from all_input_info
	output:
	set species, file("*_subsetted_annot.gtf"), file("*_subsetted_overlap_info.txt") into annotations_overlap_info
	set species, file("*_subsetted_annot.gtf") into (annotations_4_strand, annotations_4_refphases, annotations_4_exphases, annotations_4_final_phases, annotations_4_isoforms, annotations_4_isoforms1)
	set species, file("*_subsetted_ref_exon_int_number.txt") into exint_num_info

	script:
	"""
#!/usr/bin/env python
import pandas as pd
import re

my_cluster_df = pd.read_csv("${gene_clusters}", sep="\t", header=None, index_col=False)
selected_genes = list(my_cluster_df[my_cluster_df.iloc[:,0]=="${gene_clusterID}"].iloc[:,2])
#Subset GTF
my_gtf_df = pd.read_table("${annotations}", header=None, index_col=False)
my_raw_gene_id = [part for element in list(my_gtf_df.iloc[:,my_gtf_df.shape[1]-1]) for part in element.split(";") if "gene_id" in part]
my_gene_id = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_gene_id]
my_gtf_df["GeneID"] = my_gene_id
my_gtf_sub_df = my_gtf_df.loc[my_gtf_df.GeneID.isin(selected_genes)]
my_gtf_sub_df = my_gtf_sub_df.drop(columns=["GeneID"])
my_gtf_sub_df.to_csv("${species}_subsetted_annot.gtf", header=False, index=False, sep="\t", na_rep="NA")
#Subset Overlap file
my_overlap_df = pd.read_table("${overlap_info}", header=None, index_col=False, sep="\t")
my_overlap_sub_df = my_overlap_df[my_overlap_df.iloc[:,1].isin(selected_genes)]
my_overlap_sub_df.to_csv("${species}_subsetted_overlap_info.txt", header=False, index=False, sep="\t", na_rep="NA")
#Subset exint_info
my_exintnum_df = pd.read_table("${exint_num_info_all}", sep="\t", header=0, index_col=False, names=["REF_PROT", "N_ex", "N_int"])
my_geneid = [re.sub(".*\\|", "", element) for element in list(my_exintnum_df["REF_PROT"])]
my_exintnum_df["GeneID"] = my_geneid
my_exintnum_sub_df = my_exintnum_df[my_exintnum_df["GeneID"].isin(selected_genes)]
my_exintnum_sub_df = my_exintnum_sub_df.drop(columns=["GeneID"])
my_exintnum_sub_df.to_csv("${species}_subsetted_ref_exon_int_number.txt", sep="\t", header=True, index=False, na_rep="NA")
#Subset Best hits
#my_besthits_df = pd.read_table("${best_hits_input}", header=0, index_col=False, sep="\t")
#my_besthits_sub_df = my_besthits_df[my_besthits_df.CID=="${gene_clusterID}"]
#my_besthits_sub_df.to_csv("Best_hits_subsetted", sep="\t", header=True, index=False, na_rep="NA")
	"""
}

//Depending on how big this file is, we might need different RAM requirements.
best_hits_input = file(params.bestscores)
process subset_best_hits {
	tag {"${gene_clusterID}"}
	label 'big_mem'
	publishDir "${params.output}/${params.geneID}/inputs", mode: 'symlink'
	input:
	val(gene_clusterID) from gene_clusterID1
	file(best_hits_input)
	output:
	file(Best_hits_subsetted) into best_hits
	script:
	"""
#!/usr/bin/env python

import pandas as pd
my_besthits_df = pd.read_table("${best_hits_input}", header=0, index_col=False, sep="\t")
my_besthits_sub_df = my_besthits_df[my_besthits_df.CID=="${gene_clusterID}"]
my_besthits_sub_df.to_csv("Best_hits_subsetted", sep="\t", header=True, index=False, na_rep="NA")
	"""
}

//Still to adapt to the single gene view.
/*
 * Select only one exon from each group of overlapping exons
 */

process overlapping_exons_info {
    //The tag directive associates each process execution with a custom label.
    //It is used to generate the logs.
    //From what I see, the tag is usually some variable defined in the input?
    tag {"${species}"} 
    publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink' //in theory simlink is the default

	input:
	set species, file(GTF), file(overlap) from annotations_overlap_info

	output:
    	set val(species), file("*_exons_overlap_info.tab") into (processed_overlap_info, processed_overlap_info1, processed_overlap_info2, processed_overlap_info3)
    	//The header of the resulting file is: GeneID, OverlapID, Exon_coords, exon_frequency (in isoforms), exon length 

	script:
	"""
#!/usr/bin/env python

import pandas as pd
import re
import collections

my_overlap_df = pd.read_table("${overlap}", sep="\t", header=None, names=["OverlapID", "GeneID", "ExCoords"])
my_gtf = pd.read_table("${GTF}", sep="\t", header=None)
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
my_gtf_exons = my_gtf[my_gtf.iloc[:,2]=="exon"] #Create a dictionary with key=coords, value=freq
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
my_overlap_df.to_csv("${species}_exons_overlap_info.tab", sep="\t", header=False, index=False, na_rep="NA")
	"""
}

process overlapping_exons_selection {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'

	input:
	set species, file(processed_overlap) from processed_overlap_info
	output:
	set species, file("*_exons_overlap_selected.tab") into (selected_exons_4_strand, selected_exons_4_isoforms, selected_exons_4_strand1, selected_exons_4_strand2, selected_exons_4_strand3)

	script:
	"""
	python ${baseDir}/bin/select_overlapping_exons.py -i ${processed_overlap} -o ${species}_exons_overlap_selected.tab
	"""
}

/*
 * Isolate exons phase phases
 */

//Generate input channel
annotations_4_refphases.join(exint_num_info).set{ref_phases_input}

//Isolate phases from the reference protein
//The header of the output file is: exonID (chr:start-stop), phase.
process ref_protein_phases {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'

	input:
	set species, file(GTF), file(int_ex_num) from ref_phases_input
	output:
	set species, file("*_exons_RefProt_phases.tab") into ref_prot_phases

	script:
	"""
#!/usr/bin/env python

import pandas as pd
import re

my_gtf = pd.read_table("${GTF}", sep="\t", header=None)
my_ex_int_num_df = pd.read_table("${int_ex_num}", sep="\t", header=0)

#select only protein coding
my_gtf = my_gtf[my_gtf.iloc[:,my_gtf.shape[1]-1].str.contains("protein_id")]
my_gtf_subset = my_gtf.iloc[:,my_gtf.shape[1]-1]
my_gtf_subset = my_gtf_subset[my_gtf_subset.str.contains("protein_id")]
my_raw_prot_id = [part for element in list(my_gtf_subset) for part in element.split(";") if "protein_id" in part]
my_gtf["proteinID"] = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_prot_id]

ref_proteins_list = [re.sub("\\|.*", "", element) for element in list(my_ex_int_num_df["REF_PROT"])]
my_int_gtf = my_gtf.loc[my_gtf["proteinID"].isin(ref_proteins_list)] 
my_coords = [str(element)+":"+str(element1)+"-"+str(element2) for element, element1, element2 in zip(list(my_int_gtf.iloc[:,0]), list(my_int_gtf.iloc[:,3]), list(my_int_gtf.iloc[:,4]))]
my_final_gtf = pd.concat([pd.Series(my_coords), pd.Series(list(my_int_gtf.iloc[:,7]))], axis=1) #get a dataframe with exonID, exonPhase
my_final_gtf.to_csv("${species}_exons_RefProt_phases.tab", sep="\t", header=False, index=False, na_rep="NA")
	"""
}

//Generate input channel
annotations_4_exphases.join(ref_prot_phases).set{intron_phases_input}


//Associate phases to each intron/exon
//In the cases where an exon is annotated with different phases depending on the isoform, I select the phase that the exon has in the reference protein.
//This is what the previous file is useful for.
//The header of the output file is: exonID (chr:start-stop), phase.
process isolate_intron_phases {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'

	input:
	set species, file(GTF), file(ref_phases) from intron_phases_input
	output:
	set species, file("*_GTF_intron_phases") into GTF_phases

	script:
	"""
#!/usr/bin/env python

import pandas as pd
import re
import collections

my_gtf = pd.read_table("${GTF}", sep="\t", header=None)
my_phases = pd.read_table("${ref_phases}", sep="\t", header=None, names=["Coords", "Phase"])

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
my_final_phases.to_csv("${species}_GTF_intron_phases", sep="\t", header=False, index=False, na_rep="NA")
	"""
}



/*
 * Add strand and exon phases
 */
//Generate input channel
selected_exons_4_strand.join(annotations_4_strand).join(GTF_phases).set{add_strand_phases_input}
//As before, I am removing all those genes whose isoforms are annotated on different strands (for some reasons)
//The header of the output is: GeneID, ExonID, Strand
process add_exon_strand_phases {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
	input:
	set species, file(selected_exons), file(GTF), file(exon_phases) from add_strand_phases_input
	output:
	set species, file("*_exons_overlap_selected_phases.tab") into added_exon_phases
	script:
	"""
#!/usr/bin/env python

import pandas as pd
import re
import collections

my_selected_exons = pd.read_table("${selected_exons}", sep="\t", header=0)
my_gtf = pd.read_table("${GTF}", sep="\t", header=None)
my_exon_phases_df = pd.read_table("${exon_phases}", sep="\t", header=None, names=["Coords", "Phase"]) 

my_gtf_subset = my_gtf.iloc[:,my_gtf.shape[1]-1] #select the last columns of the GTF
my_raw_gene_id = [part for element in list(my_gtf_subset) for part in element.split(";") if "gene_id" in part] #select geneID
my_gtf["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in my_raw_gene_id]

#Remove genes with exons annotated on different strands
geneID_strand_df = my_gtf.iloc[:,[6,my_gtf.shape[1]-1]].drop_duplicates() #If a gene has exons annotated on both strands, the geneID will be duplicated. 
selected_geneIDs = [item for item, count in collections.Counter(list(geneID_strand_df["geneID"])).items() if count == 1]
my_gtf = my_gtf.loc[my_gtf["geneID"].isin(selected_geneIDs)]
#my_gtf["Coords"] = [str(element1)+"-"+str(element2) for element1, element2 in zip(list(my_gtf.iloc[:,3]), list(my_gtf.iloc[:,4]))]
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
#my_coords_dict = pd.Series(my_gtf_filt.CompleteCoords.values, index=my_gtf_filt.Coords).to_dict()
#my_selected_exons["ExonID"] = my_selected_exons["ExonID"].map(my_coords_dict)
my_final_df = my_selected_exons[["geneID", "ExonID", "Phase", "Strand"]]
my_final_df.to_csv("${species}_exons_overlap_selected_phases.tab", sep="\t", header=False, index=False, na_rep="NA")
	"""
}

//Distinguish between upstream and downstream phases
process add_updown_phases {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
	input:
	set species, file(exons) from added_exon_phases
	output:
	set species, file("*_exons_UpDown_phases.tab") into updown_phases

	script:
	"""
	python ${baseDir}/bin/add_updown_phases.py -i ${exons} -o ${species}_exons_UpDown_phases.tab
	"""
}

/*
 * Add exon cluster information
 */
clusterfile = file(params.exonclusters)
process add_exon_cluster_info {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
	input:
	file(clusterfile)
	set species, file(exons) from updown_phases
	output:
	set species, file("*_exons_cluster_info.tab") into ex_cluster_info

	script:
	"""
	python ${baseDir}/bin/add_exon_clusters_info.py -i ${exons} -c ${clusterfile} -o ${species}_exons_cluster_info.tab -s ${species}
	"""
}

//These are all the files we will use for when each species is a query species.
/*
 * Add introns and fake coords (for plotting) to query species
 */

process add_introns_and_fake_coords {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/processed_tables", mode: 'copy'
	input:
	set species, file(cluster_info) from ex_cluster_info
	output:
	file("*_exons-introns_cluster_info.tab") into ex_int_cluster_info

	script:
	"""
	python ${baseDir}/bin/add_introns_and_fake_coords.py -i ${cluster_info} -o ${species}_exons-introns_cluster_info.tab
	"""
}

/*
 * Add information from best hits
 */
//We need this to eventually select a hit by sequence when no orthology is found.
//With the next 2 processes, I am generating temporary files which should be removed
//But the best hits scores files are too big to be read every single time

my_geneID = "${params.geneID}"
gene_clusters = file(params.geneclusters)
process isolate_query_species {
	tag {"${geneID}"}
	input:
	val(my_geneID)
	file(gene_clusters)
	output:
	stdout into (query_species, query_species1)
	script:
	"""
	#!/usr/bin/env python

	import pandas as pd
	
	my_df = pd.read_csv("${gene_clusters}", sep="\t", header=None, index_col=False)
	query_species = list(my_df[my_df.iloc[:,2]=="${my_geneID}"].iloc[:,1])[0]
	print(str(query_species), end='')
	"""
}

//Generate a channel with all combinations of species query with each species target
//Take the species from GTFs name
def my_query = "${query_species1}".toString()
Channel
    .fromFilePairs( params.annotations, size: 1).map{it[0]}.flatMap()
    .toList().map{[it, it].combinations().findAll{a,b -> a!=b && a=="Hs2"}}
    .flatMap()
    .map{"${it[0]}_${it[1]}".toString()}
    .into{all_species_pairs; all_species_pairs1; all_species_pairs2}


process break_besthits_speciespair {
	tag {"${species_pair}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
	input:
	file(best_hits) from best_hits
	val species_pair from all_species_pairs
	output:
	set species_pair, file("*-best_scores_hits_exons.txt") into best_hits_speciespairs
	
	script:
	def single_species = "${species_pair}".split("_")
	"""
	awk 'NR==1' ${best_hits} > ${species_pair}-best_scores_hits_exons.txt
	cat ${best_hits} | awk '\$16=="${single_species[0]}" && \$17=="${single_species[1]}"' >> ${species_pair}-best_scores_hits_exons.txt	
	"""
}

//Generate input channel. This part is quite ugly.
//Channel with lists [species_pair, species1, species2]
all_species_pairs1.map{"${it}".split("_")}.map{["${it[0]}_${it[1]}".toString(), "${it[0]}", "${it[1]}"]}.into{first_channel; first_channel1}
//join with besthits file on the species_pair
first_channel.join(best_hits_speciespairs).map{[it[1], it[0], it[2], it[3]]}.into{second_channel; second_channel1}
second_channel.combine(processed_overlap_info1).filter{it[0]==it[4]}.into{intermediate; intermediate1}
intermediate.combine(selected_exons_4_strand1).filter{it[0]==it[6]}.map{[it[2], it[0], it[1], it[3], it[5], it[7]]}.into{third_channel; third_channel1}
third_channel.combine(processed_overlap_info2).filter{it[0]==it[6]}.combine(selected_exons_4_strand2).filter{it[0]==it[8]}.map{[it[2], it[3], it[4], it[5], it[7], it[9]]}.into{add_besthits_input; add_besthits_input1}

//Add the besthits scores
process add_best_hits_scores {
	tag {"${species_pair}"}
	publishDir "${params.output}/${params.geneID}/processed_tables", mode: 'copy'
	input:
	set species_pair, file(best_hits_scores), file(over_query), file(over_target), file(selected_over_query), file(selected_over_target) from add_besthits_input
	output:
	file("*-best_scores_hits_exons.translated_coords") into (translated_coords)
	script:
	"""
	python ${baseDir}/bin/translate_overlapping_coords.py -b $best_hits_scores -q $over_query -t $over_target -sq $selected_over_query -st $selected_over_target -o ${species_pair}-best_scores_hits_exons.translated_coords
	"""	
}

/*
 * Isoform information section
 */
//I merged two rules from the snakemake here
process isoform_exon_numbers {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
	input:
	set species, file(GTF) from annotations_4_isoforms
	output:
	set species, file("*_exon_number_by_isoform") into isoform_tot_ex
	set species, file("*_exon_number_in_isoform") into isoform_ex_index
	script:
	"""
	python ${baseDir}/bin/isoform_info.py -i $GTF -o1 "${species}"_exon_number_by_isoform -o2 "${species}"_exon_number_in_isoform
	"""
}


//There was a small bug before. Now it should be fixed.
geneclusters = file(params.geneclusters)
process first_last_internal_ex_info {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
	input:
	file(geneclusters)
	set species, file(tot_ex), file(ex_index), file(selected_ex) from isoform_tot_ex.join(isoform_ex_index).join(selected_exons_4_isoforms)
	output:
	//set species, file("*_first_exons.tab") into first_exons
	//set species, file("*_last_exons.tab") into last_exons
	file("*_all_exons_positions.tab") into all_ex_positions

	script:
	"""
#!/usr/bin/env python

import pandas as pd
import collections

my_ex_index_df = pd.read_table("${ex_index}", sep="\t", header=None, names=["Coords", "Pos", "GeneID", "TranscriptID"])
my_ex_tot_df = pd.read_table("${tot_ex}", sep="\t", header=None, names=["GeneID", "TranscriptID", "TotEx"])
my_first_ex_list = list(my_ex_index_df.loc[my_ex_index_df["Pos"]==1]["Coords"])
my_last_ex_dict = pd.Series(my_ex_tot_df.TotEx.values, index=my_ex_tot_df.TranscriptID).to_dict()
my_ex_index_df["TotEx"] = my_ex_index_df["TranscriptID"].map(my_last_ex_dict)
my_last_ex_list = list(my_ex_index_df.loc[my_ex_index_df["Pos"]==my_ex_index_df["TotEx"]]["Coords"])
my_internal_ex_list = [element for element in list(my_ex_index_df["Coords"]) if element not in my_first_ex_list+my_last_ex_list]
my_first_last_ex_raw = list(set(my_first_ex_list))+list(set(my_last_ex_list))
my_first_last_ex_list = [element for element, count in collections.Counter(my_first_last_ex_raw).items() if count > 1]

#create dictionary for exon status. This works because update will replace the last value.
my_first_ex_dict = {element : "first" for element in my_first_ex_list}
my_last_ex_dict = {element : "last" for element in my_last_ex_list}
my_firstlast_ex_dict = {element : "first;last" for element in my_first_last_ex_list}
my_internal_ex_dict = {element : "Internal" for element in my_internal_ex_list}
my_all_ex_dict = {**my_first_ex_dict, **my_last_ex_dict, **my_firstlast_ex_dict, **my_internal_ex_dict}

#create dictionary for gene clusters
my_clusters = pd.read_table("${geneclusters}", sep="\t", header=None)
my_gene_clusters_dict = pd.Series(list(my_clusters.iloc[:,0]), index=list(my_clusters.iloc[:,2])).to_dict()

#Join results in final table
my_ex_index_df["Status"] = my_ex_index_df["Coords"].map(my_all_ex_dict)
my_ex_index_df["ClusterID"] = my_ex_index_df["GeneID"].map(my_gene_clusters_dict)
my_final_df = my_ex_index_df[["GeneID", "ClusterID", "Coords", "Status"]].drop_duplicates()
my_final_df.to_csv("${species}_all_exons_positions.tab", sep="\t", header=False, index=False, na_rep="NA")
	"""
}


//I don't care about the order here, I will not really be calling these files.
ex_int_cluster_info.mix(all_ex_positions).collect().mix(translated_coords).collect().set{plot_input}
//Get the order of the species for the plot. It can be either provided by the user or computed from the gene cluster file.
if (params.ordered_species) {
	Channel.value("${params.ordered_species}").set{ordered_target}
} 
else {
	Channel
	    .fromFilePairs(params.annotations, size: 1)
	    .map{"${it[0]}".toString()}.flatMap().collect().map{it -> it.join(",")}.set{all_species}
	gene_clusters = file(params.geneclusters)
	process derived_ordered_species {
		input:
		file(gene_clusters)
		val(all_species)
		val(gene_clusterID) from gene_clusterID2
		output:
		stdout into ordered_target
		script:
		"""
		#!/usr/bin/env python
		import pandas as pd
	
		my_df = pd.read_csv("${gene_clusters}", sep="\t", header=None, index_col=False)
		species_list = list(my_df[my_df.iloc[:,0]=="${gene_clusterID}"].iloc[:,1])
		interesting_species = list(set([element for element in species_list if element in str("${all_species}").split(",")]))
		interesting_species.insert(0, interesting_species.pop(interesting_species.index("Hs2")))
		final_species_list = ",".join(str(element) for element in interesting_species)
		print(final_species_list)
		"""
	}
}

//R script to actually make the plot.
my_geneID = "${params.geneID}"
gene_clusters = file(params.geneclusters) 
process plot_exint {
	tag{"${species}"}
	publishDir "${params.output}/${params.geneID}", mode: 'copy'
	input:
	val(my_geneID)
	val(my_query_species) from query_species1 
	val(ordered_target)
	file(gene_clusters)
	file("*") from plot_input
	output:
	"${baseDir}/exint_plots"
	script:
	"""
	Rscript $baseDir/bin/exint_plotter.R ${my_geneID} ${my_query_species} ${params.output}/${params.geneID}/ ${baseDir}/bin ${gene_clusters} ${ordered_target} 
	"""
}

//###############################
//###### TMP backup #############
//###############################

//process overlapping_exons_info {
//    //The tag directive associates each process execution with a custom label.
//    //It is used to generate the logs.
//    //From what I see, the tag is usually some variable defined in the input?
//    tag {"${species}"} 
//    publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink' //in theory simlink is the default
//
//	input:
//	set species, file(GTF), file(overlap) from annotations_overlap_info
//
//	output:
//    	set val(species), file("*_exons_overlap_info.tab") into (processed_overlap_info, processed_overlap_info1, processed_overlap_info2, processed_overlap_info3)
//    	//The header of the resulting file is: GeneID, OverlapID, Exon_coords, exon_frequency (in isoforms), exon length 
//
//	script:
//	"""
//	cat ${overlap} | awk -v OFS="\t" '{print \$1,\$2,\$3,\$3}' | sed 's/-/\t/' \
//| awk -v OFS="\t" '{print \$1,\$2,\$5,\$4-\$3}' | translate -a -r -k <(cat ${GTF} | sed 's/;/\t/; s/"//g; s/gene_id //' \
//| awk -v OFS="\t" '{print \$9,\$1";"\$7}' | sort | uniq | filter_1col 1 <(cat ${GTF} | sed 's/;/\t/; s/"//g; s/gene_id //' \
//| awk -v OFS="\t" '{print \$9,\$1";"\$7}' | sort | uniq | cut -f1 | sort | uniq -u)) 2 | sed 's/;/\t/' \
//| awk -v OFS="\t" '{print \$2,\$1,\$5":"\$3,\$4}' | translate -a -v -e 0 <(cat ${GTF} | awk '\$3=="exon"' \
//| awk -v OFS="\t" '{print \$1":"\$4"-"\$5}' | sort | uniq -c | sed 's/^[ \t]*//; s/ /\t/' \
//| awk -v OFS="\t" '{print \$2,\$1}') 3 > ${species}_exons_overlap_info.tab
//	"""
//	}

//Isolate phases from the reference protein
//process ref_protein_phases {
//	tag {"${species}"}
//	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
//
//	input:
//	set species, file(GTF), file(int_ex_num) from ref_phases_input
//	output:
//	set species, file("*_exons_RefProt_phases.tab") into ref_prot_phases
//
//	script:
//	"""
//cat $GTF | sed 's/ /_/g' | awk -v OFS="\t" '{print \$9,\$1":"\$4"-"\$5,\$8}' | grep protein_id \
//| sed 's/.*protein_id_//; s/"//g; s/;/\t/' | awk -v OFS="\t" '{print \$(NF-1),\$NF,\$1}' \
//| filter_1col 3 <(cat $int_ex_num | tail -n+2 | cut -f1 | sed 's/|.*//') | awk -v OFS="\t" '{print \$2,\$1}' \
//| sort | uniq | sort -k2,2 | uniq -u -f1 | awk -v OFS="\t" '{print \$2,\$1}' > ${species}_exons_RefProt_phases.tab 
//	"""
//}

//Associate phases to each intron/exon
//process isolate_intron_phases {
//	tag {"${species}"}
//	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
//
//	input:
//	set species, file(GTF), file(ref_phases) from intron_phases_input
//	output:
//	set species, file("*_GTF_intron_phases") into GTF_phases
//
//	script:
//	"""
//	cat <(cat $GTF | awk -v OFS="\t" '\$3=="CDS" {print \$8, \$1":"\$4"-"\$5}' | sort | uniq | sort -k2,2 | uniq -d -f1 \
//| cut -f2 | translate -v -e "NA" -a -r <(cat $ref_phases) 1) <(cat $GTF | awk -v OFS="\t" '\$3=="CDS" {print \$8, \$1":"\$4"-"\$5}' \
//| sort | uniq | sort -k2,2 | uniq -u -f1 | awk -v OFS="\t" '{print \$2,\$1}') > ${species}_GTF_intron_phases  
//	"""
//}

//process add_exon_strand {
//	tag {"${species}"}
//	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
//	input:
//	set species, file(selected_exons), file(GTF) from add_strand_input
//	output:
//	set species, file("*_exons_overlap_selected_strand.tab") into selected_exons_strand
//	script:
//	"""
//	cat $selected_exons | tail -n+2 | cut -f1,3 | translate -a -r <(cat $GTF | sed 's/;/\t/; s/"//g; s/gene_id //' \
//| awk -v OFS="\t" '{print \$9,\$7}' | sort | uniq | filter_1col 1 <(cat $GTF | sed 's/;/\t/; s/"//g; s/gene_id //' | awk -v OFS="\t" '{print \$9,\$7}' \
//| sort | uniq | cut -f1 | sort | uniq -u)) 1 > ${species}_exons_overlap_selected_strand.tab
//	"""
//}

//Generate input channel
//selected_exons_strand.join(GTF_phases).set{phases_exons_comb}
//Add phases to the selected exons file
//The snakemake was working, but I am doing something horrible here.
//I will need to recheck all the code step by step
//process add_intron_phases {
//	tag {"${species}"}
//	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
//	input:
//	set species, file(selected_exons), file(phases) from phases_exons_comb
//	output:
//	set species, file("*_exons_overlap_selected_phases.tab") into added_exon_phases
//
//	script:
//	"""
//#!/usr/bin/env python
//
//import pandas as pd
//
//my_selected_exons = pd.read_table("${selected_exons}", sep="\t", header=)
//
//	cat $selected_exons | translate -a -v -e 0 <(cat $phases) 2 > ${species}_exons_overlap_selected_phases.tab
//	"""
//}

//process add_intron_phases {
//	tag {"${species}"}
//	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
//	input:
//	set species, file(selected_exons), file(phases) from phases_exons_comb
//	output:
//	set species, file("*_exons_overlap_selected_phases.tab") into added_exon_phases
//
//	script:
//	"""
//	cat $selected_exons | translate -a -v -e 0 <(cat $phases) 2 > ${species}_exons_overlap_selected_phases.tab
//	"""
//}

//process isolate_first_last_ex {
//	tag {"${species}"}
//	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
//	input:
//	set species, file(tot_ex), file(ex_index) from isoform_tot_ex.join(isoform_ex_index)
//	output:
//	set species, file("*_first_exons.tab") into first_exons
//	set species, file("*_last_exons.tab") into last_exons
//
//	script:
//	"""
//	cat $ex_index | awk '\$2==1' > "${species}"_first_exons.tab
//cat $ex_index | awk -v OFS="\t" '{print \$1,\$3";"\$4";"\$2}' \
//| filter_1col 2 <(cat $tot_ex | awk '{print \$1";"\$2";"\$3}') \
//| tr ";" "\t" | awk -v OFS="\t" '{print \$1,\$4,\$3,\$2}' > "${species}"_last_exons.tab 
//	"""
//}



//Here I am merging three rules which were separated in the snakemake
//geneclusters = file(params.geneclusters)
//process first_last_internal_ex_info {
//	tag{"${species}"}
//	publishDir "${params.output}/${params.geneID}/processed_tables", mode: 'copy'
//	input:
//	file(geneclusters)
//	set species, file(first_ex), file(last_ex), file(selected_ex) from first_exons.join(last_exons).join(selected_exons_4_isoforms)
//	output:
//	file("*_all_exons_positions.tab") into all_ex_positions
//
//	script:
//	"""
//	cat $first_ex | cut -f1,4 | sort | uniq |  awk -v OFS="\t" '{print \$1,\$2,"first"}' > output.tmp
//	cat $last_ex | cut -f1,4 | sort | uniq |  awk -v OFS="\t" '{print \$1,\$2,"last"}' >> output.tmp
//	cat $selected_ex | tail -n+2 | cut -f1,3 | translate -a -d -v -e "Internal" <(cat output.tmp | cut -f1,3) 2 \
//	| translate -a -n -k <(cat $geneclusters | cut -f1,3 | awk '\$1!=""') 1 > "$species"_all_exons_positions.tab
//	rm output.tmp
//	"""
//}
