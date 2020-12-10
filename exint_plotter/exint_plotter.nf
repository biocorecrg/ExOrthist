#!/usr/bin/env nextflow

/*
 * Exint plotter for Exhortist pipeline
 *
 */


/*
 * Input parameters:
 */

log.info  """
Executing with the following parameters:

output main				: ${params.output_main}
geneID					: ${params.geneID}
isoformID				: ${params.isoformID}
relevant exons				: ${params.relevant_exs}
reclustered gene orthology file:	: ${params.sub_orthologs}
"""

/*
 * Define params starting from the main.nf output folder
 */


params.geneclusters = "${params.output_main}/gene_cluster_file.gz"
params.annotations = "${params.output_main}/*/*_annot_fake.gtf.gz"
params.overlap = "${params.output_main}/*/*_overlap_CDS_exons.txt"
params.refprot = "${params.output_main}/*/*_ref_proteins.txt"
params.exonclusters = "${params.output_main}/EX_clusters.tab"
params.bestscores = "${params.output_main}/*/best_scored_EX_matches_by_targetgene.txt" //This are all unfiltered scores. I need to identify exons matched by sequence conservation but not phased conservation.

/*
 * Create channels for input data
 */

//This channel will contain a list of the GTF files, in theory each with a key
//The key corresponds to the value assumed by the wildcard in the annotation variable (which is defined in the params.config)
//annotations  = "$baseDir/data/GTF/*_annot.gtf"
Channel 
    .fromFilePairs(params.annotations, size: 1)
    .ifEmpty{error "Cannot find any annotation matching: ${params.annotations}"}
    .set{annotations}

//The key is the species, same as for the annotations channel
Channel
    .fromFilePairs( params.overlap, size: 1)
    .ifEmpty{error "Cannot find any overlap info: ${params.overlap}"}
    .set{overlap_info}

//Create channel for files with ref proteins info
Channel
    .fromFilePairs(params.refprot, size: 1)
    .ifEmpty{error "Cannot find any overlap info: ${params.refprot}"}
    .set{refprot_info}

//Create a joint channel where each key is paired with the corresponding files
//annotations.join(overlap_info).join(refprot_info).into{all_input_info_raw; all_input_info_raw1}
annotations.join(overlap_info).join(refprot_info)map{it.flatten()}.set{all_input_info_raw}

Channel
    .fromPath(params.bestscores)
    .toList()
    .ifEmpty{error "Cannot find any overlap info: ${params.bestscores}"}
    .set{best_hits_input}

/*
 * Isolate clusterID, query_species and species with orthologs
 */
my_geneID = "${params.geneID}"
if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
process isolate_cluster_id {
	tag {"${my_geneID}"}
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

my_geneID = "${params.geneID}"
if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
process isolate_query_species {
	tag {"${my_geneID}"}
	input:
	val(my_geneID)
	file(gene_clusters)
	output:
	stdout into (query_species, query_species1, query_species2, query_species3)
	script:
	"""
	#!/usr/bin/env python

	import pandas as pd
	
	my_df = pd.read_csv("${gene_clusters}", sep="\t", header=None, index_col=False)
	query_species = list(my_df[my_df.iloc[:,2]=="${my_geneID}"].iloc[:,1])[0]
	print(str(query_species), end='')
	"""
}

/*
 * Filter only for the species actually having exons in the exon clusters.
 */
//The pipeline crushes otherwise. But there is a distinction: if a gene gets plotted with no exons, it means that its exons make it to the exon clusters. They simply do not have any ortholog exons to the query species.
//Interrupt the pipeline in case the query gene is not in the exon clusters
my_exon_clusters = file(params.exonclusters)
process select_species_with_orthologs {
	tag {"${clusterID}"}
	input:
	file(my_exon_clusters)
	val(clusterID) from gene_clusterID3
	val(my_geneID)
	output:
	stdout into selected_species
	script:
	"""
	#!/usr/bin/env python
	import pandas as pd
	import re
	import sys

	my_df = pd.read_table("${my_exon_clusters}", sep="\t", header=0, index_col=False)
	my_df["GeneClusterID"] = [re.sub("\\..*", "", element) for element in my_df["ExCluster_ID"]]
	gene_ids_list = list(my_df[my_df.GeneClusterID=="${clusterID}"]["GeneID"])
	selected_species = list(set(list(my_df[my_df.GeneClusterID=="${clusterID}"]["Species"])))
	if "${my_geneID}" in gene_ids_list and len(selected_species) >=2:
	  for element in selected_species:
	    print(str(element), end=',')
	else:
	  sys.exit("Query gene does not have any exon orthologs in the other species")
	"""
}

selected_species.map{it -> it.split(",")}.flatMap().set{single_species}
all_input_info_raw.join(single_species).set{all_input_info}

/*
 * Subset input files for genes belonging to the gene cluster of interest
 */

if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
process subset_input_files {
	tag {"${gene_clusterID}"}
	input:
	val(gene_clusterID) from gene_clusterID
	file(gene_clusters)
	set species, file(annotations), file(overlap_info), file(ref_proteins) from all_input_info
	output:
	set species, file("*_subsetted_annot.gtf"), file("*_subsetted_overlap_info.txt"), file("*_subsetted_ref_proteins.txt") into all_subsetted_inputs
	set species, file("*_subsetted_overlap_info.txt") into overlap_info_4_isoforms
	script:
	"""
	subset_inputs.py -a ${annotations} -o ${overlap_info} -r ${ref_proteins} -c ${gene_clusters} -g ${gene_clusterID} -s ${species}  
	"""
}

//Depending on how big this file is, we might need different RAM requirements.
//I am using bash because it's much faster than python here.
process subset_best_hits {
	tag {"${gene_clusterID}"}
	label 'big_mem'
	input:
	val(gene_clusterID) from gene_clusterID1
	file("best_hits_species_pairs_*") from best_hits_input
	output:
	file(Best_hits_subsetted) into best_hits
	script:
"""
cat best_hits_species_pairs_* > best_hits.tmp
awk 'NR==1' best_hits.tmp > Best_hits_subsetted
cat best_hits.tmp | awk '\$1=="${gene_clusterID}"' >> Best_hits_subsetted
rm best_hits.tmp
"""
}


/*
 * Generate exint plotter input for each species 
 */
if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
exon_clusters = file(params.exonclusters)
process generate_fake_coords_table {
    tag {"${species}"}
	input:
	set species, file(annotations), file(overlap_info), file(ref_prot) from all_subsetted_inputs
	file(gene_clusters)
	file(exon_clusters)
	output:
	set species, file("*_exons_cluster_info-fakecoords.tab") into fake_coords_tables
	set species, file("*_ExNum_in_isoform") into ExNum_number_in_isoform
	set species, file("*_overlapID_chosenID.txt") into overlapID_chosenID_4plot, overlapID_chosenID_4besthits, overlapID_chosenID_4besthits1
	script:
	"""
	generate_exint_plotter_input.py -a ${annotations} -o ${overlap_info} -c ${gene_clusters} -e ${exon_clusters} -r ${ref_prot} -s ${species} -out ${species}_exons_cluster_info-fakecoords.tab
	"""
}

/*
 * Generate best hits input files
 */

//Generate a channel with all combinations of species query with each species target
Channel
    .fromFilePairs( params.annotations, size: 1).map{it[0]}
    .toList()
    .map{[it, it].combinations().findAll{a,b -> a!=b}}
    .flatMap()
    .groupTuple()
    .join(query_species)
    .transpose().set{species_pairs_tmp}

//Add the needed files from the previous process
overlapID_chosenID_4besthits.cross(species_pairs_tmp).map{it.flatten()}
    .map{[it[3], it[0], it[1], it[2]]}
    .join(overlapID_chosenID_4besthits1)
    .map{["${it[3]}_${it[0]}".toString(), it[2], it[4]]}
    .into{all_species_pairs}

process break_besthits_speciespair {
	tag {"${species_pair}"}
	input:
	file(best_hits) from best_hits
	set species_pair, file(query_chosenID), file(target_chosenID) from all_species_pairs
	output:
	set species_pair, file("*-best_scores_with_overlapIDs.txt") into best_hits_speciespairs
	
	script:
	def single_species = "${species_pair}".split("_")
"""
awk 'NR==1' ${best_hits} > ${species_pair}-best_scores_hits_exons.txt
cat ${best_hits} | awk '\$16=="${single_species[0]}" && \$17=="${single_species[1]}"' >> ${species_pair}-best_scores_hits_exons.txt
add_selected_overlapID_to_besthits.py -b ${species_pair}-best_scores_hits_exons.txt -q ${query_chosenID} -t ${target_chosenID} -out ${species_pair}-best_scores_with_overlapIDs.txt
"""
}

/*
 * Facultative processes for isoform and exons highlighting
 */

fake_coords_tables.collect().mix(best_hits_speciespairs).collect().set{plot_input}
//Get the order of the species for the plot. It can be either provided by the user or computed from the gene cluster file.
if (params.ordered_species) {
	Channel.value("${params.ordered_species}").set{ordered_target}
} 
else {
	Channel
	    .fromFilePairs(params.annotations, size: 1)
	    .map{"${it[0]}".toString()}.flatMap().collect().map{it -> it.join(",")}.set{all_species}
	if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
	process derive_ordered_species {
		input:
		file(gene_clusters)
		val(all_species)
		val(query_species) from query_species3
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
		interesting_species.insert(0, interesting_species.pop(interesting_species.index("${query_species}")))
		final_species_list = ",".join(str(element) for element in interesting_species)
		print(final_species_list, end='')
		"""
	}
}

my_geneID = "${params.geneID}"
if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
if (params.isoformID) {
	isoformID = params.isoformID
	process get_isoform_exons {
		tag{"${isoformID}"}
		publishDir "${params.output}/${params.geneID}/processed_table", mode: 'copy'
		input:
		val(isoformID)
		set species, file(exons_info), file(overlap_info) from ExNum_number_in_isoform.join(overlap_info_4_isoforms).join(query_species1)
		output:
		stdout into (isoform_interesting_exs) 
		script:
		"""
		#!/usr/bin/env python
		import pandas as pd
		import re

		my_df = pd.read_csv("${exons_info}", sep="\t", header=None, index_col=False, names=["ExonID", "ExonNum", "GeneID", "ProteinID"])
		interesting_exons = list(my_df.loc[my_df["ProteinID"]=="${isoformID}"]["ExonID"])
		#select all the exons from the same overlapping group
		my_chr = [re.sub(":.*", "", element) for element in interesting_exons][0]
		my_overlap_df = pd.read_csv("${overlap_info}", sep="\t", header=None, index_col=False, names=["OverlapID", "GeneID", "Start_Stop"])
		my_overlap_df["Coords"] = [my_chr+":"+element for element in list(my_overlap_df["Start_Stop"])]
		overlapping_groups = list(my_overlap_df.loc[my_overlap_df["Coords"].isin(interesting_exons)]["OverlapID"])
		all_interesting_exons = list(my_overlap_df.loc[my_overlap_df["OverlapID"].isin(overlapping_groups)]["Coords"])
		all_interesting_exons = ["${isoformID}"] + all_interesting_exons
		print(",".join(str(element) for element in all_interesting_exons), end='')
		"""
 }
} else {
	Channel.from("None").into{isoform_interesting_exs}
}

//Rscript to actually make the plot
if (params.relevant_exs) {relevant_exons = "${params.relevant_exs}"} else {relevant_exons = "None"}
process plot_exint {
	tag{"${my_geneID}"}
	containerOptions '-B $PWD:/tmp'
	publishDir "${params.output}", mode: 'copy'
	//publishDir "${params.output}/${params.geneID}", mode: 'copy'
	input:
	val(my_geneID)
	val(my_query_species) from query_species2 
	val(ordered_target)
	val(relevant_exons)
	file(gene_clusters)
	val(isoform_interesting_exs)
	file("*") from plot_input
	output:
	file("*_exint_plot.pdf")
	script:
	"""
	Rscript $baseDir/bin/exint_plotter.R ${my_geneID} ${my_query_species} ${baseDir}/bin ${gene_clusters} ${ordered_target} ${isoform_interesting_exs} ${relevant_exons} 
	"""
}
