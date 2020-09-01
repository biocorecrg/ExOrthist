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
ref proteins info (pipeline output)	: ${params.refprot}
exon clusters (pipeline output)		: ${params.exonclusters}
exon best scores (pipeline output)	: ${params.bestscores}
gene clusters (same as pipeline input)	: ${params.geneclusters}
geneID					: ${params.geneID}
relevant exons				: ${params.relevant_exs}
reclustered gene orthology file:	: ${params.sub_orthologs}
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

//Create channel for files with ref proteins info
Channel
    .fromFilePairs(params.refprot, size: 1)
    .ifEmpty{error "Cannot find any overlap info: ${params.refprot}"}
    .set{refprot_info}

//create a joint channel where each key is paired with the corresponding files
annotations.join(overlap_info).join(refprot_info).set{all_input_info_raw}

my_geneID = "${params.geneID}"
if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
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

selected_species.map{it -> it.split(",")}.flatMap().set{single_species}
all_input_info_raw.join(single_species).set{all_input_info}


/*
 * Filter input files for genes belonging to the gene cluster of interest
 */

//gene_clusters = file(params.geneclusters)
if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
process subset_input_files {
	tag {"${gene_clusterID}"}
	publishDir "${params.output}/${params.geneID}/inputs", mode: 'symlink'
	input:
	val(gene_clusterID) from gene_clusterID
	file(gene_clusters)
	set species, file(annotations), file(overlap_info), file(ref_proteins) from all_input_info
	output:
	set species, file("*_subsetted_annot.gtf"), file("*_subsetted_overlap_info.txt") into annotations_overlap_info
	set species, file("*_subsetted_annot.gtf") into (annotations_4_strand, annotations_4_refphases, annotations_4_exphases, annotations_4_isoforms, annotations_4_annots)
	set species, file("*_subsetted_ref_proteins.txt") into ref_prot_info

	script:
	"""
	python ${baseDir}/bin/subset_inputs.py -a ${annotations} -o ${overlap_info} -r ${ref_proteins} -c ${gene_clusters} -g ${gene_clusterID} -s ${species}  
	"""
}

//Depending on how big this file is, we might need different RAM requirements.
//I am using bash because it's much faster than python here.
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
awk 'NR==1' ${best_hits_input} > Best_hits_subsetted
cat ${best_hits_input} | awk '\$1=="${gene_clusterID}"' >> Best_hits_subsetted
"""
}


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
	set species, file(annotations), file(overlap_info) from annotations_overlap_info

	output:
    	set val(species), file("*_exons_overlap_info.tab") into (processed_overlap_info, processed_overlap_info1, processed_overlap_info2)
    	//The header of the resulting file is: GeneID, OverlapID, Exon_coords, exon_frequency (in isoforms), exon length 
	script:
	"""
	python ${baseDir}/bin/overlapping_exons_info.py -a ${annotations} -o ${overlap_info} -out ${species}_exons_overlap_info.tab
	"""
}

process overlapping_exons_selection {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'

	input:
	set species, file(processed_overlap) from processed_overlap_info
	output:
	set species, file("*_exons_overlap_selected.tab") into (selected_exons_4_strand, selected_exons_4_isoforms, selected_exons_4_strand1, selected_exons_4_strand2)

	script:
	"""
	python ${baseDir}/bin/select_overlapping_exons.py -i ${processed_overlap} -out ${species}_exons_overlap_selected.tab
	"""
}

/*
 * Isolate exons phase phases
 */

//Generate input channel
annotations_4_refphases.join(ref_prot_info).set{ref_phases_input}

//Isolate phases from the reference protein
//The header of the output file is: exonID (chr:start-stop), phase.
process isolate_refprotein_exons_phases {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'

	input:
	set species, file(annotations), file(ref_prots) from ref_phases_input
	output:
	set species, file("*_exons_RefProt_phases.tab") into ref_prot_phases

	script:
	"""
	python ${baseDir}/bin/isolate_refprotein_exons_phases.py -a ${annotations} -r ${ref_prots} -out ${species}_exons_RefProt_phases.tab
	"""
}

//Generate input channel
annotations_4_exphases.join(ref_prot_phases).set{intron_phases_input}

//Associate phases to each intron/exon
//In the cases where an exon is annotated with different phases depending on the isoform, I select the phase that the exon has in the reference protein.
//This is what the previous file is useful for.
//The header of the output file is: exonID (chr:start-stop), phase.
process isolate_all_exons_phases {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'

	input:
	set species, file(annotations), file(ref_phases) from intron_phases_input
	output:
	set species, file("*_all_exons_phases") into all_phases

	script:
	"""
	python ${baseDir}/bin/isolate_all_exons_phases.py -a ${annotations} -r ${ref_phases} -out ${species}_all_exons_phases
	"""
}


/*
 * Add strand and exon phases
 */
//Generate input channel
selected_exons_4_strand.join(annotations_4_strand).join(all_phases).set{add_strand_phases_input}
//As before, I am removing all those genes whose isoforms are annotated on different strands (for some reasons)
//The header of the output is: GeneID, ExonID, Strand
process add_strand_and_phase {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
	input:
	set species, file(selected_exons), file(annotations), file(exon_phases) from add_strand_phases_input
	output:
	set species, file("*_exons_overlap_selected_phases.tab") into added_exon_phases
	script:
	"""
	python ${baseDir}/bin/add_strand_and_phase.py -e ${selected_exons} -a ${annotations} -p ${exon_phases} -out ${species}_exons_overlap_selected_phases.tab
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
	python ${baseDir}/bin/add_updown_phases.py -i ${exons} -out ${species}_exons_UpDown_phases.tab
	"""
}

/*
 * Add exon annotation status
 */
process add_exon_annotation_status {
	tag {"${$species}"}
	publishDir "${params.output}/${params.geneID}/intermediate_files", mode: 'symlink'
	input:
	set species, file(annotations), file(exons) from annotations_4_annots.join(updown_phases)
	output:
	set species, file("*_exons_annotation_status") into annotation_status
	script:
	"""
	python ${baseDir}/bin/add_annotation_status.py -a ${annotations} -e ${exons} -out ${species}_exons_annotation_status
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
	//set species, file(exons) from updown_phases
	set species, file(exons) from annotation_status
	output:
	set species, file("*_exons_cluster_info.tab") into ex_cluster_info

	script:
	"""
	python ${baseDir}/bin/add_exon_clusters_info.py -i ${exons} -c ${clusterfile} -out ${species}_exons_cluster_info.tab -s ${species}
	"""
}


//These are all the files we will use for when each species is a query species.
/*
 * Add introns and fake coords (for plotting) to query species
 */

process add_fake_coords {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/processed_tables", mode: 'copy'
	input:
	set species, file(cluster_info) from ex_cluster_info
	output:
	file("*_exons_cluster_info-fakecoords.tab") into ex_int_cluster_info

	script:
	"""
	python ${baseDir}/bin/add_fake_coords.py -i ${cluster_info} -out ${species}_exons_cluster_info-fakecoords.tab
	"""
}

/*
 * Add information from best hits
 */
//We need this to eventually select a hit by sequence when no orthology is found.
//With the next 2 processes, I am generating temporary files which should be removed
//But the best hits scores files are too big to be read every single time

my_geneID = "${params.geneID}"
//gene_clusters = file(params.geneclusters)
if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
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
//def my_query = "${query_species1}".toString()
Channel
    .fromFilePairs( params.annotations, size: 1).map{it[0]}.flatMap()
    //.toList().map{  def my_query = "${query_species1}" [it, it].combinations().findAll{a,b -> a!=b && a==my_query}}
    .toList().map{[it, it].combinations().findAll{a,b -> a!=b && a=="Hs2"}}
    .flatMap()
    .map{"${it[0]}_${it[1]}".toString()}
    .into{all_species_pairs; all_species_pairs1}


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
process isoform_info {
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
//gene_clusters = file(params.geneclusters)
if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
process exon_position_info {
	tag {"${species}"}
	publishDir "${params.output}/${params.geneID}/processed_tables", mode: 'copy'
	input:
	file(gene_clusters)
	set species, file(tot_ex), file(ex_index), file(selected_ex) from isoform_tot_ex.join(isoform_ex_index).join(selected_exons_4_isoforms)
	output:
	file("*_all_exons_positions.tab") into all_ex_positions

	script:
	"""
	python ${baseDir}/bin/exon_position_info.py -t ${tot_ex} -i ${ex_index} -g ${gene_clusters} -out ${species}_all_exons_positions.tab
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
	//gene_clusters = file(params.geneclusters)
	if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
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
		print(final_species_list, end='')
		"""
	}
}

//R script to actually make the plot.
my_geneID = "${params.geneID}"
//gene_clusters = file(params.geneclusters)
if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}
if (params.relevant_exs) {relevant_exons = file("${params.relevant_exs}")} else {
	relevant_exons = ""
}
 
process plot_exint {
	tag{"${species}"}
	publishDir "${params.output}/${params.geneID}", mode: 'copy'
	input:
	val(my_geneID)
	val(my_query_species) from query_species1 
	val(ordered_target)
	file(gene_clusters)
	file(relevant_exons)
	file("*") from plot_input
	output:
	"${baseDir}/exint_plots"
	script:
	"""
	Rscript $baseDir/bin/exint_plotter.R ${my_geneID} ${my_query_species} ${params.output}/${params.geneID}/ ${baseDir}/bin ${gene_clusters} ${ordered_target} "${relevant_exons}" 
	"""
}

//process isolate_cluster_id {
//	tag {"${geneID}"}
//	input:
//	val(my_geneID)
//	file(gene_clusters)
//	output:
//	stdout into (gene_clusterID, gene_clusterID1, gene_clusterID2, gene_clusterID3)
//	script:
//	"""
//	#!/usr/bin/env python
//
//	import pandas as pd
//	
//	my_df = pd.read_csv("${gene_clusters}", sep="\t", header=None, index_col=False)
//	clusterID = list(my_df[my_df.iloc[:,2]=="${my_geneID}"].iloc[:,0])[0]
//	print(clusterID, end='')
//	"""
//}
