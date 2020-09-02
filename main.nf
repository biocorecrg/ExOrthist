#!/usr/bin/env nextflow


/*
 * Copyright (c) 2019-2020, Centre for Genomic Regulation (CRG)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


/*
===========================================================
ExOrthist pipeline for Bioinformatics Core @ CRG

 @authors
 Luca Cozzuto <lucacozzuto@gmail.com>
===========================================================
*/

version = '0.1'

/*
 * Input parameters:
*/

params.help            = false
params.resume          = false


log.info """

╔╦╗┬ ┬┌─┐  ╔═╗─┐ ┬╔═╗┬─┐┌┬┐┬ ┬┬┌─┐┌┬┐
 ║ ├─┤├┤   ║╣ ┌┴┬┘║ ║├┬┘ │ ├─┤│└─┐ │ 
 ╩ ┴ ┴└─┘  ╚═╝┴ └─╚═╝┴└─ ┴ ┴ ┴┴└─┘ ┴ 

==============================================================================
annotations (GTF files)          : ${params.annotations}
genomes (fasta files)            : ${params.genomes}
cluster file (txt files)         : ${params.cluster}
intcons (1 or 2)                 : ${params.intcons}
exsim (from 0 to 1)            : ${params.exsim}
maxsize                          : ${params.maxsize}
clusternum (number of clusters)  : ${params.clusternum}
extraexons (e.g. from VastDB)    : ${params.extraexons}
liftover                         : ${params.liftover}
orthofolder                      : ${params.orthofolder}
vastdb                           : ${params.vastdb}
output (output folder)           : ${params.output}
email for notification           : ${params.email}

INFORMATION ABOUT OPTIONS:
- intcons (1 or 2): Whether to consider one or two introns 
     bodering the exon when filtering by conservation.
- exsim (from 0 to 1): Minimum % of similarity between the
     pair of exons and their corresponding upstream and 
     downstream exons.
- maxsize: Maximum size difference between the two exons 
     (as a fraction of either exon).
     
     
"""

if (params.help) {
    log.info """ExOrthist v0.0.1.beta"""
    log.info """ExOrthist is a Nextflow-based pipeline to obtain groups of exon orthologous at all evolutionary timescales.\n"""
    exit 1
}
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

clusterfile       = file(params.cluster)
outputQC          = "${params.output}/QC"
blosumfile        = file("${baseDir}/files/blosum62.txt")

if ( !blosumfile.exists() ) exit 1, "Missing blosum file: ${blosumfile}!"
if ( !clusterfile.exists() ) exit 1, "Missing clusterfile file: ${clusterfile}!"

/*
 * Create channels for sequences data
 */
Channel
    .fromFilePairs( params.genomes, size: 1)
    .ifEmpty { error "Cannot find any genome matching: ${params.genomes}" }
    .set {genomes}

Channel
    .fromFilePairs( params.annotations, size: 1)
    .ifEmpty { error "Cannot find any annotation matching: ${params.annotations}" }
    .set {annotations}


/*
 * Create channels for optional extra exons
 */
if (params.extraexons) {
	Channel
   	 .fromFilePairs( params.extraexons, size: 1)
   	 .ifEmpty { print "Not using extra exons" }
   	 .set {extraexons}
	genomes.join(annotations).into{data_to_annotation_raw; pipe_data}
	data_to_annotation_raw.join(extraexons).set{data_to_annotation}
} else {
    genomes.join(annotations).into{pipe_data; data_to_annotation}
}


/*
 * Generate annotations
 */
if (params.extraexons) {
	process generate_annotations_with_extra_exons {
		tag { genomeid }
		label 'big_cpus'
		publishDir "${params.output}/", mode: 'copy'

		input:
		set genomeid, file(genome), file(annotation), file(extraexons) from data_to_annotation

		output:
		set val(genomeid), file (genomeid) into idfolders

		script:
		"""
		A1_generate_annotations.pl -GTF ${annotation} -G ${genome} -sp ${genomeid} -add_exons ${extraexons}
		"""
	}
} else {
	process generate_annotations {
		tag { genomeid }
		publishDir "${params.output}/", mode: 'copy'

		input:
		set genomeid, file(genome), file(annotation) from data_to_annotation

		output:
		set val(genomeid), file (genomeid) into idfolders

		script:
		"""
		A1_generate_annotations.pl -GTF ${annotation} -G ${genome} -sp ${genomeid}
		"""
	}
}

/*
 * split cluster file
 */
//Copy the gene cluster file to output to use for the exint plotter
process split_clusters_by_species_pairs {
    tag { clusterfile }
    publishDir "${params.output}/", mode: 'copy', pattern: "gene_cluster_file.gz" 

    input:
    file(clusterfile)

    output:
    file "*.cls.tab"  into cls_tab_files, cls_tab_file_4_clustering
    file("gene_cluster_file.gz")

	script:
	"""
   if [ `echo ${clusterfile} | grep ".gz"` ]; then
       zcat ${clusterfile} > gene_cluster_file
       A2_split_clusters_by_species_pairs.pl -f gene_cluster_file
       gzip gene_cluster_file
       #rm cluster_file
    else
       cat ${clusterfile} > gene_cluster_file
       #A2_split_clusters_by_species_pairs.pl -f ${clusterfile}
       A2_split_clusters_by_species_pairs.pl -f gene_cluster_file
       gzip gene_cluster_file
    fi
	"""
}

idfolders
  .toList().map{ [it, it] .combinations().findAll{ a, b -> a[0] < b[0]} }
  .flatMap()
  .map { ["${it[0][0]}-${it[1][0]}".toString(), it[0][1], it[1][1]] }
  .into{cluster_2_split; anno_2_score_ex_int; species_to_recluster_genes}



/*
 * split clusters
*/

process split_clusters_in_chunks {
    tag { id_comb }

    input:
    file(cls_tab_files).collect()
    set id_comb, file(idfolder_A), file(idfolder_B) from cluster_2_split

    output:
    set file(idfolder_A), file(idfolder_B), file("${idfolder_A}_${idfolder_B}/*.cls.tab-part_*") into cls_files_2_align

	script:
	"""
		A3_split_clusters_in_chunks.pl --sp1 ${idfolder_A} --sp2 ${idfolder_B} --expath ./ --project_dir ./ --N_split ${params.clusternum} --gene_cluster ${id_comb}.cls.tab
	"""
}

cls_files_2_align.transpose().set{cls_files_2_align_t}


/*
 * Align pairs
 */

process parse_IPA_prot_aln {
    tag { "${cls_part_file}" }
    label 'big_cpus'

    input:
    file(blosumfile)
    set file(sp1), file(sp2), file(cls_part_file) from cls_files_2_align_t

    output:
   // set val("${sp1}-${sp2}"), file(sp1), file(sp2), file("${sp1}-${sp2}_*/exons_to_realign_part_*.txt") into aligned_subclusters_4_realign_A
   // set val("${sp1}-${sp2}"), file("${sp1}-${sp2}_*/exons_to_realign_part_*.txt") into file_4_realign_A
   // set val("${sp1}-${sp2}"), file("${sp1}-${sp2}_*") into aligned_subclusters_4_merge
	set val("${sp1}-${sp2}"), file(sp1), file(sp2), file("${sp1}-${sp2}_*"), file("${sp1}-${sp2}_*/EXs_to_split_part_*.txt") into aligned_subclusters_4_splitting

	script:
    def cls_parts = "${cls_part_file}".split("_")
	"""
	echo ciao
		B1_parse_IPA_prot_aln.pl ${sp1} ${sp2} ${cls_part_file} \
${sp1}/${sp1}_annot_exons_prot_ids.txt ${sp2}/${sp2}_annot_exons_prot_ids.txt \
${sp1}/${sp1}_protein_ids_exons_pos.txt ${sp2}/${sp2}_protein_ids_exons_pos.txt \
${sp1}/${sp1}_protein_ids_intron_pos_CDS.txt ${sp2}/${sp2}_protein_ids_intron_pos_CDS.txt \
${sp1}/${sp1}.exint ${sp2}/${sp2}.exint ${cls_parts[1]} ${blosumfile} ${sp1}-${sp2}_${cls_parts[1]} ${task.cpus}
	"""
}



/*
 * Split exons pairs to realign
 */

//10 times as many exon alignments as gene clusters in part
//Channel.from("${params.clusternum}").toInteger().map{it*10}.set{ex_aln_per_part}
process split_EX_pairs_to_realign {
    tag { "${folders}" }

    input:
    set comp_id, file(sp1), file(sp2), path(folders), val(req_file) from aligned_subclusters_4_splitting
    //val(ex_aln_per_part)
    output:
    set comp_id, file(sp1), file(sp2), path(folders), file("${folders}/EXs_to_realign_part_*.txt") into aligned_subclusters_4_realign
    script:
    """
	#B2_split_EX_pairs_to_realign.pl ${folders} ${params.clusternum}
	B2_split_EX_pairs_to_realign.pl ${folders} 1000
    """
}




/*
 * Realign exons pairs (with multiple hits)
 */

process realign_EX_pairs {
    tag { "${aligned_folder}/${exons_to_realign.simpleName}" }
    label 'incr_time_cpus'

    input:
    file(blosumfile)
    set val(comp_id), file(sp1), file(sp2), path(aligned_folder), val(exons_to_realign) from aligned_subclusters_4_realign.transpose()

    output:
    set val(comp_id), path(aligned_folder), file("${aligned_folder}/realigned_*") into realigned_exons_4_merge

	script:
    def cls_parts = "${aligned_folder.simpleName}".split("_")
    def infile = exons_to_realign.getFileName()
    def realigned_file = "realigned_" + exons_to_realign.getBaseName() + "_${cls_parts[1]}.txt"
	"""
		B3_realign_EX_pairs.pl ${sp1} ${sp2} ${aligned_folder}/${infile} \
		${sp1}/${sp1}.exint ${sp2}/${sp2}.exint ${cls_parts[1]} ${aligned_folder}/${realigned_file} \
		${sp1}_${sp2} ${blosumfile} ${task.cpus}
	"""
}

realigned_exons_4_merge.map{
	[it[1].getFileName(), it[1]]
}.groupTuple().map{
    def pieces = "${it[0]}".split("_")
	[pieces[0], it[1].first()]
}.groupTuple().set{
	data_4_merge
}

/*
 * Join scores
 */

process join_scores {
    tag { "${comp_id}" }
	stageInMode = 'copy'

    input:
    set comp_id, file("FOLDERS_*") from data_4_merge

    output:
    set val(comp_id), file("${comp_id}/") into folder_jscores

	script:
	"""
	    mkdir ${comp_id}
	    rm FOLDERS_*/exons_to_realign*
	    rm FOLDERS_*/tmp.txt
	    mv FOLDERS_*/* ${comp_id}

    	join_score_files.pl ${comp_id}
        get_best_score_ex_pair.pl ${comp_id}/Score_all_exons.txt ${comp_id}/Best_scores_pair_exons.txt
	"""
}

folder_jscores.join(anno_2_score_ex_int).map{
   [it[0], it[1..-1] ]
}.set{data_to_score}


/*
 * Score exon and introns together
 */
 
process get_all_scores_exon_introns {
    tag { "${comp_id}" }
    label('big_mem')

    input:
    set val(comp_id), file("*") from data_to_score

    output:
    set val(comp_id), file("${comp_id}/Best_score_hits_exons.txt") into bestscore_per_filt

	script:
    def species = comp_id.split("-")
	"""
    get_scores_exons_introns.pl ${species[0]} ${species[1]} \
    ${comp_id}/Aligned_proteins.txt ${comp_id}/Best_scores_pair_exons.txt ${comp_id}/Score_all_introns.txt \
    ${species[0]}/${species[0]}.exint ${species[1]}/${species[1]}.exint \
    ${species[0]}/${species[0]}_protein_ids_intron_pos_CDS.txt ${species[1]}/${species[1]}_protein_ids_intron_pos_CDS.txt \
    ${comp_id}/Final_aln_scores_${comp_id}.txt;
    get_score_by_exon.pl ${comp_id}/Final_aln_scores_${comp_id}.txt ${comp_id}
    #filter_exons_by_score.pl -b ${comp_id}/Best_score_hits_exons.txt -sps ${species[0]},${species[1]} -int ${params.intcons} -id ${params.exsim} -max_size ${params.maxsize}
    """
}

/*
 * Score exon and introns together
 */
process filter_scores {
    tag { "${comp_id}" }

    input:
    set val(comp_id), file(best_score) from bestscore_per_filt

    output:
    file("*.tab") into filterscore_per_joining

	script:
    def species = comp_id.split("-")
	"""
    filter_exons_by_score.pl -b ${best_score} -sps ${species[0]},${species[1]} -int ${params.intcons} -id ${params.exsim} -max_size ${params.maxsize}
    """
}

/*
 * join best scores
 */
process join_best_filtered_scores {
    publishDir "${params.output}/", mode: 'copy'

    input:
    file ("best_score_*") from filterscore_per_joining.collect()

    output:
    file("Best_score_exon_hits_filtered_${params.maxsize}-${params.intcons}-${params.exsim}.tab") into filtered_all_scores

	script:
	"""
    #cat best_score_* >> Best_score_exon_hits_filtered_${params.maxsize}-${params.intcons}-${params.exsim}.tab
    echo "GeneID_sp1\tExon_coords_sp1\tGeneID_sp2\tExon_coords_sp2\tSp1\tSp2" > Best_score_exon_hits_filtered_${params.maxsize}-${params.intcons}-${params.exsim}.tab;
    for i in best_score_*; do cat \$i | tail -n+2 >> Best_score_exon_hits_filtered_${params.maxsize}-${params.intcons}-${params.exsim}.tab; done
    """
}

/*
 * Removing redundant hits
 */
process filter_redundant {

    input:
    file(scores) from filtered_all_scores

    output:
    file("Best_score_exon_hits_pairs.txt") into score_exon_hits_pairs

	script:
	liftfile = ""
	if (params.liftover != "") {
		liftfile = file(params.liftover)
		if ( !liftfile.exists() ) exit 1, "Missing liftover file: ${liftfile}!"

	}
	"""
    get_count_exons.pl ${scores} Exon_count_hits_by_sp.tab ${liftfile}
    get_overlap_exons.pl -i Exon_count_hits_by_sp.tab -o Overlap_exons_by_sp.tab
    Filter_final_exons_pair.pl Overlap_exons_by_sp.tab ${scores} Best_score_exon_hits_pairs.txt ${liftfile}
    """
}

/*
* Split the file of exon pairs
*/
clusterfile = file(params.cluster)
process get_pre_cluster_exons {
    input:
    file(score_exon_hits_pairs)
    //file(cls_tab_file_4_clustering)
    file(clusterfile)

    output:
    file("PART_*/*.tab") into cluster_parts

	script:
	"""
   if [ `echo ${clusterfile} | grep ".gz"` ]; then
       zcat ${clusterfile} > cluster_file
       Get_pre_cluster_exons.pl cluster_file ${score_exon_hits_pairs} 500 out.txt
       rm cluster_file
    else
       Get_pre_cluster_exons.pl ${clusterfile} ${score_exon_hits_pairs} 500 out.txt
    fi
    #Get_pre_cluster_exons.pl ${clusterfile} ${score_exon_hits_pairs} 500 out.txt
	"""
}

/*
* Split the file of exon pairs
*/
process cluster_with_R {

    input:
    file(cluster_part) from cluster_parts.flatten()

    output:
    file("EX${cluster_part}") into ex_clusters

	script:
	"""
    cluster.R ${cluster_part} EX${cluster_part}
    """
}

/*
* Re-clustering of genes
*/
process recluster_genes_species {
    publishDir "${params.output}/", mode: 'copy'
    tag { "${combid}" }

	when:
	params.orthofolder != ""

    input:
    set combid, file(folderA), file(folderB) from species_to_recluster_genes
    file(clusterfile)

    output:
    set combid, file("Recluster_genes_*.tab") into recl_genes_for_rec_exons
    file("Table_recl_exons_*.tab") optional true

	script:
	def combid1 = combid.replace("-", "_")
	def combid2 = combid.replace("-", ",")
	def species_out_file = "Recluster_genes_${combid1}.tab"
	def table_out_file = "Table_recl_exons_${combid1}.tab"
	def vastbcmd = ""
	if (params.vastdb!= "") {
		vastbcmd = "Add_vastids_excls.pl ${params.vastdb} ${species_out_file} ${table_out_file}"
	}
	"""
 		Recluster_genes_sps_pair.pl --cl_file ${clusterfile} --sps ${combid2} --outfile ${species_out_file} --of_path ${params.orthofolder}
    	${vastbcmd}
    """
}


/*
* Final output
*/
process make_final_output {
    publishDir "${params.output}/", mode: 'copy'

    input:
    file("*") from ex_clusters.collect()

    output:
    file("Exon_Clusters.tab") into exon_cluster_for_reclustering
    file("Exon_Clusters_Info.tab.gz")
    file("Table_exon_clusters.tab") optional true

	script:
	def vastbcmd = ""
	if (params.vastdb!= "") {
		vastbcmd = "Add_vastids_excls.pl ${params.vastdb} Exon_Clusters.tab Table_exon_clusters.tab"
	}
	"""
    Get_pre_table_clusters.pl
    ${vastbcmd}
    """
}

/*
* Re-clustering of exons
*/

process recluster_exons_species {
    publishDir "${params.output}/", mode: 'copy'
    tag { "${combid}" }

	when:
	params.orthofolder != ""

    input:
    set combid, file(recl_genes) from recl_genes_for_rec_exons
    file(exon_cluster_for_reclustering)

    output:
    file("Recluster_exons_*.tab")

	script:
	def combid1 = combid.replace("-", "_")
	def species_out_file = "Recluster_exons_${combid1}.tab"

	"""
 		Recluster_exons_sps_pair.pl ${recl_genes} ${exon_cluster_for_reclustering} ${species_out_file}
    """
}

/*
* functions
*/

def getFolderName(sample) {
   folder_info = sample.toString().tokenize("/")
   return folder_info[-2]
}

// make named pipe
def unzipBash(filename) {
    def cmd = filename.toString()
    if (cmd[-3..-1] == ".gz") {
    	cmd = "<(zcat ${filename})"
    }
    return cmd
}

/*
 * Mail notification
 */

if (params.email == "yourmail@yourdomain" || params.email == "") {
    log.info 'Skipping the email\n'
}
else {
    log.info "Sending the email to ${params.email}\n"

    workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
        .stripIndent()

        sendMail(to: params.email, subject: "VectorQC execution", body: msg,  attach: "${outputMultiQC}/multiqc_report.html")
    }
}

workflow.onComplete {
    println "--- Pipeline BIOCORE@CRG ExOrthist ---"
    println "Started at  $workflow.start"
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

