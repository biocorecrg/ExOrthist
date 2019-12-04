#!/usr/bin/env nextflow


/* 
 * Copyright (c) 2019, Centre for Genomic Regulation (CRG) 
 * 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */


/*
===========================================================
vectorQC pipeline for Bioinformatics Core @ CRG

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
Biocore@CRG Yamile's pipeline - N F  ~  version ${version}

╦ ╦╔═╗╔╦╗╦╦  ╔═╗╔═╗  ┌─┐┬┌─┐┌─┐┬  ┬┌┐┌┌─┐
╚╦╝╠═╣║║║║║  ║╣ ╚═╗  ├─┘│├─┘├┤ │  ││││├┤ 
 ╩ ╩ ╩╩ ╩╩╩═╝╚═╝╚═╝  ┴  ┴┴  └─┘┴─┘┴┘└┘└─┘
                                                                                       
==============================================================================
annotations (GTF files)          : ${params.annotations}
genomes (fasta files)     	     : ${params.genomes}
cluster file (txt files)         : ${params.cluster}
intcons (1 or 2)                 : ${params.intcons}
Whether to consider one or two introns bodering the exon 
when filtering by conservation.
idexons (from 0 to 1)            : ${params.idexons}
Minimum % of similarity between the pair of exons and 
their corresponding upstream and downstream exons. 
maxsize                          : ${params.maxsize}
% of maximum size ... 
cluster (txt files)              : ${params.cluster}
clusternum (number of clusters)  : ${params.clusternum}
output (output folder)           : ${params.output}
email for notification           : ${params.email}

"""

if (params.help) {
    log.info """This is the pipeline"""
    log.info """Please write some description here\n"""
    exit 1
}
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

clusterfile       = file(params.cluster)
outputQC          = "${params.output}/QC"
blosumfile        = file("${baseDir}/files/blosum62.txt")

if ( !blosumfile.exists() ) exit 1, "Missing blosum file: ${blosumfile}!"

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

genomes.join(annotations).into{pipe_data; data_to_annotation}



/*
 * Generate annotations
 */
process generate_annotations {
    tag { genomeid }
    publishDir "${params.output}/", mode: 'copy'	  

    input:
    set genomeid, file(genome), file(annotation) from data_to_annotation
    
    output:
    set val(genomeid), file (genomeid) into idfolders
	
	script:
	"""
	generate_annotations_lc.pl -GTF ${annotation} -G ${genome} -sp ${genomeid} 
	"""
}

/*
 * split cluster file
 */
process split_cluster_file_per_specie {
    tag { clusterfile }

    input:
    file(clusterfile)
    
    output:
    file("*.cls.tab") into cls_tab_files
    
	script:
	"""
   if [ `echo ${clusterfile} | grep ".gz"` ]; then 
       zcat ${clusterfile} > cluster_file			
       get_gcl_sp_pair.pl -f cluster_file
       rm cluster_file
    else
       get_gcl_sp_pair.pl -f ${clusterfile}
    fi
	"""
}


idfolders
  .toList().map{ [it, it] .combinations().findAll{ a, b -> a[0] < b[0]} }
  .flatMap()
  .map { ["${it[0][0]}-${it[1][0]}".toString(), it[0][1], it[1][1]] }
  .into{cluster_2_split; anno_2_score_ex_int; limorte}


/*
 * split clusters
 */
 
process split_clusters {
    tag { id_comb }

    input:
    file(cls_tab_files).collect()
    set id_comb, file(idfolder_A), file(idfolder_B) from cluster_2_split
        
    output:
    set file(idfolder_A), file(idfolder_B), file("${idfolder_A}_${idfolder_B}/*.cls.tab-part_*") into cls_files_2_align
    
	script:
	"""
		Prepare_and_Submit_Aln_sp_pair_V2.pl --sp1 ${idfolder_A} --sp2 ${idfolder_B} --expath ./ --project_dir ./ --N_split ${params.clusternum} --gene_cluster ${id_comb}.cls.tab
	"""
}

cls_files_2_align.transpose().set{cls_files_2_align_t}

/*
 * Align pairs
 */

process score_exons_introns_pair {
    tag { "${cls_part_file}" }
    label 'big_cpus'
    
    input:
    file(blosumfile)
    set file(sp1), file(sp2), file(cls_part_file) from cls_files_2_align_t

    output:
    set file(sp1), file(sp2), file("${sp1}-${sp2}_*") into aligned_subclusters_4_realign
    set val("${sp1}-${sp2}"), file("${sp1}-${sp2}_*") into aligned_subclusters_4_merge

	script:
    def cls_parts = "${cls_part_file}".split("_")
	"""
		Score_exons_introns_pair_sp.pl ${sp1} ${sp2} ${cls_part_file} \
${sp1}/${sp1}_annot_exons_prot_ids.txt ${sp2}/${sp2}_annot_exons_prot_ids.txt \
${sp1}/${sp1}_protein_ids_exons_pos.txt ${sp2}/${sp2}_protein_ids_exons_pos.txt \
${sp1}/${sp1}_protein_ids_intron_pos_CDS.txt ${sp2}/${sp2}_protein_ids_intron_pos_CDS.txt \
${sp1}/${sp1}.exint ${sp2}/${sp2}.exint ${cls_parts[1]} ${blosumfile} ${sp1}-${sp2}_${cls_parts[1]} ${task.cpus}
	"""
}

/*
 * Realign target exons
 */

process realign_exons {
    tag { "${aligned_output}" }
    label 'incr_time_cpus'
    
    input:
    file(blosumfile)
    set file(sp1), file(sp2), file(aligned_output) from aligned_subclusters_4_realign
    
    output:
    set val("${sp1}-${sp2}"), file("realigned_exons_${sp1}_${sp2}_*.txt") into realigned_exons_4_merge
   
	script:
    def cls_parts = "${aligned_output}".split("_")
	"""
		Realign_target_exons_by_part.pl ${sp1} ${sp2} ${aligned_output}/exons_to_realign_part_${cls_parts[1]}.txt \
${sp1}/${sp1}.exint ${sp2}/${sp2}.exint ${cls_parts[1]} realigned_exons_${sp1}_${sp2}_${cls_parts[1]}.txt \
${sp1}_${sp2} ${blosumfile} ${task.cpus}
	"""
}

aligned_subclusters_4_merge.groupTuple().join(realigned_exons_4_merge.groupTuple()).set{
	data_4_merge
}


/*
 * Join scores 
 */
 
process join_scores {
    tag { "${comp_id}" }
    
    input:
    set comp_id, file(folders), file(aligned_output) from data_4_merge
    
    output:
    set val(comp_id), file("${comp_id}/") into folder_jscores
   
	script:
	"""
	    mkdir ${comp_id}
	    ln -s \$PWD/realigned_exons_* ${comp_id}
	    ln -s \$PWD/${comp_id}_*/* ${comp_id}
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
    filter_exons_by_score.pl -b ${comp_id}/Best_score_hits_exons.txt -sps ${species[0]},${species[1]} -int ${params.intcons} -id ${params.idexons} -max_size ${params.maxsize}
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
    filter_exons_by_score.pl -b ${best_score} -sps ${species[0]},${species[1]} -int ${params.intcons} -id ${params.idexons} -max_size ${params.maxsize}
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
    file("Best_score_exon_hits_filtered_${params.maxsize}-${params.intcons}-${params.idexons}.tab") into filtered_all_scores

	script:
	"""
    cat best_score_* >> Best_score_exon_hits_filtered_${params.maxsize}-${params.intcons}-${params.idexons}.tab
    #for i in best_score_*; do cat \$i >> Best_score_exon_hits_filtered_${params.maxsize}-${params.intcons}-${params.idexons}.tab; done
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
	"""
    get_count_exons.pl ${scores} Exon_count_hits_by_sp.tab
    get_overlap_exons.pl -i Exon_count_hits_by_sp.tab -o Overlap_exons_by_sp.tab
    Filter_final_exons_pair.pl Overlap_exons_by_sp.tab ${scores} Best_score_exon_hits_pairs.txt
    """
}

/*
* Split the file of exon pairs 

process get_pre_cluster_exons {

    input:
    file (score_exon_hits_pairs)
    
    output:
    file("Best_score_exon_hits_filtered_${params.maxsize}-${params.intcons}-${params.idexons}.tab") into filtered_all_scores

	script:
	"""
    Get_Pre_cluster_exons.pl Consensus_cl_16sp-V5_2B.tab Best_score_exon_hits_liftover_pairs.txt 500 out.txt
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
    println "Pipeline BIOCORE@CRG YAMILE'S PIPELINE!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

