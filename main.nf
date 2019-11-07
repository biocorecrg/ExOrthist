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
annotations (GTF files)      : ${params.annotations}
genomes (fasta files)     	 : ${params.genomes}
cluster file (txt files)     : ${params.cluster}
output (output folder)       : ${params.output}
email for notification       : ${params.email}
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
	generate_annotations_lc.pl -GTF ${annotation} -G ${genomeid} -sp ${genomeid} 
	"""
}

/*
 * split cluster file
 */
process split_cluster_file {
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
  .toList().map{ [it, it].combinations().findAll{ a, b -> a[1] > b[1]} }
  .flatMap()
  .map { ["${it[0][0]}-${it[1][0]}", it[0][1], it[1][1]] }
  .set{cluster_2_split}

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
		Prepare_and_Submit_Aln_sp_pair.pl --sp1 ${idfolder_A} --sp2 ${idfolder_B} --expath ./ --project_dir ./ --gene_cluster ${id_comb}.cls.tab
	"""
}

cls_files_2_align.transpose().set{cls_files_2_align_t}

/*
 * Align pairs
 */
 
process align_pairs {
    tag { "${cls_part_file}" }
    label 'big_cpus'
    
    input:
    file(blosumfile)
    set file(sp1), file(sp2), file(cls_part_file) from cls_files_2_align_t

    output:
    set file(sp1), file(sp2), file("${sp1}-${sp2}_*") into aligned_subclusters

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
    label 'big_cpus'
    
    input:
    file(blosumfile)
    set file(sp1), file(sp2), file(aligned_output) from aligned_subclusters
    
    output:
    file("realigned_exons_${sp1}_${sp2}_*.txt") into realigned_exons
   
	script:
    def cls_parts = "${aligned_output}".split("_")
	"""
		Realign_target_exons_by_part.pl ${sp1} ${sp2} ${aligned_output}/exons_to_realign_part_${cls_parts[1]}.txt \
${sp1}/${sp1}.exint ${sp2}/${sp2}.exint ${cls_parts[1]} realigned_exons_${sp1}_${sp2}_${cls_parts[1]}.txt \
${sp1}_${sp2} ${blosumfile} ${task.cpus}
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
    println "Pipeline BIOCORE@CRG completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

