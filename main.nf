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
                                                                                       
====================================================
annotations (folder with GTF files)                : ${params.annotations}
genomes (folder with fasta files)                  : ${params.genomes}
output (output folder)           : ${params.output}
email for notification           : ${params.email}
"""

if (params.help) {
    log.info """This is the pipeline"""
    log.info """Please write some description here\n"""
    exit 1
}
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

annofolder        = getFolderName(params.annotations)
seqfolder         = getFolderName(params.genomes)
outputQC          = "${params.output}/QC"


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
pipe_data.println()

/*
 * Generate annotations
 */
process generate_annotations {
    publishDir "${params.output}/", mode: 'copy'	  

    input:
    set id, file(genome), file(annotation) from data_to_annotation
    
    output:
    file("${id}") into idfolder
    
	script:
	"""
	generate_annotations.pl -GTF ./ -G ./ -sp ${id} 
	"""

}

/*
* functions
*/

def getFolderName(sample) {
   folder_info = sample.toString().tokenize("/")
   return folder_info[-2]
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

