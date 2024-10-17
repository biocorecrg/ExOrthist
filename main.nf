#!/usr/bin/env nextflow


/*
 * Copyright (c) 2019-2024, Centre for Genomic Regulation (CRG)
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
 Federica Mantica <federica.mantica93@gmail.com>
 Toni Hermoso Pulido <toni.hermoso@crg.eu>
===========================================================
*/

nextflow.enable.dsl=2

version = '2.0.0'

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
pairwise evo distances           : ${params.evodists}
long distance parameters         : ${params.long_dist}
medium distance parameters       : ${params.medium_dist}
short distance parameters        : ${params.short_dist}
pre-computed alignments		 : ${params.prevaln}
alignment number		 : ${params.alignmentnum}
orthogroup number		 : ${params.orthogroupnum}
extraexons (e.g. from VastDB)    : ${params.extraexons}
bona fide orthologous exon pairs  : ${params.bonafide_pairs}
orthopairs                       : ${params.orthopairs}
output (output folder)           : ${params.output}
email for notification           : ${params.email}

INFORMATION ABOUT OPTIONS:
The long, medium, short distance cut-offs are in the format: "int_num;ex_seq;ex_len;prot_sim".
Only exon matches respecting all cut-offs are considered homologous.
- int_num (0,1,2): Number of surrounding intron positions required to be conserved.
- ex_seq (from 0 to 1): Minimum sequence similarity % between a
     pair of homologous exons and their corresponding upstream and
     downstream exons.
- ex_len (from 0 to 1): Maximum size difference between two homologous exons
     (as a fraction of either exon).
- prot_sim (from 0 to 1): Minimum sequence similarity over the entire pairwise alignment
     for a pair of protein isoforms to be considered for comparison.

See online README at https://github.com/biocorecrg/ExOrthist for further information about the options.
"""

if (params.help) {
    log.info """ExOrthist v2.0.0"""
    log.info """ExOrthist is a Nextflow-based pipeline to obtain groups of exon orthologous at all evolutionary timescales.\n"""
    exit 1
}
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

if( !workflow.resume ) {
    println "Removing the output folder"
	  new File("${params.output}").delete()
}

// clusterfile       = file(params.cluster)
outputQC          = "${params.output}/QC"
blosumfile        = file("${baseDir}/files/blosum62.txt")
// evodisfile	  = file(params.evodists)

if ( !blosumfile.exists() ) exit 1, "Missing blosum file: ${blosumfile}!"
// if ( !clusterfile.exists() ) exit 1, "Missing clusterfile file: ${clusterfile}!"
// if ( !evodisfile.exists() ) exit 1, "Missing evodists file: ${evodisfile}!"
//

LOCAL_MODULES='./modules/local/exorthist'
LOCAL_SUBWORKFLOWS='./subworkflows/local/exorthist'

include { ALIGN } from "${LOCAL_SUBWORKFLOWS}/align.nf"
include { CLUSTER } from "${LOCAL_SUBWORKFLOWS}/cluster.nf"
include { PREPARE } from "${LOCAL_SUBWORKFLOWS}/prepare.nf"
include { SCORE } from "${LOCAL_SUBWORKFLOWS}/score.nf"

/*
 * Validate input and print log file
 */
//Prepare input channels
workflow {
    gtfs = Channel.fromPath(params.annotations).collect()
    fastas = Channel.fromPath(params.genomes).collect()

    // TODO: Review this in an easier way
    gtfs_suffix = Channel.fromFilePairs(params.annotations, size: 1).flatten().collate(2).map{[it[1].getName().toString().split(it[0].toString())[1]]}.unique().flatten()
    fastas_suffix = Channel.fromFilePairs(params.genomes, size: 1).flatten().collate(2).map{[it[1].getName().toString().split(it[0].toString())[1]]}.unique().flatten()

    // Channels for sequences of data
    genomes = Channel
        .fromFilePairs(params.genomes, size: 1)
        .ifEmpty { error "Cannot find any genome matching: ${params.genomes}" }

    annotations = Channel
        .fromFilePairs(params.annotations, size: 1)
        .ifEmpty { error "Cannot find any annotation matching: ${params.annotations}" }

    extraexons = params.extraexons ?
        Channel.fromFilePairs(params.extraexons, checkIfExists: true, size: 1)
        .ifEmpty { error "Extra exons not found" } :
        Channel.empty()

    // We join channels. If no extraexons, then it's empty, so no problem
    data_to_annotation_raw = genomes.join(annotations)
    data_to_annotation = data_to_annotation_raw.join(extraexons, remainder: true)

    evodists_ch = Channel.fromPath(params.evodists)
    clusterfile_ch = Channel.fromPath(params.cluster)
    if ( params.orthopairs ) {
        orthopairs_ch = file(params.orthopairs)
    } else {
        orthopairs_ch = file("/path/to/NO_FILE")
    }

    PREPARE(
        evodists_ch,
        clusterfile_ch,
        gtfs,
        fastas,
        gtfs_suffix,
        fastas_suffix,
        params.long_dist,
        params.medium_dist,
        params.short_dist,
        data_to_annotation,
        params.extraexons
    )

    ALIGN(blosumfile, PREPARE.out.alignment_input, PREPARE.out.clusters_split_ch)
    SCORE(ALIGN.out.folder_jscores, PREPARE.out.clusters_split_ch, PREPARE.out.dist_ranges_ch, params.bonafide_pairs)
    CLUSTER(SCORE.out.score_exon_hits_pairs, file(params.cluster), PREPARE.out.clusters_split_ch, clusterfile_ch, orthopairs_ch)
}

workflow.onComplete {
    println "--- Pipeline BIOCORE@CRG ExOrthist ---"
    println "Started at  $workflow.start"
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

