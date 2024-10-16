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

include { CHECK_INPUT } from "${LOCAL_MODULES}/check_input.nf"
include { CLUSTER_EXS } from "${LOCAL_MODULES}/cluster_exons.nf"
include { COLLAPSE_OVERLAPPING_MATCHES } from "${LOCAL_MODULES}/collapse_matches.nf"
include { FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE } from "${LOCAL_MODULES}/filter_matches.nf"
include { FORMAT_EX_CLUSTERS_INPUT } from "${LOCAL_MODULES}/format_input.nf"
include { FORMAT_EX_CLUSTERS_OUTPUT } from "${LOCAL_MODULES}/format_output.nf"
include { GENERATE_ANNOTATIONS } from "${LOCAL_MODULES}/generate_annotations.nf"
include { JOIN_FILTERED_EX_MATCHES } from "${LOCAL_MODULES}/join_matches.nf"
include { MERGE_PROT_EX_INT_ALN_INFO } from "${LOCAL_MODULES}/merge_aligns.nf"
include { PARSE_IPA_PROT_ALN } from "${LOCAL_MODULES}/align_pairs.nf"
include { REALIGN_EX_PAIRS } from "${LOCAL_MODULES}/realign_pairs.nf"
include { RECLUSTER_EXS_BY_SPECIES_PAIR } from "${LOCAL_MODULES}/recluster_exs.nf"
include { RECLUSTER_GENES_BY_SPECIES_PAIR } from "${LOCAL_MODULES}/recluster_genes.nf"
include { SCORE_EX_MATCHES } from "${LOCAL_MODULES}/score_matches.nf"
include { SPLIT_CLUSTERS_IN_CHUNKS } from "${LOCAL_MODULES}/split_clusters_chunks.nf"
include { SPLIT_CLUSTERS_BY_SPECIES_PAIRS } from "${LOCAL_MODULES}/split_clusters_species.nf"
include { SPLIT_EX_PAIRS_TO_REALIGN } from "${LOCAL_MODULES}/split_pairs.nf"

include { PREPARE_INPUT } from "${LOCAL_SUBWORKFLOWS}/prepare_input.nf"

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

    PREPARE_INPUT(
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

    // the last argument is the protein similarity alignment.
    // if a prevaln folder is provided, the protein alignments present in each species pair subfolder will not be repeated.
    PARSE_IPA_PROT_ALN(blosumfile, PREPARE_INPUT.out.alignment_input)

    // Collapse EXs_to_split in batches of 500 files
    EXs_to_split = PARSE_IPA_PROT_ALN.out.EXs_to_split
    EXs_to_split_batches = EXs_to_split.toSortedList().flatten().buffer(size : 500, remainder: true)
    // Split exons pairs to realign
    SPLIT_EX_PAIRS_TO_REALIGN(EXs_to_split_batches)
    EXs_to_realign_batches = SPLIT_EX_PAIRS_TO_REALIGN.out.EXs_to_realign_batches
    // Flatten the results from the previous batch run and combine with sp1 and sp2 information, using sp1-sp2 as key.
    EXs_to_realign = EXs_to_realign_batches.flatten().map{[it.getName().toString().split("_")[0],it]}.groupTuple().join(PREPARE_INPUT.out.clusters_split_ch).transpose()
    //  Realign exons pairs (with multiple hits)
    REALIGN_EX_PAIRS(blosumfile, EXs_to_realign)
    // Combine all the aln_info with the realigned_exon_info for each species pair
    aligned_subclusters_4_splitting = PARSE_IPA_PROT_ALN.out.aligned_subclusters_4_splitting
    realigned_exons_4_merge = REALIGN_EX_PAIRS.out.realigned_exons_4_merge
    data_4_merge = aligned_subclusters_4_splitting.groupTuple().join(realigned_exons_4_merge.groupTuple())
    // Merge alignments information
    MERGE_PROT_EX_INT_ALN_INFO(data_4_merge)
    folder_jscores = MERGE_PROT_EX_INT_ALN_INFO.out.folder_jscores
    data_to_score = folder_jscores.join(PREPARE_INPUT.out.clusters_split_ch).map{ [it[0], it[1..-1] ]}
    // Score EX matches from aln info
    SCORE_EX_MATCHES(data_to_score)
    // Filter the best matches above score cutoffs by target gene.
    all_scores_to_filt_ch = SCORE_EX_MATCHES.out.all_scores_to_filt
    FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE(all_scores_to_filt_ch.join(PREPARE_INPUT.out.dist_ranges_ch))
    //  Join filtered scored EX matches
    filterscore_per_joining_ch = FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE.out.filterscore_per_joining
    JOIN_FILTERED_EX_MATCHES(filterscore_per_joining_ch.collect())
    //  Removing matches from overlapping exons

    filtered_all_scores = JOIN_FILTERED_EX_MATCHES.out.filtered_all_scores
    // Sic: https://nextflow-io.github.io/patterns/optional-input/
    if ( params.bonafide_pairs ) {
        COLLAPSE_OVERLAPPING_MATCHES(filtered_all_scores, params.bonafide_pairs)
    } else {
        COLLAPSE_OVERLAPPING_MATCHES(filtered_all_scores, file("/path/to/NO_FILE"))
    }
    score_exon_hits_pairs = COLLAPSE_OVERLAPPING_MATCHES.out.score_exon_hits_pairs
    FORMAT_EX_CLUSTERS_INPUT(score_exon_hits_pairs, file(params.cluster))

    // Split the file of exon pairs
    // Unclustered are the exons ending up in single-exon clusters
    cluster_parts = FORMAT_EX_CLUSTERS_INPUT.out.cluster_parts.flatten()
    CLUSTER_EXS(cluster_parts)
    FORMAT_EX_CLUSTERS_OUTPUT(CLUSTER_EXS.out.ex_clusters.collect(), CLUSTER_EXS.out.unclustered_exs.collect())

    // Re-clustering of genes
    RECLUSTER_GENES_BY_SPECIES_PAIR(
        PREPARE_INPUT.out.clusters_split_ch,
        clusterfile_ch,
        orthopairs_ch
    )

    // Re-clustering of exons
    RECLUSTER_EXS_BY_SPECIES_PAIR(
        RECLUSTER_GENES_BY_SPECIES_PAIR.out.recl_genes_for_rec_exons,
        FORMAT_EX_CLUSTERS_OUTPUT.out.exon_cluster_for_reclustering,
        COLLAPSE_OVERLAPPING_MATCHES.out.score_exon_hits_pairs,
        orthopairs_ch
    )

}

workflow.onComplete {
    println "--- Pipeline BIOCORE@CRG ExOrthist ---"
    println "Started at  $workflow.start"
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

