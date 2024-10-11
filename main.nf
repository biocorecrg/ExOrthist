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

include { CHECK_INPUT } from "${LOCAL_MODULES}/check_input.nf"
include { CLUSTER_EXS } from "${LOCAL_MODULES}/cluster_exons.nf"
include { COLLAPSE_OVERLAPPING_MATCHES } from "${LOCAL_MODULES}/collapse_matches.nf"
include { FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE } from "${LOCAL_MODULES}/filter_matches.nf"
include { FORMAT_EX_CLUSTERS_INPUT } from "${LOCAL_MODULES}/format_clusters.nf"
include { GENERATE_ANNOTATIONS } from "${LOCAL_MODULES}/generate_annotations.nf"
include { JOIN_FILTERED_EX_MATCHES } from "${LOCAL_MODULES}/join_matches.nf"
include { MERGE_PROT_EX_INT_ALN_INFO } from "${LOCAL_MODULES}/merge_aligns.nf"
include { PARSE_IPA_PROT_ALN } from "${LOCAL_MODULES}/align_pairs.nf"
include { REALIGN_EX_PAIRS } from "${LOCAL_MODULES}/realign_pairs.nf"
include { SCORE_EX_MATCHES } from "${LOCAL_MODULES}/score_matches.nf"
include { SPLIT_CLUSTERS_IN_CHUNKS } from "${LOCAL_MODULES}/split_clusters_chunks.nf"
include { SPLIT_CLUSTERS_BY_SPECIES_PAIRS } from "${LOCAL_MODULES}/split_clusters_species.nf"
include { SPLIT_EX_PAIRS_TO_REALIGN } from "${LOCAL_MODULES}/split_pairs.nf"
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
    pipe_data = data_to_annotation_raw
    data_to_annotation = data_to_annotation_raw.join(extraexons, remainder: true)

    // Print contents of each channel
    gtfs.view { "GTF file: $it" }
    fastas.view { "FASTA file: $it" }
    gtfs_suffix.view { "GTF suffix: $it" }
    fastas_suffix.view { "FASTA suffix: $it" }

    genomes.view { "Genome file: $it" }
    annotations.view { "Genome file: $it" }
    extraexons.view { "Extra: $it" }
    data_to_annotation_raw.view { "DAR: $it"}
    data_to_annotation.view { "DA: $it"}
    pipe_data.view { "PD: $it"}

    evodists_ch = Channel.fromPath(params.evodists)
    clusterfile_ch = Channel.fromPath(params.cluster)

    CHECK_INPUT(
        evodists_ch,
        clusterfile_ch,
        gtfs,
        fastas,
        gtfs_suffix,
        fastas_suffix,
        params.long_dist,
        params.medium_dist,
        params.short_dist
    )

    // Sic: https://nextflow-io.github.io/patterns/optional-input/
    if ( params.extraexons ) {
         GENERATE_ANNOTATIONS(data_to_annotation, extraexons)
    } else {
         GENERATE_ANNOTATIONS(data_to_annotation, file("/path/to/NO_FILE"))
    }

    clusters_split_ch = GENERATE_ANNOTATIONS.out.idfolders.toList().map{ [it, it].combinations().findAll{ a, b -> a[0] < b[0]} }
        .flatMap()
        .map { ["${it[0][0]}-${it[1][0]}".toString(), it[0][1], it[1][1]] }


    // Copy the gene cluster file to output to use for the exint_plotter and compare_exon_sets modules
    SPLIT_CLUSTERS_BY_SPECIES_PAIRS(clusterfile_ch)

    // Split clusters
    cls_tab_files_ch = SPLIT_CLUSTERS_BY_SPECIES_PAIRS.out.cls_tab_files
    SPLIT_CLUSTERS_IN_CHUNKS(cls_tab_files_ch.collect(), clusters_split_ch)

    cls_files_2_align = SPLIT_CLUSTERS_IN_CHUNKS.out.cls_files_2_align
    cls_files_2_align_t = cls_files_2_align.transpose().map{[it[0].getFileName().toString()+"-"+it[1].getFileName().toString(), it[0], it[1], it[2]]}

    //Create a channel for the evo distances
    sp1_sp2_dist = Channel
     .fromPath("${params.evodists}")
     .splitText()
     .map{"${it}".trim().split("\t")}.map{[it[0]+"-"+it[1], it[2]]}

    sp2_sp1_dist = Channel
     .fromPath("${params.evodists}")
     .splitText()
     .map{"${it}".trim().split("\t")}.map{[it[1]+"-"+it[0], it[2]]}

    species_pairs_dist = sp1_sp2_dist.concat(sp2_sp1_dist)
    //Only the species pairs with a common index will be kept
    dist_ranges_ch = clusters_split_ch.join(species_pairs_dist).map{[it[0], it[3]]}
    alignment_input = cls_files_2_align_t.groupTuple().join(dist_ranges_ch).transpose()

    //the last argument is the protein similarity alignment.
    //if a prevaln folder is provided, the protein alignments present in each species pair subfolder will not be repeated.
    PARSE_IPA_PROT_ALN(blosumfile, alignment_input)

    // Collapse EXs_to_split in batches of 500 files
    EXs_to_split = PARSE_IPA_PROT_ALN.out.EXs_to_split
    EXs_to_split_batches = EXs_to_split.toSortedList().flatten().buffer(size : 500, remainder: true)
    // Split exons pairs to realign
    SPLIT_EX_PAIRS_TO_REALIGN(EXs_to_split_batches)
    EXs_to_realign_batches = SPLIT_EX_PAIRS_TO_REALIGN.out.EXs_to_realign_batches
    // Flatten the results from the previous batch run and combine with sp1 and sp2 information, using sp1-sp2 as key.
    EXs_to_realign = EXs_to_realign_batches.flatten().map{[it.getName().toString().split("_")[0],it]}.groupTuple().join(clusters_split_ch).transpose()
    //  Realign exons pairs (with multiple hits)
    REALIGN_EX_PAIRS(blosumfile, EXs_to_realign)
    //Combine all the aln_info with the realigned_exon_info for each species pair
    aligned_subclusters_4_splitting = PARSE_IPA_PROT_ALN.out.aligned_subclusters_4_splitting
    realigned_exons_4_merge = REALIGN_EX_PAIRS.out.realigned_exons_4_merge
    data_4_merge = aligned_subclusters_4_splitting.groupTuple().join(realigned_exons_4_merge.groupTuple())
    // Merge alignments information
    MERGE_PROT_EX_INT_ALN_INFO(data_4_merge)
    folder_jscores = MERGE_PROT_EX_INT_ALN_INFO.out.folder_jscores
    data_to_score = folder_jscores.join(clusters_split_ch).map{ [it[0], it[1..-1] ]}
    // Score EX matches from aln info
    SCORE_EX_MATCHES(data_to_score)
    // Filter the best matches above score cutoffs by target gene.
    all_scores_to_filt_ch = SCORE_EX_MATCHES.out.all_scores_to_filt
    FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE(all_scores_to_filt_ch.join(dist_ranges_ch))
    //  Join filtered scored EX matches
    filterscore_per_joining_ch = FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE.out.filterscore_per_joining
    JOIN_FILTERED_EX_MATCHES(filterscore_per_joining_ch)
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
    cluster_parts = FORMAT_EX_CLUSTERS_INPUT.out.cluster_parts
    CLUSTER_EXS(cluster_parts)

    // Review outputs below
    CHECK_INPUT.out.run_info.view()
    GENERATE_ANNOTATIONS.out.idfolders.view { "ANN: $it" }
    SPLIT_CLUSTERS_BY_SPECIES_PAIRS.out.cls_tab_files.view()
    SPLIT_CLUSTERS_BY_SPECIES_PAIRS.out.gene_cluster_file.view()
    clusters_split_ch.view { "CL: $it" }
    dist_ranges_ch.view { "DR: $it" }
    alignment_input.view { "AL: $it" }
    PARSE_IPA_PROT_ALN.out.aligned_subclusters_4_splitting.view { "SC: $it" }
    PARSE_IPA_PROT_ALN.out.EXs_to_split.view { "EX: $it" }
    EXs_to_realign.view { "EXR: $it" }
    REALIGN_EX_PAIRS.out.realigned_exons_4_merge.view{ "RER: $it" }
    MERGE_PROT_EX_INT_ALN_INFO.out.folder_jscores.view()
    MERGE_PROT_EX_INT_ALN_INFO.out.aln_features.view()
    MERGE_PROT_EX_INT_ALN_INFO.out.exint_aln.view()
    SCORE_EX_MATCHES.out.all_features.view()
    SCORE_EX_MATCHES.out.all_scores_to_filt.view()
    FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE.out.filterscore_per_joining.view()
    FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE.out.best_scored_matches.view()
    JOIN_FILTERED_EX_MATCHES.out.filtered_all_scores.view()
    COLLAPSE_OVERLAPPING_MATCHES.out.score_exon_hits_pairs.view()
    COLLAPSE_OVERLAPPING_MATCHES.out.overlapping_exs.view()
    FORMAT_EX_CLUSTERS_INPUT.out.cluster_parts.view()
    CLUSTER_EXS.out.unclustered_exs.view()
    CLUSTER_EXS.out.ex_clusters.view()
}

//
// /*
// * Final exon clusters
// */
// process format_EX_clusters_output {
//     publishDir "${params.output}/", mode: 'copy'
//
//     input:
//     file("*") from ex_clusters.collect()
//     file("*") from unclustered_exs.collect()
//
//     output:
//     file("EX_clusters.tab") into exon_cluster_for_reclustering
//     file("EX_clusters_info.tab.gz")
//     file("unclustered_EXs.txt")
//
// 	script:
// 	"""
// 	D3_format_EX_clusters_output.pl
// 	"""
// }
//
// /*
// * Re-clustering of genes
// */
//
// process recluster_genes_by_species_pair {
//     publishDir "${params.output}/reclustering", mode: 'copy'
//     tag { "${combid}" }
//
//     when:
//     params.orthopairs != ""
//
//     input:
//     set combid, file(folderA), file(folderB) from species_to_recluster_genes
//     file(clusterfile)
//     //file(orthopairs) from orthopairs_file
//
//     output:
//     set combid, file("reclustered_genes_*.tab") into recl_genes_for_rec_exons
//
// 	script:
// 	def species = "${combid}".split("-")
// 	def orthopairs = file("${params.orthopairs}")
// 	"""
//  	D3.1_recluster_genes_by_species_pair.py -og ${clusterfile} -op ${orthopairs} --species1 ${species[0]} --species2 ${species[1]} -out reclustered_genes_${combid}.tab
//     	"""
// }
//
//
// /*
// * Re-clustering of exons
// */
//
// process recluster_EXs_by_species_pair {
//     publishDir "${params.output}/reclustering", mode: 'copy'
//     tag { "${combid}" }
//
// 	when:
// 	params.orthopairs != ""
//
//     input:
//     set combid, file(recl_genes) from recl_genes_for_rec_exons
//     file(exon_clusters) from exon_cluster_for_reclustering
//     file(exon_pairs) from exon_pairs_for_reclustering
//
//     output:
//     file("reclustered_EXs_*.tab")
//
// 	script:
// 	def combid1 = combid.replace("-", "_")
// 	def species = "${combid1}".split("-")
// 	"""
// 	D3.2_recluster_EXs_by_species_pair.py -ep ${exon_pairs} -rg ${recl_genes} -ec ${exon_clusters} -sp1 ${species[0]} -sp2 ${species[1]} -out reclustered_EXs_${combid1}.tab
//     	"""
// }
//
// /*
//  * Mail notification
//  */
//
// if (params.email == "yourmail@yourdomain" || params.email == "") {
//     log.info 'Skipping the email\n'
// }
// else {
//     log.info "Sending the email to ${params.email}\n"
//
//     workflow.onComplete {
//
//     def msg = """\
//         Pipeline execution summary
//         ---------------------------
//         Completed at: ${workflow.complete}
//         Duration    : ${workflow.duration}
//         Success     : ${workflow.success}
//         workDir     : ${workflow.workDir}
//         exit status : ${workflow.exitStatus}
//         Error report: ${workflow.errorReport ?: '-'}
//         """
//         .stripIndent()
//
//         sendMail(to: params.email, subject: "VectorQC execution", body: msg,  attach: "${outputMultiQC}/multiqc_report.html")
//     }
// }
//
// workflow.onComplete {
//     println "--- Pipeline BIOCORE@CRG ExOrthist ---"
//     println "Started at  $workflow.start"
//     println "Finished at $workflow.complete"
//     println "Time elapsed: $workflow.duration"
//     println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
// }
//
