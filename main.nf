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
include { GENERATE_ANNOTATIONS } from "${LOCAL_MODULES}/generate_annotations.nf"
include { PARSE_IPA_PROT_ALN } from "${LOCAL_MODULES}/align_pairs.nf"
include { REALIGN_EX_PAIRS } from "${LOCAL_MODULES}/realign_pairs.nf"
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
}


//
// //Combine all the aln_info with the realigned_exon_info for each species pair
// aligned_subclusters_4_splitting.groupTuple().join(realigned_exons_4_merge.groupTuple()).set{data_4_merge}
//
//
// /*
//  * Merge alignments information
//  */
//
// process merge_PROT_EX_INT_aln_info {
//     tag { "${comp_id}" }
//     label 'incr_time_cpus'
//
//     stageInMode = 'copy'
//     //this matches all_PROT_aln_features.txt, all_EX_aln_features.txt, all_INT_aln_features.txt, Exint_Alignments.aln.gz
//     publishDir "${params.output}", mode: "copy", pattern: "${comp_id}/all_*_aln_features.txt"
//     publishDir "${params.output}", mode: "copy", pattern: "${comp_id}/EXINT_aln.gz"
//
//     input:
//     set comp_id, file("FOLDERS_*"), file("*") from data_4_merge //05/03/21
//
//     output:
//     set val(comp_id), file("${comp_id}/") into folder_jscores
//     file("${comp_id}/all_*_aln_features.txt")
//     file("${comp_id}/EXINT_aln.gz")
// 	script:
// 	"""
// 	    mkdir ${comp_id}
// 	    mv FOLDERS_*/* ${comp_id}/
// 	    mv realigned_* ${comp_id}/
//
//     	B4_merge_PROT_EX_INT_aln_info.pl ${comp_id}
// 	"""
// }
//
// folder_jscores.join(anno_2_score_ex_int).map{
//    [it[0], it[1..-1] ]
// }.set{data_to_score}
//
//
// /*
//  * Score EX matches from aln info
//  */
//
// process score_EX_matches {
//     tag { "${comp_id}" }
//     label('big_mem_retry')
//     //I need to modify the name so that it has the species pair in the output.
//     storeDir "${params.output}/${comp_id}"
//
//     input:
//     set val(comp_id), file("*") from data_to_score
//
//     output:
//     file("all_PROT_EX_INT_aln_features_*")
//     set val(comp_id), file("all_scored_EX_matches.txt") into all_scores_to_filt
//
// 	script:
//     def species = comp_id.split("-")
// 	"""
//     B5_format_aln_info_by_best_isoform_match.pl ${species[0]} ${species[1]} \
//     ${comp_id}/all_PROT_aln_features.txt ${comp_id}/all_EX_aln_features.txt ${comp_id}/all_INT_aln_features.txt \
//     ${species[0]}/${species[0]}.exint ${species[1]}/${species[1]}.exint \
//     ${species[0]}/${species[0]}_protein_ids_intron_pos_CDS.txt ${species[1]}/${species[1]}_protein_ids_intron_pos_CDS.txt \
//     all_PROT_EX_INT_aln_features_${comp_id}.txt;
//     C1_score_EX_matches.pl all_PROT_EX_INT_aln_features_${comp_id}.txt .
//     """
// }
//
//
// /*
//  * Filter the best matches above score cutoffs by target gene.
//  */
//
// process filter_and_select_best_EX_matches_by_targetgene {
//     tag { "${comp_id}" }
//     publishDir "${params.output}", mode: "copy", pattern: "best_scored_EX_matches_by_targetgene.txt", saveAs: { filename -> "${comp_id}/foo_$filename" }
//
//     input:
//     set val(comp_id), file(all_scores), val(dist_range) from all_scores_to_filt.join(dist_ranges_ch)
//
//     output:
//     file("*.tab") into filterscore_per_joining
//     file("best_scored_EX_matches_by_targetgene.txt")
//
//     script:
//     def species = comp_id.split("-")
//     if (dist_range == "long")
// 	dist_range_par = "${params.long_dist}".split(",")
//     if (dist_range == "medium")
// 	dist_range_par = "${params.medium_dist}".split(",")
//     if (dist_range == "short")
// 	dist_range_par = "${params.short_dist}".split(",")
//     """
//     C2_filter_and_select_best_EX_matches_by_targetgene.pl -b ${all_scores} -sps ${species[0]},${species[1]} -int ${dist_range_par[0]} -id ${dist_range_par[1]} -max_size ${dist_range_par[2]}
//     """
// }
//
// /*
//  * join filtered scored EX matches
//  */
// process join_filtered_EX_matches {
//     publishDir "${params.output}/", mode: 'copy'
//
//     input:
//     file ("filtered_best_scored-*") from filterscore_per_joining.collect()
//
//     output:
//     file("filtered_best_scored_EX_matches_by_targetgene.tab") into filtered_all_scores
//
// 	script:
// 	"""
//     echo "GeneID_sp1\tExon_coords_sp1\tGeneID_sp2\tExon_coords_sp2\tSp1\tSp2" > filtered_best_scored_EX_matches_by_targetgene.tab;
//     for file in \$(ls filtered_best_scored-*); do cat \$file | tail -n+2 >> filtered_best_scored_EX_matches_by_targetgene.tab; done
//     """
// }
//
// /*
//  * Removing matches from overlapping exons
//  */
// process collapse_overlapping_matches {
//     publishDir "${params.output}/", mode: "copy"
//     input:
//     file(scores) from filtered_all_scores
//
//     output:
//     file("filtered_best_scored_EX_matches_by_targetgene-NoOverlap.tab") into (score_exon_hits_pairs, exon_pairs_for_reclustering)
//     file("overlapping_EXs_by_species.tab")
//
// 	script:
// 	bonafide = ""
// 	if (params.bonafide_pairs != "") {
// 		bonafide = file(params.bonafide_pairs)
// 		if ( !bonafide.exists() ) exit 1, "Missing file with bonafide pairs: ${bonafide}!"
//
// 	}
// 	"""
//     C3_count_matches_by_EX.pl ${scores} EX_matches_count_by_species.tab ${bonafide}
//     C4_get_overlapping_EXs.pl -i EX_matches_count_by_species.tab -o overlapping_EXs_by_species.tab
//     C5_collapse_overlapping_matches.pl overlapping_EXs_by_species.tab ${scores} filtered_best_scored_EX_matches_by_targetgene-NoOverlap.tab ${bonafide}
//     """
// }
//
// /*
// * Split the file of exon pairs
// */
// clusterfile = file(params.cluster)
// process format_EX_clusters_input {
//     input:
//     file(score_exon_hits_pairs)
//     file(clusterfile)
//
//     output:
//     file("PART_*-cluster_input.tab") into cluster_parts
//
// 	script:
// 	"""
//    if [ `echo ${clusterfile} | grep ".gz"` ]; then
//        zcat ${clusterfile} > cluster_file
//        D1_format_EX_clusters_input.pl cluster_file ${score_exon_hits_pairs} ${params.orthogroupnum}
//        rm cluster_file
//     else
//        D1_format_EX_clusters_input.pl ${clusterfile} ${score_exon_hits_pairs} ${params.orthogroupnum}
//     fi
// 	"""
// }
//
// /*
// * Split the file of exon pairs
// */
// //Unclustered are the exons ending up in single-exon clusters
// process cluster_EXs {
//
//     input:
//     file(cluster_part) from cluster_parts.flatten()
//
//     output:
//     file("EXs_${cluster_part}") into ex_clusters
//     file("unclustered_EXs_${cluster_part}") into unclustered_exs
//
// 	script:
// 	"""
//     D2_cluster_EXs.R ${cluster_part} EXs_${cluster_part} unclustered_EXs_${cluster_part}
//     """
// }
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
