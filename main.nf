#!/usr/bin/env nextflow

/*
 * Copyright (c) 2019-2024, Centre for Genomic Regulation (CRG)
 *
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

// if( !workflow.resume ) {
//     println "Removing the output folder"
// 	  new File("${params.output}").delete()
// }

LOCAL_SUBWORKFLOWS='./subworkflows/local/exorthist'
WORKFLOWS='./workflows/'

include { ALIGN } from "${LOCAL_SUBWORKFLOWS}/align.nf"
include { CLUSTER } from "${LOCAL_SUBWORKFLOWS}/cluster.nf"

include { PIPELINE_COMPLETION; PIPELINE_INITIALISATION } from "${LOCAL_SUBWORKFLOWS}/util.nf"

include { PREPARE } from "${LOCAL_SUBWORKFLOWS}/prepare.nf"
include { SCORE } from "${LOCAL_SUBWORKFLOWS}/score.nf"

include { PLOT } from "${WORKFLOWS}/plot.nf"

workflow {

    PIPELINE_INITIALISATION(
        params,
        args
    )

    if (params.wf == "plot" ) {

        PLOT(
            params.output,
            params.geneID,
            params.relevant_exs,
            params.ordered_species,
            params.isoformID,
            params.sub_orthologs
        )

    } else {

        PREPARE(
            params.evodists,
            params.cluster,
            params.genomes,
            params.annotations,
            params.long_dist,
            params.medium_dist,
            params.short_dist,
            params.extraexons,
            params.alignmentnum
        )

        ALIGN(
            "${projectDir}/files/blosum62.txt",
            PREPARE.out.alignment_input,
            PREPARE.out.clusters_split_ch,
            params.long_dist,
            params.medium_dist,
            params.short_dist,
            params.alignmentnum,
            params.prevaln
        )

        SCORE(
            ALIGN.out.folder_jscores,
            PREPARE.out.clusters_split_ch,
            PREPARE.out.dist_ranges_ch,
            params.bonafide_pairs,
            params.long_dist,
            params.medium_dist,
            params.short_dist
        )

        CLUSTER(
            SCORE.out.score_exon_hits_pairs,
            PREPARE.out.clusters_split_ch,
            params.cluster,
            params.orthopairs,
            params.orthogroupnum
        )
    }

    PIPELINE_COMPLETION(
        params.wf,
        params.email,
        params.hook_url
    )

}
