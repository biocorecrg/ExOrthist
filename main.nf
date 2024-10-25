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

include { PIPELINE_COMPLETION; PIPELINE_INITIALISATION } from "${LOCAL_SUBWORKFLOWS}/util.nf"

include { MAIN; PLOT } from "${WORKFLOWS}/plot.nf"


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
        MAIN()
    }

    PIPELINE_COMPLETION(
        params.wf,
        params.email,
        params.hook_url
    )

}
