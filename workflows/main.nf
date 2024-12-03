LOCAL_SUBWORKFLOWS='../subworkflows/local/exorthist'

include { ALIGN } from "${LOCAL_SUBWORKFLOWS}/align.nf"
include { CLUSTER } from "${LOCAL_SUBWORKFLOWS}/cluster.nf"
include { PREPARE } from "${LOCAL_SUBWORKFLOWS}/prepare.nf"
include { SCORE } from "${LOCAL_SUBWORKFLOWS}/score.nf"

def blosum = "${projectDir}/files/blosum62.txt"
// https://nextflow-io.github.io/patterns/optional-input/
def nofile = "${projectDir}/files/NO_FILE"

workflow MAIN {

    main:

    PREPARE(
        params.evodists,
        params.cluster,
        params.genomes,
        params.annotations,
        params.long_dist,
        params.medium_dist,
        params.short_dist,
        params.extraexons,
        params.alignmentnum,
        nofile
    )

    ALIGN(
        blosum,
        PREPARE.out.alignment_input,
        PREPARE.out.clusters_split_ch,
        params.long_dist,
        params.medium_dist,
        params.short_dist,
        params.alignmentnum,
        params.prevaln,
        nofile
    )

    SCORE(
        ALIGN.out.folder_jscores,
        PREPARE.out.clusters_split_ch,
        PREPARE.out.dist_ranges_ch,
        params.bonafide_pairs,
        params.long_dist,
        params.medium_dist,
        params.short_dist,
        nofile
    )

    CLUSTER(
        SCORE.out.score_exon_hits_pairs,
        PREPARE.out.clusters_split_ch,
        params.cluster,
        params.orthopairs,
        params.orthogroupnum,
        nofile
    )

    emit:
    exon_cluster_for_reclustering = CLUSTER.out.exon_cluster_for_reclustering

}


