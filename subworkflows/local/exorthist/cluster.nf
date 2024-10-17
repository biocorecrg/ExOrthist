LOCAL_MODULES='../../../modules/local/exorthist'

include { CLUSTER_EXS } from "${LOCAL_MODULES}/cluster_exons.nf"
include { FORMAT_EX_CLUSTERS_INPUT } from "${LOCAL_MODULES}/format_input.nf"
include { FORMAT_EX_CLUSTERS_OUTPUT } from "${LOCAL_MODULES}/format_output.nf"
include { RECLUSTER_EXS_BY_SPECIES_PAIR } from "${LOCAL_MODULES}/recluster_exs.nf"
include { RECLUSTER_GENES_BY_SPECIES_PAIR } from "${LOCAL_MODULES}/recluster_genes.nf"

workflow CLUSTER {

    take:
    score_exon_hits_pairs
    cluster
    clusters_split_ch
    clusterfile_ch
    orthopairs_ch

    main:

    FORMAT_EX_CLUSTERS_INPUT(score_exon_hits_pairs, cluster)

    // Split the file of exon pairs
    // Unclustered are the exons ending up in single-exon clusters
    cluster_parts = FORMAT_EX_CLUSTERS_INPUT.out.cluster_parts.flatten()
    CLUSTER_EXS(cluster_parts)
    FORMAT_EX_CLUSTERS_OUTPUT(CLUSTER_EXS.out.ex_clusters.collect(), CLUSTER_EXS.out.unclustered_exs.collect())

    // Re-clustering of genes
    RECLUSTER_GENES_BY_SPECIES_PAIR(
        clusters_split_ch,
        clusterfile_ch,
        orthopairs_ch
    )

    // Re-clustering of exons
    RECLUSTER_EXS_BY_SPECIES_PAIR(
        RECLUSTER_GENES_BY_SPECIES_PAIR.out.recl_genes_for_rec_exons,
        FORMAT_EX_CLUSTERS_OUTPUT.out.exon_cluster_for_reclustering,
        score_exon_hits_pairs,
        orthopairs_ch
    )

    emit:
    recl_exs = RECLUSTER_EXS_BY_SPECIES_PAIR.out.recl_exs
    recl_genes_for_rec_exons = RECLUSTER_GENES_BY_SPECIES_PAIR.out.recl_genes_for_rec_exons

}
