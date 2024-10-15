process FORMAT_EX_CLUSTERS_OUTPUT {
    publishDir "${params.output}/", mode: 'copy'

    input:
    path ex_clusters
    path unclustered_exs

    output:
    path "EX_clusters.tab", emit: exon_cluster_for_reclustering
    path "EX_clusters_info.tab.gz"
    path "unclustered_EXs.txt"

    script:
    """
    D3_format_EX_clusters_output.pl
    """
}
