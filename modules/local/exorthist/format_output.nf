process FORMAT_EX_CLUSTERS_OUTPUT {
    label 'publish'

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
