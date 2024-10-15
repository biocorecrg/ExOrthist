process FORMAT_EX_CLUSTERS_OUTPUT {
    publishDir "${params.output}/", mode: 'copy'

    input:
    path "*", stageAs: 'ex_clusters_*'
    path "*", stageAs: 'unclustered_exs_*'

    output:
    path "EX_clusters.tab", emit: exon_cluster_for_reclustering
    path "EX_clusters_info.tab.gz"
    path "unclustered_EXs.txt"

    script:
    """
    D3_format_EX_clusters_output.pl
    """
}
