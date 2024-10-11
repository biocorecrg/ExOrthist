process FORMAT_EX_CLUSTERS_INPUT {
    input:
    path score_exon_hits_pairs
    path clusterfile

    output:
    path "PART_*-cluster_input.tab", emit: cluster_parts

    script:
    """
    if [[ "${clusterfile}" == *.gz ]]; then
        zcat ${clusterfile} > cluster_file
        D1_format_EX_clusters_input.pl cluster_file ${score_exon_hits_pairs} ${params.orthogroupnum}
        rm cluster_file
    else
        D1_format_EX_clusters_input.pl ${clusterfile} ${score_exon_hits_pairs} ${params.orthogroupnum}
    fi
    """
}
