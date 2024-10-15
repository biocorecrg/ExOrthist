process CLUSTER_EXS {
    label 'rscript'
    input:
    path cluster_parts

    output:
    path "EXs_${cluster_parts.name}", emit: ex_clusters
    path "unclustered_EXs_${cluster_parts.name}", emit: unclustered_exs

    script:
    """
    D2_cluster_EXs.R ${cluster_parts} EXs_${cluster_parts.name} unclustered_EXs_${cluster_parts.name}
    """
}
