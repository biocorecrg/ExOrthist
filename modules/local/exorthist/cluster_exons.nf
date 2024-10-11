process CLUSTER_EXS {
    label 'rscript'
    input:
    path cluster_part

    output:
    path "EXs_${cluster_part.name}", emit: ex_clusters
    path "unclustered_EXs_${cluster_part.name}", emit: unclustered_exs

    script:
    """
    D2_cluster_EXs.R ${cluster_part} EXs_${cluster_part.name} unclustered_EXs_${cluster_part.name}
    """
}
