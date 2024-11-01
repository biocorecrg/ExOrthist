process SPLIT_CLUSTERS_BY_SPECIES_PAIRS {
    tag { clusterfile.name }

    input:
    path clusterfile

    output:
    path "*.cls.tab", emit: cls_tab_files
    path "gene_cluster_file.gz", emit: gene_cluster_file

    script:
    """
    if [[ "${clusterfile}" == *.gz ]]; then
        zcat ${clusterfile} > gene_cluster_file
    else
        cat ${clusterfile} > gene_cluster_file
    fi

    A2_split_clusters_by_species_pairs.pl -f gene_cluster_file
    gzip gene_cluster_file
    """
}
