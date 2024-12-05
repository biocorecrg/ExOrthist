process SPLIT_CLUSTERS_IN_CHUNKS {
    tag { "${idfolder_A}_${idfolder_B}" }

    input:
    path cls_tab_files
    tuple val(id_comb), path(idfolder_A), path(idfolder_B)
    val(alignmentnum)

    output:
    tuple path(idfolder_A), path(idfolder_B), path("${idfolder_A}_${idfolder_B}/*.cls.tab-part_*"), emit: cls_files_2_align

    script:
    """
    A3_split_clusters_in_chunks.pl \
        --sp1 ${idfolder_A} \
        --sp2 ${idfolder_B} \
        --expath ./ \
        --project_dir ./ \
        --N_split ${alignmentnum} \
        --gene_cluster ${id_comb}.cls.tab
    """
}
