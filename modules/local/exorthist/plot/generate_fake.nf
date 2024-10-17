process GENERATE_FAKE_COORDS_TABLE {
    tag "${species}"
    label 'pandas'

    input:
    tuple val(species), path(annotations), path(overlap_info), path(ref_prot)
    path gene_clusters
    path exon_clusters

    output:
    tuple val(species), path("${species}_exons_cluster_info-fakecoords.tab"), emit: fake_coords_tables
    tuple val(species), path("${species}_ExNum_in_isoform"), emit: ExNum_number_in_isoform
    tuple val(species), path("${species}_overlapID_chosenID.txt"), emit: overlapID_chosenID

    script:
    """
    generate_exint_plotter_input.py \
        -a ${annotations} \
        -o ${overlap_info} \
        -c ${gene_clusters} \
        -e ${exon_clusters} \
        -r ${ref_prot} \
        -s ${species} \
        -out ${species}_exons_cluster_info-fakecoords.tab
    """
}
