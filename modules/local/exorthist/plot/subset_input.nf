process SUBSET_INPUT_FILES {
    tag "${gene_clusterID}"
    label 'pandas'

    input:
    val(gene_clusterID)
    path(gene_clusters)
    tuple val(species), path(annotations), path(overlap_info), path(ref_proteins)

    output:
    tuple val(species), path("*_subsetted_annot.gtf"), path("*_subsetted_overlap_info.txt"), path("*_subsetted_ref_proteins.txt"), emit: all_subsetted_inputs
    tuple val(species), path("*_subsetted_overlap_info.txt"), emit: overlap_info_4_isoforms

    script:
    """
    subset_inputs.py -a ${annotations} -o ${overlap_info} -r ${ref_proteins} -c ${gene_clusters} -g ${gene_clusterID} -s ${species}
    """
}
