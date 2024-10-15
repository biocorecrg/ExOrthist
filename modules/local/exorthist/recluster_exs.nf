process RECLUSTER_EXS_BY_SPECIES_PAIR {
    publishDir "${params.output}/reclustering", mode: 'copy'
    label 'pandas'
    tag { combid }

    input:
    tuple val(combid), path(recl_genes)
    path exon_clusters
    path exon_pairs
    path orthopairs

    output:
    path "reclustered_EXs_*.tab", emit: recl_exs

    when:
    orthopairs.name != 'NO_FILE'

    script:
    def combid1 = combid.replace("-", "_")
    def (species1, species2) = combid1.split("_")
    """
    D3.2_recluster_EXs_by_species_pair.py \
        -ep ${exon_pairs} \
        -rg ${recl_genes} \
        -ec ${exon_clusters} \
        -sp1 ${species1} \
        -sp2 ${species2} \
        -out reclustered_EXs_${combid1}.tab
    """
}
