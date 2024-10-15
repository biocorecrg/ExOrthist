process RECLUSTER_GENES_BY_SPECIES_PAIR {
    publishDir "${params.output}/reclustering", mode: 'copy'
    label 'pandas'
    tag { combid }

    input:
    tuple val(combid), path(folderA), path(folderB)
    path clusterfile
    path orthopairs

    output:
    tuple val(combid), path("reclustered_genes_*.tab"), emit: recl_genes_for_rec_exons

    when:
    orthopairs.name != 'NO_FILE'

    script:
    def (species1, species2) = combid.split("-")
    """
    D3.1_recluster_genes_by_species_pair.py \
        -og ${clusterfile} \
        -op ${orthopairs} \
        --species1 ${species1} \
        --species2 ${species2} \
        -out reclustered_genes_${combid}.tab
    """
}
