process SUBSET_BEST_HITS {
    tag "${gene_clusterID}"
    label 'big_mem'

    input:
    val(gene_clusterID)
    path("best_hits_species_pairs_*")

    output:
    path("Best_hits_subsetted"), emit: best_hits

    script:
    """
    cat best_hits_species_pairs_* > best_hits.tmp
    awk 'NR==1' best_hits.tmp > Best_hits_subsetted
    cat best_hits.tmp | awk '\$1=="${gene_clusterID}"' >> Best_hits_subsetted
    rm best_hits.tmp
    """
}
