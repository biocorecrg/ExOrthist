process BREAK_BESTHITS_SPECIESPAIR {
    tag "${species_pair}"
    label 'pandas'
    input:
    path best_hits
    tuple val(species_pair), path(query_chosenID), path(target_chosenID)

    output:
    tuple val(species_pair), path("${species_pair}-best_scores_with_overlapIDs.txt"), emit: best_hits_speciespairs

    script:
    def (species1, species2) = species_pair.split("_")
    """
    awk 'NR==1' ${best_hits} > ${species_pair}-best_scores_hits_exons.txt
    cat ${best_hits} | awk '\$16=="${species1}" && \$17=="${species2}"' >> ${species_pair}-best_scores_hits_exons.txt
    add_selected_overlapID_to_besthits.py \
        -b ${species_pair}-best_scores_hits_exons.txt \
        -q ${query_chosenID} \
        -t ${target_chosenID} \
        -out ${species_pair}-best_scores_with_overlapIDs.txt
    """
}
