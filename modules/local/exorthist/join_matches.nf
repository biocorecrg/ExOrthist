process JOIN_FILTERED_EX_MATCHES {
    publishDir "${params.output}/", mode: 'copy'

    input:
    path "filtered_best_scored-*"

    output:
    path "filtered_best_scored_EX_matches_by_targetgene.tab", emit: filtered_all_scores

    script:
    """
    echo "GeneID_sp1\tExon_coords_sp1\tGeneID_sp2\tExon_coords_sp2\tSp1\tSp2" > filtered_best_scored_EX_matches_by_targetgene.tab
    for file in \$(ls filtered_best_scored-*); do
        cat \$file | tail -n+2 >> filtered_best_scored_EX_matches_by_targetgene.tab
    done
    """
}
