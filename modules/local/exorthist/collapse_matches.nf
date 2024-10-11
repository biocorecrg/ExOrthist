process COLLAPSE_OVERLAPPING_MATCHES {
    publishDir "${params.output}/", mode: "copy"

    input:
    path scores
    path bonafide_pairs

    output:
    path "filtered_best_scored_EX_matches_by_targetgene-NoOverlap.tab", emit: score_exon_hits_pairs
    path "overlapping_EXs_by_species.tab", emit: overlapping_exs

    script:
    def bonafide = bonafide_pairs.name != 'NO_FILE' ? "-b ${bonafide_pairs}" : ''
    """
    C3_count_matches_by_EX.pl ${scores} EX_matches_count_by_species.tab ${bonafide}
    C4_get_overlapping_EXs.pl -i EX_matches_count_by_species.tab -o overlapping_EXs_by_species.tab
    C5_collapse_overlapping_matches.pl overlapping_EXs_by_species.tab ${scores} filtered_best_scored_EX_matches_by_targetgene-NoOverlap.tab ${bonafide}
    """
}
