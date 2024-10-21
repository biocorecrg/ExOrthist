LOCAL_MODULES='../../../modules/local/exorthist'

include { COLLAPSE_OVERLAPPING_MATCHES } from "${LOCAL_MODULES}/collapse_matches.nf"
include { FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE } from "${LOCAL_MODULES}/filter_matches.nf"
include { JOIN_FILTERED_EX_MATCHES } from "${LOCAL_MODULES}/join_matches.nf"
include { SCORE_EX_MATCHES } from "${LOCAL_MODULES}/score_matches.nf"

workflow SCORE {

    take:
    folder_jscores
    clusters_split_ch
    dist_ranges_ch
    bonafide_pairs

    main:

    data_to_score = folder_jscores.join(clusters_split_ch).map{ [it[0], it[1..-1] ]}
    // Score EX matches from aln info
    SCORE_EX_MATCHES(data_to_score)
    // Filter the best matches above score cutoffs by target gene.
    all_scores_to_filt_ch = SCORE_EX_MATCHES.out.all_scores_to_filt
    FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE(all_scores_to_filt_ch.join(dist_ranges_ch))
    //  Join filtered scored EX matches
    filterscore_per_joining_ch = FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE.out.filterscore_per_joining
    JOIN_FILTERED_EX_MATCHES(filterscore_per_joining_ch.collect())
    //  Removing matches from overlapping exons

    filtered_all_scores = JOIN_FILTERED_EX_MATCHES.out.filtered_all_scores
    // Sic: https://nextflow-io.github.io/patterns/optional-input/
    if ( bonafide_pairs ) {
        COLLAPSE_OVERLAPPING_MATCHES(filtered_all_scores, bonafide_pairs)
    } else {
        COLLAPSE_OVERLAPPING_MATCHES(filtered_all_scores, Channel.fromPath("/path/to/NO_FILE").collect())
    }

    emit:
    score_exon_hits_pairs = COLLAPSE_OVERLAPPING_MATCHES.out.score_exon_hits_pairs
}
