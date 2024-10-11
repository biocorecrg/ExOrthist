process FILTER_AND_SELECT_BEST_EX_MATCHES_BY_TARGETGENE {
    tag { comp_id }

    publishDir "${params.output}", mode: "copy",
        pattern: "best_scored_EX_matches_by_targetgene.txt",
        saveAs: { filename -> "${comp_id}/foo_$filename" }

    input:
    tuple val(comp_id), path(all_scores), val(dist_range)

    output:
    path "*.tab", emit: filterscore_per_joining
    path "best_scored_EX_matches_by_targetgene.txt", emit: best_scored_matches

    script:
    def (species1, species2) = comp_id.split("-")
    def dist_range_par

    switch(dist_range) {
        case "long":
            dist_range_par = params.long_dist.split(",")
            break
        case "medium":
            dist_range_par = params.medium_dist.split(",")
            break
        case "short":
            dist_range_par = params.short_dist.split(",")
            break
    }

    """
    C2_filter_and_select_best_EX_matches_by_targetgene.pl \
        -b ${all_scores} \
        -sps ${species1},${species2} \
        -int ${dist_range_par[0]} \
        -id ${dist_range_par[1]} \
        -max_size ${dist_range_par[2]}
    """
}
