process SCORE_EX_MATCHES {
    tag { comp_id }
    label 'big_mem'

    storeDir "${params.output}/${comp_id}"

    input:
    tuple val(comp_id), path("*")

    output:
    path "all_PROT_EX_INT_aln_features_*", emit: all_features
    tuple val(comp_id), path("all_scored_EX_matches.txt"), emit: all_scores_to_filt

    script:
    def (species1, species2) = comp_id.split("-")
    """
    B5_format_aln_info_by_best_isoform_match.pl ${species1} ${species2} \
    ${comp_id}/all_PROT_aln_features.txt ${comp_id}/all_EX_aln_features.txt ${comp_id}/all_INT_aln_features.txt \
    ${species1}/${species1}.exint ${species2}/${species2}.exint \
    ${species1}/${species1}_protein_ids_intron_pos_CDS.txt ${species2}/${species2}_protein_ids_intron_pos_CDS.txt \
    all_PROT_EX_INT_aln_features_${comp_id}.txt

    C1_score_EX_matches.pl all_PROT_EX_INT_aln_features_${comp_id}.txt .
    """
}
