process PARSE_IPA_PROT_ALN {
    tag { cls_part_file.name }
    label 'aligners'

    input:
    path blosumfile
    tuple val(combid), path(sp1), path(sp2), path(cls_part_file), val(dist_range)

    output:
    tuple val("${sp1.name}-${sp2.name}"), path("${sp1.name}-${sp2.name}-*"), emit: aligned_subclusters_4_splitting
    path "${sp1.name}-${sp2.name}_EXs_to_split_part_*.txt", emit: EXs_to_split

    script:
    def prev_alignments = params.prevaln ? params.prevaln : ""
    def cls_parts = cls_part_file.name.split("_")
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
    B1_parse_IPA_prot_aln.pl ${sp1.name} ${sp2.name} ${cls_part_file} \
    ${sp1}/${sp1.name}_annot_exons_prot_ids.txt ${sp2}/${sp2.name}_annot_exons_prot_ids.txt \
    ${sp1}/${sp1.name}_protein_ids_exons_pos.txt ${sp2}/${sp2.name}_protein_ids_exons_pos.txt \
    ${sp1}/${sp1.name}_protein_ids_intron_pos_CDS.txt ${sp2}/${sp2.name}_protein_ids_intron_pos_CDS.txt \
    ${sp1}/${sp1.name}.exint ${sp2}/${sp2.name}.exint ${cls_parts[1]} ${blosumfile} ${sp1.name}-${sp2.name}-${cls_parts[1]} ${dist_range_par[3]} ${task.cpus} ${prev_alignments}
    """
}
