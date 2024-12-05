process MERGE_PROT_EX_INT_ALN_INFO {
    tag { comp_id }
    label 'aligners'

    stageInMode 'copy'

    // TODO: This needs to be changed once new publish approach
    publishDir "${params.output}", mode: params.publish_dir_mode, pattern: "${comp_id}"

    input:
    tuple val(comp_id), path("FOLDERS_*"), path("*")

    output:
    // TODO: We should only publish certain files, not all directory
    tuple val(comp_id), path("${comp_id}/"), emit: folder_jscores
    path "${comp_id}/all_*_aln_features.txt", emit: aln_features
    path "${comp_id}/EXINT_aln.gz", emit: exint_aln

    script:
    """
    mkdir ${comp_id}
    mv FOLDERS_*/* ${comp_id}/
    mv realigned_* ${comp_id}/

    B4_merge_PROT_EX_INT_aln_info.pl ${comp_id}
    # We don't publish part files
    mkdir -p FOLDERS_
    mv ${comp_id}/*part_* FOLDERS_/
    """
}
