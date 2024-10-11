process SPLIT_EX_PAIRS_TO_REALIGN {
    label 'pandas'
    input:
    path '*'

    output:
    path '*EXs_to_realign_part_*', emit: EXs_to_realign_batches

    script:
    """
    for file in \$(ls *); do
        B2_split_EX_pairs_to_realign.py -i \${file} -n ${params.alignmentnum}
    done
    """
}
