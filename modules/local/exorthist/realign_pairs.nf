process REALIGN_EX_PAIRS {
    label 'incr_time_cpus'

    input:
    path blosumfile
    tuple val(comp_id), path(EXs_to_realign), path(sp1), path(sp2)

    output:
    tuple val(comp_id), path("realigned_*"), emit: realigned_exons_4_merge

    script:
    """
    B3_realign_EX_pairs.pl ${sp1.name} ${sp2.name} ${EXs_to_realign} \
    ${sp1}/${sp1.name}.exint ${sp2}/${sp2.name}.exint 1 realigned_${EXs_to_realign.name} \
    ${sp1.name}_${sp2.name} ${blosumfile} ${task.cpus}
    """
}
