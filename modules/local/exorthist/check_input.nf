process CHECK_INPUT {
    tag "Checking input files"
    label 'publish'

    input:
    path evodisfile
    path clusterfile
    path "GTF/*"
    path "FASTAS/*"
    val gtfs_suffix
    val fastas_suffix
    val long_dist
    val medium_dist
    val short_dist

    output:
    path "run_info.log", emit: run_info

    script:
    """
    echo "Evolutionary distance parameters:" > run_info.log
    echo "long distance: ${long_dist}" >> run_info.log
    echo "medium distance: ${medium_dist}" >> run_info.log
    echo -e "short distance: ${short_dist}\n" >> run_info.log
    echo "Pairwise evolutionary distances:" >> run_info.log
    cat ${evodisfile} >> run_info.log
    echo -e "\nInput files parsing:" >> run_info.log
    A0_check_input.pl -e ${evodisfile} -g GTF -gs ${gtfs_suffix} -f FASTAS -fs ${fastas_suffix} -c ${clusterfile} >> run_info.log
    """
}
