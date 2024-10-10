process GENERATE_ANNOTATIONS {
    tag { genomeid }
    label 'big_cpus'
    publishDir "${params.output}/", mode: 'copy'

    input:
    tuple val(genomeid), path(genome), path(annotation)
    path extraexons

    output:
    tuple val(genomeid), path("${genomeid}"), emit: idfolders

    script:
    // Sic. https://nextflow-io.github.io/patterns/optional-input/
    def extrapars = extraexons.name != 'NO_FILE' ? "-add_exons ${extraexons}" : ""
    """
    A1_generate_annotations.pl -GTF ${annotation} -G ${genome} -sp ${genomeid} ${extrapars}
    """
}
