process PLOT_EXINT {
    tag "${my_geneID}"
    label 'rscript'
    // containerOptions '-B $PWD:/tmp' # TODO: To consider if this should be moved to nextflow.config
    publishDir "${params.output_plot}", mode: 'copy'
    //publishDir "${params.output}/${params.geneID}", mode: 'copy'

    input:
    val(my_geneID)
    val(my_query_species)
    val(ordered_target)
    val(relevant_exons)
    path(gene_clusters)
    val(isoform_interesting_exs)
    path("*")

    output:
    path("*_exint_plot.pdf")

    script:
    """
    Rscript $baseDir/bin/exint_plotter.R ${my_geneID} ${my_query_species} ${baseDir}/bin ${gene_clusters} ${ordered_target} ${isoform_interesting_exs} ${relevant_exons}
    """
}
