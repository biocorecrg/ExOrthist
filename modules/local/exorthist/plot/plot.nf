process PLOT_EXINT {
    tag "${my_geneID}"
    label 'rscript'
    // containerOptions '-B $PWD:/tmp' # TODO: To consider if this should be moved to nextflow.config

    input:
    val(my_geneID)
    val(my_query_species)
    val(ordered_target)
    val(relevant_exons)
    path(gene_clusters)
    val(isoform_interesting_exs)
    path("*")

    output:
    path("*_exint_plot.pdf"), emit: pdf

    script:
    """
    Rscript ${projectDir}/bin/exint_plotter.R ${my_geneID} ${my_query_species} ${projectDir}/bin ${gene_clusters} ${ordered_target} ${isoform_interesting_exs} ${relevant_exons}
    """
}
