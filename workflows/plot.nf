LOCAL_MODULES='../modules/local/exorthist/plot'
include { BREAK_BESTHITS_SPECIESPAIR } from "${LOCAL_MODULES}/break_species.nf"
include { GENERATE_FAKE_COORDS_TABLE } from "${LOCAL_MODULES}/generate_fake.nf"
include { GET_ISOFORM_EXONS } from "${LOCAL_MODULES}/get_exons.nf"
include { ISOLATE_CLUSTER_ID } from "${LOCAL_MODULES}/isolate_clusters.nf"
include { ISOLATE_QUERY_SPECIES } from "${LOCAL_MODULES}/isolate_species.nf"
include { PLOT_EXINT } from "${LOCAL_MODULES}/plot.nf"
include { SELECT_SPECIES_WITH_ORTHOLOGS } from "${LOCAL_MODULES}/select_species.nf"
include { SUBSET_BEST_HITS } from "${LOCAL_MODULES}/subset_best.nf"
include { SUBSET_INPUT_FILES } from "${LOCAL_MODULES}/subset_input.nf"

workflow PLOT {
    take:
    geneID
    gene_clusters
    annotations
    all_input_info_raw
    best_hits_input
    exon_clusters
    relevant_exons
    ordered_species
    isoformID

    main:
    /*
     * Create channels for input data
     */

    ISOLATE_CLUSTER_ID(geneID, gene_clusters)
    ISOLATE_QUERY_SPECIES(geneID, gene_clusters)

    /*
     * Filter only for the species actually having exons in the exon clusters.
     */
    //The pipeline crushes otherwise. But there is a distinction: if a gene gets plotted with no exons, it means that its exons make it to the exon clusters. They simply do not have any ortholog exons to the query species.
    //Interrupt the pipeline in case the query gene is not in the exon clusters
    SELECT_SPECIES_WITH_ORTHOLOGS(
        exon_clusters,
        ISOLATE_CLUSTER_ID.out.gene_clusterID,
        geneID
    )
    single_species = SELECT_SPECIES_WITH_ORTHOLOGS.out.selected_species.map{it -> it.split(",")}.flatMap()
    all_input_info = all_input_info_raw.join(single_species)

    // /*
    //  * Subset input files for genes belonging to the gene cluster of interest
    //  */

    SUBSET_INPUT_FILES(
        ISOLATE_CLUSTER_ID.out.gene_clusterID,
        gene_clusters,
        all_input_info
    )

    SUBSET_BEST_HITS(
        ISOLATE_CLUSTER_ID.out.gene_clusterID,
        best_hits_input
    )

    /*
     * Generate exint plotter input for each species
     */

    GENERATE_FAKE_COORDS_TABLE(
        SUBSET_INPUT_FILES.out.all_subsetted_inputs,
        gene_clusters,
        exon_clusters
    )

    // /*
    //  * Generate best hits input files
    //  */
    //Generate a channel with all combinations of species query with each species target
    species_pairs_tmp = annotations.map{it[0]}
        .toList()
        .map{[it, it].combinations().findAll{a,b -> a!=b}}
        .flatMap()
        .groupTuple()
        .join(ISOLATE_QUERY_SPECIES.out.query_species)
        .transpose()

    //Add the needed files from the previous process
    overlapID_chosenID = GENERATE_FAKE_COORDS_TABLE.out.overlapID_chosenID

    all_species_pairs = overlapID_chosenID.cross(species_pairs_tmp).map{it.flatten()}
        .map{[it[3], it[0], it[1], it[2]]}
        .join(overlapID_chosenID)
        .map{["${it[3]}_${it[0]}".toString(), it[2], it[4]]}

    BREAK_BESTHITS_SPECIESPAIR(
        SUBSET_BEST_HITS.out.best_hits,
        all_species_pairs
    )

    /*
     * Facultative processes for isoform and exons highlighting
     */
    // We remove species of the channel
    plot_input = GENERATE_FAKE_COORDS_TABLE.out.fake_coords_tables.collect()
    .mix(BREAK_BESTHITS_SPECIESPAIR.out.best_hits_speciespairs.collect())
    .flatten()
    .filter { it instanceof Path }
    .collect()

    //Get the order of the species for the plot. It can be either provided by the user or computed from the gene cluster file.
    if (ordered_species) {
      ordered_target = Channel.value("${ordered_species}")
    }
    else {
        all_species = Channel
          .fromFilePairs(annotations, size: 1)
          .map{"${it[0]}".toString()}.flatMap().collect().map{it -> it.join(",")}

        DERIVE_ORDERED_SPECIES(
            gene_clusters,
            all_species,
            ISOLATE_QUERY_SPECIES.out.query_species,
            ISOLATE_CLUSTER_ID.out.gene_clusterID
        )
        ordered_target = DERIVE_ORDERED_SPECIES.out.ordered_target
    }

    if (isoformID) {
      GET_ISOFORM_EXONS(
        isoformID,
        GENERATE_FAKE_COORDS_TABLE.out.ExNum_number_in_isoform
            .join(SUBSET_INPUT_FILES.out.overlap_info_4_isoforms)
            .join(ISOLATE_QUERY_SPECIES.out.query_species),
        geneID
      )
    } else {
      isoform_interesting_exs = Channel.from("None")
    }
    //Rscript to actually make the plot

    PLOT_EXINT(
        geneID,
        ISOLATE_QUERY_SPECIES.out.query_species,
        ordered_target,
        relevant_exons,
        gene_clusters,
        GET_ISOFORM_EXONS.out.isoform_interesting_exs,
        plot_input
    )

}
