#!/usr/bin/env nextflow

/*
 * Exint plotter for Exhortist pipeline
 *
 */


/*
 * Input parameters:
 */

log.info  """
Executing with the following parameters:

output main				: ${params.output_main}
geneID					: ${params.geneID}
isoformID				: ${params.isoformID}
relevant exons				: ${params.relevant_exs}
reclustered gene orthology file:	: ${params.sub_orthologs}
"""

/*
 * Define params starting from the main.nf output folder
 */


params.geneclusters = "${params.output_main}/gene_cluster_file.gz"
params.annotations = "${params.output_main}/*/*_annot_fake.gtf.gz"
params.overlap = "${params.output_main}/*/*_overlap_CDS_exons.txt"
params.refprot = "${params.output_main}/*/*_ref_proteins.txt"
params.exonclusters = "${params.output_main}/EX_clusters.tab"
params.bestscores = "${params.output_main}/*/best_scored_EX_matches_by_targetgene.txt" //This are all unfiltered scores. I need to identify exons matched by sequence conservation but not phased conservation.

LOCAL_MODULES='../modules/local/exorthist'
include { BREAK_BESTHITS_SPECIESPAIR } from "${LOCAL_MODULES}/exint_break_species.nf"
include { GENERATE_FAKE_COORDS_TABLE } from "${LOCAL_MODULES}/exint_generate_fake.nf"
include { GET_ISOFORM_EXONS } from "${LOCAL_MODULES}/exint_get_exons.nf"
include { ISOLATE_CLUSTER_ID } from "${LOCAL_MODULES}/exint_isolate_clusters.nf"
include { ISOLATE_QUERY_SPECIES } from "${LOCAL_MODULES}/exint_isolate_species.nf"
include { PLOT_EXINT } from "${LOCAL_MODULES}/exint_plot.nf"
include { SELECT_SPECIES_WITH_ORTHOLOGS } from "${LOCAL_MODULES}/exint_select_species.nf"
include { SUBSET_BEST_HITS } from "${LOCAL_MODULES}/exint_subset_best.nf"
include { SUBSET_INPUT_FILES } from "${LOCAL_MODULES}/exint_subset_input.nf"

workflow {
    /*
     * Create channels for input data
     */

    //This channel will contain a list of the GTF files, in theory each with a key
    //The key corresponds to the value assumed by the wildcard in the annotation variable (which is defined in the params.config)
    //annotations  = "$baseDir/data/GTF/*_annot.gtf"
    annotations = Channel.fromFilePairs(params.annotations, size: 1)
        .ifEmpty{error "Cannot find any annotation matching: ${params.annotations}"}

    //The key is the species, same as for the annotations channel
    overlap_info = Channel.fromFilePairs( params.overlap, size: 1)
        .ifEmpty{error "Cannot find any overlap info: ${params.overlap}"}

    //Create channel for files with ref proteins info
    refprot_info = Channel.fromFilePairs(params.refprot, size: 1)
        .ifEmpty{error "Cannot find any overlap info: ${params.refprot}"}

    //Create a joint channel where each key is paired with the corresponding files
    //annotations.join(overlap_info).join(refprot_info).into{all_input_info_raw; all_input_info_raw1}
    all_input_info_raw = annotations.join(overlap_info).join(refprot_info)map{it.flatten()}

    best_hits_input = Channel.fromPath(params.bestscores).toList()
        .ifEmpty{error "Cannot find any overlap info: ${params.bestscores}"}

    /*
     * Isolate clusterID, query_species and species with orthologs
     */
    my_geneID = "${params.geneID}"

    if (params.sub_orthologs) {gene_clusters = file(params.sub_orthologs)} else {gene_clusters = file(params.geneclusters)}

    ISOLATE_CLUSTER_ID(my_geneID, gene_clusters)
    ISOLATE_QUERY_SPECIES(my_geneID, gene_clusters)

    /*
     * Filter only for the species actually having exons in the exon clusters.
     */
    //The pipeline crushes otherwise. But there is a distinction: if a gene gets plotted with no exons, it means that its exons make it to the exon clusters. They simply do not have any ortholog exons to the query species.
    //Interrupt the pipeline in case the query gene is not in the exon clusters
    exon_clusters = file(params.exonclusters)
    SELECT_SPECIES_WITH_ORTHOLOGS(
        exon_clusters,
        ISOLATE_CLUSTER_ID.out.gene_clusterID,
        my_geneID
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
    species_pairs_tmp = Channel
        .fromFilePairs( params.annotations, size: 1).map{it[0]}
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
    GENERATE_FAKE_COORDS_TABLE.out.fake_coords_tables.collect().view()
    BREAK_BESTHITS_SPECIESPAIR.out.best_hits_speciespairs.view()

    // plot_input = GENERATE_FAKE_COORDS_TABLE.out.fake_coords_tables.collect().mix(BREAK_BESTHITS_SPECIESPAIR.out.best_hits_speciespairs).collect()
    // We remove species of the channel
    plot_input = GENERATE_FAKE_COORDS_TABLE.out.fake_coords_tables.collect()
    .mix(BREAK_BESTHITS_SPECIESPAIR.out.best_hits_speciespairs.collect())
    .flatten()
    .filter { it instanceof Path }
    .collect()

    //Get the order of the species for the plot. It can be either provided by the user or computed from the gene cluster file.
    if (params.ordered_species) {
      ordered_target = Channel.value("${params.ordered_species}")
    }
    else {
        all_species = Channel
          .fromFilePairs(params.annotations, size: 1)
          .map{"${it[0]}".toString()}.flatMap().collect().map{it -> it.join(",")}

        DERIVE_ORDERED_SPECIES(
            gene_clusters,
            all_species,
            ISOLATE_QUERY_SPECIES.out.query_species,
            ISOLATE_CLUSTER_ID.out.gene_clusterID
        )
        ordered_target = DERIVE_ORDERED_SPECIES.out.ordered_target
    }

    if (params.isoformID) {
      isoformID = params.isoformID
      GET_ISOFORM_EXONS(
        isoformID,
        GENERATE_FAKE_COORDS_TABLE.out.ExNum_number_in_isoform
            .join(SUBSET_INPUT_FILES.out.overlap_info_4_isoforms)
            .join(ISOLATE_QUERY_SPECIES.out.query_species),
        )
    } else {
      isoform_interesting_exs = Channel.from("None")
    }
    //Rscript to actually make the plot
    if (params.relevant_exs) {relevant_exons = "${params.relevant_exs}"} else {relevant_exons = "None"}

    PLOT_EXINT(
        my_geneID,
        ISOLATE_QUERY_SPECIES.out.query_species,
        ordered_target,
        relevant_exons,
        gene_clusters,
        GET_ISOFORM_EXONS.out.isoform_interesting_exs,
        plot_input
    )
}
