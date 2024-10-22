LOCAL_MODULES='../../../modules/local/exorthist'

include { CHECK_INPUT } from "${LOCAL_MODULES}/check_input.nf"
include { GENERATE_ANNOTATIONS } from "${LOCAL_MODULES}/generate_annotations.nf"
include { SPLIT_CLUSTERS_IN_CHUNKS } from "${LOCAL_MODULES}/split_clusters_chunks.nf"
include { SPLIT_CLUSTERS_BY_SPECIES_PAIRS } from "${LOCAL_MODULES}/split_clusters_species.nf"

workflow PREPARE {

    take:
    evodists
    clusterfile_ch
    gtfs
    fastas
    gtfs_suffix
    fastas_suffix
    long_dist
    medium_dist
    short_dist
    genomes
    annotations
    extraexons
    alignmentnum

    main:

    evodists_ch = Channel.fromPath(evodists, checkIfExists: true).collect()

    CHECK_INPUT(
        evodists_ch,
        clusterfile_ch,
        gtfs,
        fastas,
        gtfs_suffix,
        fastas_suffix,
        long_dist,
        medium_dist,
        short_dist
    )

    extraexons_ch = extraexons ?
        Channel.fromFilePairs(extraexons, checkIfExists: true, size: 1)
        .ifEmpty { error "Extra exons not found" } :
        Channel.empty()

    // We join channels. If no extraexons, then it's empty, so no problem
    data_to_annotation = genomes.join(annotations)
    // TODO: To check what's goiing on here - check extraexons_ch
    // data_to_annotation = data_to_annotation.join(extraexons_ch, remainder: true)

    if (extraexons) {
        GENERATE_ANNOTATIONS(data_to_annotation, extraexons_ch)
    } else {
        // Sic: https://nextflow-io.github.io/patterns/optional-input/
        GENERATE_ANNOTATIONS(data_to_annotation, Channel.fromPath("/path/to/NO_FILE").collect())
    }

    clusters_split_ch = GENERATE_ANNOTATIONS.out.idfolders.toList().map{ [it, it].combinations().findAll{ a, b -> a[0] < b[0]} }
        .flatMap()
        .map { ["${it[0][0]}-${it[1][0]}".toString(), it[0][1], it[1][1]] }

    // Copy the gene cluster file to output to use for the exint_plotter and compare_exon_sets modules
    SPLIT_CLUSTERS_BY_SPECIES_PAIRS(clusterfile_ch)

    // Split clusters
    cls_tab_files_ch = SPLIT_CLUSTERS_BY_SPECIES_PAIRS.out.cls_tab_files
    SPLIT_CLUSTERS_IN_CHUNKS(cls_tab_files_ch.collect(), clusters_split_ch, alignmentnum)

    cls_files_2_align = SPLIT_CLUSTERS_IN_CHUNKS.out.cls_files_2_align
    cls_files_2_align_t = cls_files_2_align.transpose().map{[it[0].getFileName().toString()+"-"+it[1].getFileName().toString(), it[0], it[1], it[2]]}

    //Create a channel for the evo distances
    sp1_sp2_dist = evodists_ch.flatten()
     .splitText()
     .map{"${it}".trim().split("\t")}.map{[it[0]+"-"+it[1], it[2]]}

    sp2_sp1_dist = evodists_ch.flatten()
     .splitText()
     .map{"${it}".trim().split("\t")}.map{[it[1]+"-"+it[0], it[2]]}

    species_pairs_dist = sp1_sp2_dist.concat(sp2_sp1_dist)

    //Only the species pairs with a common index will be kept
    dist_ranges_ch = clusters_split_ch.join(species_pairs_dist).map{[it[0], it[3]]}
    alignment_input = cls_files_2_align_t.groupTuple().join(dist_ranges_ch).transpose()

    emit:
    clusters_split_ch
    dist_ranges_ch
    alignment_input
}
