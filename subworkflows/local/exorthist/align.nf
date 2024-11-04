LOCAL_MODULES='../../../modules/local/exorthist'

include { MERGE_PROT_EX_INT_ALN_INFO } from "${LOCAL_MODULES}/merge_aligns.nf"
include { PARSE_IPA_PROT_ALN } from "${LOCAL_MODULES}/align_pairs.nf"
include { REALIGN_EX_PAIRS } from "${LOCAL_MODULES}/realign_pairs.nf"
include { SPLIT_EX_PAIRS_TO_REALIGN } from "${LOCAL_MODULES}/split_pairs.nf"

workflow ALIGN {

    take:
    blosumfile
    alignment_input
    clusters_split_ch
    long_dist
    medium_dist
    short_dist
    alignmentnum
    prevaln
    outdir

    main:

    blosumfile_ch = Channel.fromPath(blosumfile, checkIfExists: true).collect()

    if (prevaln) {
        prevaln_ch = Channel.fromPath(prevaln, type: 'dir', checkIfExists: true).collect()
    } else {
        prevaln_ch = Channel.fromPath("/path/to/NO_FILE").collect()
    }

    // the last argument is the protein similarity alignment.
    // if a prevaln folder is provided, the protein alignments present in each species pair subfolder will not be repeated.
    PARSE_IPA_PROT_ALN(blosumfile_ch, alignment_input, long_dist, medium_dist, short_dist, prevaln_ch)

    // Collapse EXs_to_split in batches of 500 files
    EXs_to_split = PARSE_IPA_PROT_ALN.out.EXs_to_split
    EXs_to_split_batches = EXs_to_split.toSortedList().flatten().buffer(size : 500, remainder: true)
    // Split exons pairs to realign
    SPLIT_EX_PAIRS_TO_REALIGN(EXs_to_split_batches, alignmentnum)
    EXs_to_realign_batches = SPLIT_EX_PAIRS_TO_REALIGN.out.EXs_to_realign_batches
    // Flatten the results from the previous batch run and combine with sp1 and sp2 information, using sp1-sp2 as key.
    EXs_to_realign = EXs_to_realign_batches.flatten().map{[it.getName().toString().split("_")[0],it]}.groupTuple().join(clusters_split_ch).transpose()
    //  Realign exons pairs (with multiple hits)
    REALIGN_EX_PAIRS(blosumfile_ch, EXs_to_realign)
    // Combine all the aln_info with the realigned_exon_info for each species pair
    aligned_subclusters_4_splitting = PARSE_IPA_PROT_ALN.out.aligned_subclusters_4_splitting
    realigned_exons_4_merge = REALIGN_EX_PAIRS.out.realigned_exons_4_merge
    data_4_merge = aligned_subclusters_4_splitting.groupTuple().join(realigned_exons_4_merge.groupTuple())
    // Merge alignments information
    MERGE_PROT_EX_INT_ALN_INFO(data_4_merge, outdir)

    emit:
    folder_jscores = MERGE_PROT_EX_INT_ALN_INFO.out.folder_jscores
    aln_features = MERGE_PROT_EX_INT_ALN_INFO.out.aln_features
    exint_aln = MERGE_PROT_EX_INT_ALN_INFO.out.exint_aln

}
