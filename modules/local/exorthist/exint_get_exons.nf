process GET_ISOFORM_EXONS {
    tag "${isoformID}"
    label 'pandas'
    publishDir "${params.output}/${params.geneID}/processed_table", mode: 'copy'

    input:
    val(isoformID)
    tuple val(species), path(exons_info), path(overlap_info)

    output:
    stdout emit: isoform_interesting_exs

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import re

    my_df = pd.read_csv("${exons_info}", sep="\t", header=None, index_col=False, names=["ExonID", "ExonNum", "GeneID", "ProteinID"])
    interesting_exons = list(my_df.loc[my_df["ProteinID"]=="${isoformID}"]["ExonID"])
    #select all the exons from the same overlapping group
    my_chr = [re.sub(":.*", "", element) for element in interesting_exons][0]
    my_overlap_df = pd.read_csv("${overlap_info}", sep="\t", header=None, index_col=False, names=["OverlapID", "GeneID", "Start_Stop"])
    my_overlap_df["Coords"] = [my_chr+":"+element for element in list(my_overlap_df["Start_Stop"])]
    overlapping_groups = list(my_overlap_df.loc[my_overlap_df["Coords"].isin(interesting_exons)]["OverlapID"])
    all_interesting_exons = list(my_overlap_df.loc[my_overlap_df["OverlapID"].isin(overlapping_groups)]["Coords"])
    all_interesting_exons = ["${isoformID}"] + all_interesting_exons
    print(",".join(str(element) for element in all_interesting_exons), end='')
    """
}
