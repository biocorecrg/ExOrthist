process SELECT_SPECIES_WITH_ORTHOLOGS {
    tag { clusterID }
    label 'pandas'

    input:
    path my_exon_clusters
    val clusterID
    val my_geneID

    output:
    stdout emit: selected_species

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import re
    import sys

    my_df = pd.read_table("${my_exon_clusters}", sep="\t", header=0, index_col=False)
    my_df["GeneClusterID"] = [re.sub("\\..*", "", element) for element in my_df["ExCluster_ID"]]
    gene_ids_list = list(my_df[my_df.GeneClusterID=="${clusterID}"]["GeneID"])
    selected_species = list(set(list(my_df[my_df.GeneClusterID=="${clusterID}"]["Species"])))
    if "${my_geneID}" in gene_ids_list and len(selected_species) >=2:
        print(','.join(map(str, selected_species)), end='')
    else:
        sys.exit("Query gene does not have any exon orthologs in the other species")
    """
}
