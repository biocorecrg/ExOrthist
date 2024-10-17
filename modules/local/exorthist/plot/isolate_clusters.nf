process ISOLATE_CLUSTER_ID {
    tag { my_geneID }
    label 'pandas'

    input:
    val(my_geneID)
    path gene_clusters

    output:
    stdout emit: gene_clusterID

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    my_df = pd.read_csv("${gene_clusters}", sep="\t", header=None, index_col=False)
    clusterID = list(my_df[my_df.iloc[:,2]=="${my_geneID}"].iloc[:,0])[0]
    print(clusterID, end='')
    """
}
