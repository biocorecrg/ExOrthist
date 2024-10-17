process ISOLATE_QUERY_SPECIES {
    tag { my_geneID }
    label 'pandas'

    input:
    val(my_geneID)
    path gene_clusters

    output:
    stdout emit: query_species

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    my_df = pd.read_csv("${gene_clusters}", sep="\t", header=None, index_col=False)
    query_species = list(my_df[my_df.iloc[:,2]=="${my_geneID}"].iloc[:,1])[0]
    print(str(query_species), end='')
    """
}
