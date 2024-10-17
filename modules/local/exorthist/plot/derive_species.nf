process DERIVE_ORDERED_SPECIES {
    label 'pandas'

    input:
    path gene_clusters
    val all_species
    val query_species
    val gene_clusterID

    output:
    stdout emit: ordered_target

    script:
    """
    #!/usr/bin/env python
    import pandas as pd

    my_df = pd.read_csv("${gene_clusters}", sep="\t", header=None, index_col=False)
    species_list = list(my_df[my_df.iloc[:,0]=="${gene_clusterID}"].iloc[:,1])
    interesting_species = list(set([element for element in species_list if element in str("${all_species}").split(",")]))
    interesting_species.insert(0, interesting_species.pop(interesting_species.index("${query_species}")))
    final_species_list = ",".join(str(element) for element in interesting_species)
    print(final_species_list, end='')
    """
}
