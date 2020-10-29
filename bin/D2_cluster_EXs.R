#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
if (length(args)<2) {stop("[USAGE] Rscript --vanilla cluster.R <file1> <file2> ")}
file1 <- args[1]
file2 <- args[2]

library("igraph", quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library("hashmap", quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)

my_table = read.table(file1, col.names=c("ID1", "ID2", "Best_reciprocal"))
######## Generate graph #############
my_matrix = as.matrix(my_table[,c("ID1", "ID2")]) #select columns with IDs
my_edges = graph_from_edgelist(my_matrix, directed=TRUE) #create directed graph starting from input.
my_clusters = cluster_edge_betweenness(my_edges, directed=FALSE) #Do not consider directionality in deriving the clusters. But the directionality is considered for everything else.
my_final_clusters = as.data.frame(as.matrix(membership(my_clusters))); colnames(my_final_clusters) = "ClusterID" #communities

######## Exon re-intronduction #############
#Consider cases of exons which get in single-exon clusters but actually present a reciprocal best with an exon in another cluster.
#Add this exon to the cluster of (one of) its reciprocal best(s).
single_exon_cluster_IDs = names(table(my_final_clusters)[table(my_final_clusters)==1])
for (my_clusterID in single_exon_cluster_IDs) {
  exon_id = rownames(subset(my_final_clusters, ClusterID==my_clusterID))
  exon_hits = subset(my_table, ID1 == exon_id | ID2 == exon_id & Best_reciprocal == "TRUE")
  if (nrow(exon_hits) >= 1) { #take the first best_hit (if there are more in different exon clusters there's not much we can do)
    if (exon_hits[1,"ID1"] == exon_id) {new_cluster_ID = my_final_clusters[as.vector(exon_hits[1,"ID2"]), "ClusterID"]} #take the ID of the matching exon
    if (exon_hits[1,"ID2"] == exon_id) {new_cluster_ID = my_final_clusters[as.vector(exon_hits[1,"ID1"]), "ClusterID"]} #take the ID of the matching exon
    my_final_clusters[exon_id,"ClusterID"] = new_cluster_ID
    } else {
      my_final_clusters$names = rownames(my_final_clusters)
      my_final_clusters = my_final_clusters[!(rownames(my_final_clusters) == exon_id),] #remove the single-exon cluster from the final graph.
      rownames(my_final_clusters) =  my_final_clusters$names;  my_final_clusters$names = NULL; 
  }
}

######## Compute values for membership score #############
#Compute the total number of genes in clusters.
all_IDs = unique(c(as.vector(my_table$ID1), as.vector(my_table$ID2)))
TOT_genes_in_cluster = length(unique(sub(";.*", "", sub("\\|", ";", sub(".*;", "", sub("\\|", ";", all_IDs))))))

#Create a dictionary with exonID, number of TRUE Best_reciprocal
#I select the number of best reciprocals within the same exon cluster
#Create a dictionary with exonID -> number of species genes in the exon cluster.
best_reciprocal_matrix = vector()
gene_species_count_matrix = vector()
for (my_exon in rownames(my_final_clusters)) {
  my_cluster_id = my_final_clusters[my_exon,"ClusterID"]
  my_species = sub("\\|", ";", my_exon)
  all_exons_in_cluster = rownames(subset(my_final_clusters, ClusterID ==  my_cluster_id))
  #build best hits matrix
  all_rec_hits = unique(c(as.vector(subset(my_table, ID1==my_exon & Best_reciprocal == "TRUE")$ID2), as.vector(subset(my_table, ID2==my_exon & Best_reciprocal == "TRUE")$ID1)))
  rec_hits_in_cluster = length(all_rec_hits[all_rec_hits %in% all_exons_in_cluster])
  best_reciprocal_matrix = rbind(best_reciprocal_matrix, c(my_exon, rec_hits_in_cluster))
  #build species count matrix
  all_exons_in_species = all_exons_in_cluster[sub("\\|", ";", all_exons_in_cluster) == my_species] #select only the exons_in_cluster belonging to species
  gene_species_vector = unique(sub("\\|.*\\|", ";", sub(".*;", "", sub("\\|", ";", all_exons_in_species)))) #"ENSG00000103266;Hs2" #select all the genes containing the exons_in_cluster belonging to species
  cluster_gene_species_count_df = as.data.frame(table(sub(".*;", "", gene_species_vector))); colnames(cluster_gene_species_count_df) = c("Species", "Freq") 
  cluster_gene_species_count_df$ExonID = rep(my_exon, nrow(cluster_gene_species_count_df))
  gene_species_count_matrix = rbind(gene_species_count_matrix, cluster_gene_species_count_df[,c("ExonID", "Freq")])
}
best_reciprocal_df = as.data.frame(best_reciprocal_matrix); colnames(best_reciprocal_df) = c("ExonID", "Best_reciprocal_num")
best_reciprocal_df$ExonID = as.vector(best_reciprocal_df$ExonID); best_reciprocal_df$Best_reciprocal_num = as.numeric(as.vector(best_reciprocal_df$Best_reciprocal_num))
gene_species_count_df = as.data.frame(gene_species_count_matrix); gene_species_count_df$ExonID = as.vector(gene_species_count_df$ExonID)


####Compute the values necessary for the membership score.
my_sizes_df = as.data.frame(sizes(my_clusters)); colnames(my_sizes_df) = c("ClusterID", "Freq") #number of exons in each community

#Create subgraphs to count the in and out degree ONLY within the same exon cluster.
#Create a dataframe with ExonID, in_degree, out_degree.
in_degree_df = data.frame()
out_degree_df = data.frame()
for (my_cluster_id in unique(as.vector(my_final_clusters$ClusterID))) {
  all_cluster_exons = rownames(subset(my_final_clusters, ClusterID==my_cluster_id))
  my_subgraph = induced_subgraph(my_edges, v=all_cluster_exons) #subset the original graph
  my_in_degrees_df = as.data.frame(as.matrix(degree(my_subgraph, mode = "in"))) #direct IN connections of each exon IN THE SAME EXON CLUSTER
  my_out_degrees_df = as.data.frame(as.matrix(degree(my_subgraph, mode = "out"))) #direct OUT connections of each exon IN THE SAME EXON CLUSTER
  in_degree_df = rbind(in_degree_df, my_in_degrees_df)
  out_degree_df = rbind(out_degree_df, my_out_degrees_df)
}
in_degree_df$V1 = as.vector(in_degree_df$V1); out_degree_df$V1 = as.vector(out_degree_df$V1)

#Compute the number of exons for each species in each community
all_ids = paste0(gsub(".*\\|", "", rownames(my_final_clusters)), "_", my_final_clusters$ClusterID)
all_species_counts_df = as.data.frame(table(gsub(".*\\|", "", all_ids))); colnames(all_species_counts_df) = c("Species_ClusterID", "Freq"); 
all_species_counts_df$Species_ClusterID = as.vector(all_species_counts_df$Species_ClusterID) #transform factor to vector for hashmap to work

#Create hashes for translation
gene_species_count_dict = hashmap(gene_species_count_df$ExonID, gene_species_count_df$Freq)
best_reciprocal_dict = hashmap(best_reciprocal_df$ExonID, best_reciprocal_df$Best_reciprocal_num)
out_degree_dict = hashmap(rownames(out_degree_df), out_degree_df$V1)
in_degree_dict = hashmap(rownames(in_degree_df), in_degree_df$V1)
species_counts_dict = hashmap(all_species_counts_df$Species, all_species_counts_df$Freq)
sizes_dict = hashmap(my_sizes_df$ClusterID, my_sizes_df$Freq)

#Operate translations from hashes
my_final_clusters$ExonID = rownames(my_final_clusters)
print(head(my_final_clusters))
my_final_clusters$Out_degree = out_degree_dict$find(my_final_clusters$ExonID)
my_final_clusters$In_degree = in_degree_dict$find(my_final_clusters$ExonID)
my_final_clusters$SPECIES_exons_in_cluster = species_counts_dict$find(paste0(gsub(".*\\|", "", my_final_clusters$ExonID), "_", my_final_clusters$ClusterID))
my_final_clusters$TOT_exons_in_cluster = sizes_dict$find(my_final_clusters$ClusterID)
my_final_clusters$N_reciprocals = best_reciprocal_dict$find(my_final_clusters$ExonID)
my_final_clusters$SPECIES_genes_in_cluster = gene_species_count_dict$find(my_final_clusters$ExonID)

######## Compute membership score #############
my_final_clusters$membership_score = (my_final_clusters$In_degree + my_final_clusters$Out_degree + my_final_clusters$N_reciprocals)/(2*(my_final_clusters$TOT_exons_in_cluster - my_final_clusters$SPECIES_exons_in_cluster) + 
                                                                                                                                       (TOT_genes_in_cluster - my_final_clusters$SPECIES_genes_in_cluster))

my_final_clusters$membership_score = round(my_final_clusters$membership_score, 2)
#reorder columns
my_final_clusters = my_final_clusters[,c("ExonID", "ClusterID", "Out_degree", "In_degree", "SPECIES_exons_in_cluster", "TOT_exons_in_cluster", "N_reciprocals", "membership_score")]

#save to file
write.table(my_final_clusters, file=file2, quote=FALSE, sep="\t", row.names = FALSE, col.names=FALSE)
