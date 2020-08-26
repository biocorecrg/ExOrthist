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
my_clusters = cluster_edge_betweenness(my_edges) #Modularity does not work with directed graphs.
my_final_clusters = as.data.frame(as.matrix(membership(my_clusters))); colnames(my_final_clusters) = "ClusterID" #communities

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

#Create a dictionary with species -> number of species genes in the cluster.
#Create a dictionary with exonID -> number of species genes in the exon cluster.
# gene_species_vector = unique(sub("\\|.*\\|", ";", sub(".*;", "", sub("\\|", ";", all_IDs))))
# gene_species_count_df = as.data.frame(table(sub(".*;", "", gene_species_vector))); colnames(gene_species_count_df) = c("Species", "Freq")
# gene_species_count_df$Species = as.vector(gene_species_count_df$Species)

####Compute the values necessary for the membership score.
my_sizes_df = as.data.frame(sizes(my_clusters)); colnames(my_sizes_df) = c("ClusterID", "Freq") #number of exons in each community
my_out_degrees_df = as.data.frame(as.matrix(degree(my_edges, mode = "out"))) #direct connections of each exon
my_in_degrees_df = as.data.frame(as.matrix(degree(my_edges, mode = "in"))) #reciprocal hits
#Compute the number of exons for each species in each community
all_ids = paste0(gsub(".*\\|", "", rownames(my_final_clusters)), "_", my_final_clusters$ClusterID)
all_species_counts_df = as.data.frame(table(gsub(".*\\|", "", all_ids))); colnames(all_species_counts_df) = c("Species_ClusterID", "Freq"); 
all_species_counts_df$Species_ClusterID = as.vector(all_species_counts_df$Species_ClusterID) #transform factor to vector for hashmap to work

#Create hashes for translation
#gene_species_count_dict = hashmap(gene_species_count_df$Species, gene_species_count_df$Freq)
gene_species_count_dict = hashmap(gene_species_count_df$ExonID, gene_species_count_df$Freq)
best_reciprocal_dict = hashmap(best_reciprocal_df$ExonID, best_reciprocal_df$Best_reciprocal_num)
out_degree_dict = hashmap(rownames(my_out_degrees_df), my_out_degrees_df$V1)
in_degree_dict = hashmap(rownames(my_in_degrees_df), my_in_degrees_df$V1)
species_counts_dict = hashmap(all_species_counts_df$Species, all_species_counts_df$Freq)
sizes_dict = hashmap(my_sizes_df$ClusterID, my_sizes_df$Freq)

#Operate translations from hashes
my_final_clusters$ExonID = rownames(my_final_clusters)
my_final_clusters$Out_degree = out_degree_dict$find(my_final_clusters$ExonID)
my_final_clusters$In_degree = in_degree_dict$find(my_final_clusters$ExonID)
my_final_clusters$SPECIES_exons_in_cluster = species_counts_dict$find(paste0(gsub(".*\\|", "", my_final_clusters$ExonID), "_", my_final_clusters$ClusterID))
my_final_clusters$TOT_exons_in_cluster = sizes_dict$find(my_final_clusters$ClusterID)
my_final_clusters$N_reciprocals = best_reciprocal_dict$find(my_final_clusters$ExonID)
#my_final_clusters$SPECIES_genes_in_cluster = gene_species_count_dict$find(sub(".*\\|", "", my_final_clusters$ExonID))
my_final_clusters$SPECIES_genes_in_cluster = gene_species_count_dict$find(my_final_clusters$ExonID)

######## Compute membership score #############
my_final_clusters$membership_score = (my_final_clusters$In_degree + my_final_clusters$Out_degree + my_final_clusters$N_reciprocals)/(2*(my_final_clusters$TOT_exons_in_cluster - my_final_clusters$SPECIES_exons_in_cluster) + 
                                                                                                                                       (TOT_genes_in_cluster - my_final_clusters$SPECIES_genes_in_cluster))
#reorder columns
my_final_clusters = my_final_clusters[,c("ExonID", "ClusterID", "Out_degree", "In_degree", "SPECIES_exons_in_cluster", "TOT_exons_in_cluster", "N_reciprocals", "membership_score")]

#save to file
write.table(my_final_clusters, file=file2, quote=FALSE, sep="\t", row.names = FALSE, col.names=FALSE)
