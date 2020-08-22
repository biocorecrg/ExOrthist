#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
if (length(args)<2) {stop("[USAGE] Rscript --vanilla cluster.R <file1> <file2> ")}
file1 <- args[1]
file2 <- args[2]

library(igraph, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(hashmap, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)

my_table = read.table(file1)
my_matrix = as.matrix(my_table)
my_edges = graph_from_edgelist(my_matrix, directed=TRUE) #create directed graph starting from input.
my_clusters = cluster_edge_betweenness(my_edges) #graph with the best modularity score.
my_final_clusters = as.data.frame(as.matrix(membership(my_clusters))); colnames(my_final_clusters) = "ClusterID" #communities

#compute degrees necessary for the membership score.
my_sizes_df = as.data.frame(sizes(my_clusters)); colnames(my_sizes_df) = c("ClusterID", "Freq") #number of exons in each community
my_out_degrees_df = as.data.frame(as.matrix(degree(my_edges, mode = "out"))) #direct connections of each exon
my_in_degrees_df = as.data.frame(as.matrix(degree(my_edges, mode = "in"))) #reciprocal hits
#compute the number of exons for each species in each community
all_ids = paste0(gsub(".*\\|", "", rownames(my_final_clusters)), "_", my_final_clusters$ClusterID)
all_species_counts_df = as.data.frame(table(gsub(".*\\|", "", all_ids))); colnames(all_species_counts_df) = c("Species_ClusterID", "Freq"); 
all_species_counts_df$Species_ClusterID = as.vector(all_species_counts_df$Species_ClusterID) #transform factor to vector for hashmap to work

#create hashes for translation
out_degree_dict = hashmap(rownames(my_out_degrees_df), my_out_degrees_df$V1)
in_degree_dict = hashmap(rownames(my_in_degrees_df), my_in_degrees_df$V1)
species_counts_dict = hashmap(all_species_counts_df$Species, all_species_counts_df$Freq)
sizes_dict = hashmap(my_sizes_df$ClusterID, my_sizes_df$Freq)

#compute membership score and save to final file.
my_final_clusters$out_degree = out_degree_dict$find(rownames(my_final_clusters))
my_final_clusters$in_degree = in_degree_dict$find(rownames(my_final_clusters))
my_final_clusters$SPECIES_exons_in_cluster = species_counts_dict$find(paste0(gsub(".*\\|", "", rownames(my_final_clusters)), "_", my_final_clusters$ClusterID))
my_final_clusters$TOT_exons_in_cluster = sizes_dict$find(my_final_clusters$ClusterID)
my_final_clusters$membership_score = (my_final_clusters$in_degree + my_final_clusters$out_degree)/(2*(my_final_clusters$TOT_exons_in_cluster - my_final_clusters$SPECIES_exons_in_cluster))
my_final_clusters$ExonID = rownames(my_final_clusters)

#reorder columns
my_final_clusters = my_final_clusters[,c("ExonID", "ClusterID", "out_degree", "in_degree", "SPECIES_exons_in_cluster", "TOT_exons_in_cluster", "membership_score")]

write.table(my_final_clusters, file=file2, quote=FALSE, sep="\t", row.names = FALSE)