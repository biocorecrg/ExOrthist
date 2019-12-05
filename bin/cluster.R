#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
if (length(args)<2) {stop("[USAGE] Rscript --vanilla cluster.R <file1> <file2> ")}
file1 <- args[1]
file2 <- args[2]

library(igraph)
e <- read.table(file1)
m <- as.matrix (e)
g2 <- graph_from_edgelist(m, directed = FALSE)
memb <- cluster_edge_betweenness(g2, weights =NULL)
e <- membership (memb)
write.table(as.matrix(e), file =file2, sep="\t")
