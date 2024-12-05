################################################
########## EXINT PLOTTER FUNCTIONS #############
################################################

#This Rscript will be source by the main one.

################################
##### Define functions #########
################################

position_unmatched_exons = function(unmatched_target, my_gene_target_table, all_matched_exons_indexes) {
  #All these differences are based on indexes.
  unmatched_target_index = match(unmatched_target, as.vector(my_gene_target_table$ExonID)) #get the index of the unmatched exon in the gene.
  differences = unmatched_target_index - all_matched_exons_indexes
  
  #all the comments refer to the case where the strand is positive.
  min_positive_difference = min(differences[differences > 0]) #directly upstream exon (difference in indexes)
  min_positive_difference_index = match(min_positive_difference, differences) #select the index of the upstream exon
  min_negative_difference = max(differences[differences < 0]) #maximum negative difference it's the minimum negative difference in absolute values
  min_negative_difference_index  =  match(min_negative_difference, differences)  #select the index of the downstream exon
  
  #the selection is different depending if the missing exons are on the positive or negative strand
  if (unique(as.vector(my_gene_target_table$Strand)) == "+") {
    upstream_exon = all_matched_exons[min_positive_difference_index]  #select the directly upstream exon
    downstream_exon = all_matched_exons[min_negative_difference_index]  #directly downstream exon
  }
  if (unique(as.vector(my_gene_target_table$Strand)) == "-") {
    upstream_exon = all_matched_exons[min_negative_difference_index]  #select the directly upstream exon
    downstream_exon = all_matched_exons[min_positive_difference_index]  #directly downstream exon
  }
  if (!(is.na(upstream_exon))) { #if there is an upstream exon.
    my_gene_target_table$ClusterID[my_gene_target_table$ExonID == unmatched_target] = paste0(upstream_exon, ";", downstream_exon)
    my_FakeStart = my_gene_target_table$FakeStop[my_gene_target_table$ExonID == upstream_exon] + 5 #stop of the previous exon + 5
    my_FakeStop = my_FakeStart+5 #These exons will only have length 5.
  }
  if (is.na(upstream_exon) & !(is.na(downstream_exon))) { #if there is no upstream exons but there is a downstream exon (exons at the end of the gene)
    my_gene_target_table$ClusterID[my_gene_target_table$ExonID == unmatched_target] = paste0("NA;", downstream_exon)
    my_FakeStart = -10
    my_FakeStop = -5
  }
  my_gene_target_table$FakeStart[my_gene_target_table$ExonID == unmatched_target] = my_FakeStart
  my_gene_target_table$FakeStop[my_gene_target_table$ExonID == unmatched_target] = my_FakeStop
  my_gene_target_table$FakeCoords[my_gene_target_table$ExonID == unmatched_target] = paste0(my_FakeStart,";", my_FakeStop)
  return(my_gene_target_table)
}

########
generate_NA_coordinates = function(unmatched_target, my_gene_target_table)  {
  my_gene_target_table$FakeStart[my_gene_target_table$ExonID == unmatched_target] = NA
  my_gene_target_table$FakeStop[my_gene_target_table$ExonID == unmatched_target] = NA
  my_gene_target_table$FakeCoords[my_gene_target_table$ExonID == unmatched_target] = "NA;NA"
  return(my_gene_target_table)
}

##########
get_last_ex_df = function(plotting_table) {
  last_ex_df = subset(plotting_table, ExPosition=="last")
  last_ex_df$first_point = paste0(last_ex_df$FakeStart, ";", 
                                   as.character(last_ex_df$Order+0.25))
  last_ex_df$second_point = paste0(as.character(last_ex_df$FakeStop), 
                                   ";", as.character(last_ex_df$Order+0.5))
  last_ex_df$third_point = paste0(as.character(last_ex_df$FakeStop), 
                                    ";", as.character(last_ex_df$Order))

  last_ex_df$my_ID = paste0(last_ex_df$ExonID,  ";", last_ex_df$State, ";", last_ex_df$Filling_status, ";", last_ex_df$AnnotStatus, ";", last_ex_df$IsoformExs)
  last_ex_df = last_ex_df[,c("my_ID", "first_point", "second_point",  "third_point")]
  last_ex_df_long = melt(last_ex_df, id="my_ID")
  last_ex_df_long$ExonID = sapply(strsplit(last_ex_df_long $my_ID,";"), `[`, 1)
  last_ex_df_long$State = sapply(strsplit(last_ex_df_long $my_ID,";"), `[`, 2)
  last_ex_df_long$Filling_status = sapply(strsplit(last_ex_df_long $my_ID,";"), `[`, 3)
  last_ex_df_long$AnnotStatus = sapply(strsplit(last_ex_df_long $my_ID,";"), `[`, 4)
  last_ex_df_long$IsoformExs = sapply(strsplit(last_ex_df_long $my_ID,";"), `[`, 5)
  #last_ex_df_long$ExonID = sub(";.*", "", last_ex_df_long$my_ID)
  #last_ex_df_long$State  = sub(".*;", "", last_ex_df_long$my_ID)
  last_ex_df_long$x = as.numeric(sub(";.*", "", last_ex_df_long$value))
  last_ex_df_long$y = as.numeric(sub(".*;", "", last_ex_df_long$value))
  last_ex_df_long$my_ID =  NULL
  last_ex_df_long$variable = NULL
  last_ex_df_long$value = NULL
  return(last_ex_df_long)
}

get_first_ex_df = function(plotting_table) {
  first_ex_df = subset(plotting_table, ExPosition=="first")
  first_ex_df$first_point = paste0(first_ex_df$FakeStart, ";", 
                                  as.character(first_ex_df$Order))
  first_ex_df$second_point = paste0(first_ex_df$FakeStart, ";", 
                                   as.character(first_ex_df$Order+0.5))

  first_ex_df$third_point = paste0(as.character(first_ex_df$FakeStop), 
                                   ";", as.character(first_ex_df$Order+0.25))
  
  first_ex_df$my_ID = paste0(first_ex_df$ExonID,  ";", first_ex_df$State, ";", first_ex_df$Filling_status, ";", first_ex_df$AnnotStatus, ";", first_ex_df$IsoformExs)
  first_ex_df = first_ex_df[,c("my_ID", "first_point", "second_point", "third_point")]
  first_ex_df_long = melt(first_ex_df, id="my_ID")
  first_ex_df_long$ExonID = sapply(strsplit(first_ex_df_long$my_ID,";"), `[`, 1)
  first_ex_df_long$State = sapply(strsplit(first_ex_df_long$my_ID,";"), `[`, 2)
  first_ex_df_long$Filling_status = sapply(strsplit(first_ex_df_long$my_ID,";"), `[`, 3)
  first_ex_df_long$AnnotStatus = sapply(strsplit(first_ex_df_long$my_ID,";"), `[`, 4)
  first_ex_df_long$IsoformExs = sapply(strsplit(first_ex_df_long$my_ID,";"), `[`, 5)
  #first_ex_df_long$ExonID = sub(";.*", "", first_ex_df_long$my_ID)
  #first_ex_df_long$State  = sub(".*;", "", first_ex_df_long$my_ID)
  first_ex_df_long$x = as.numeric(sub(";.*", "", first_ex_df_long$value))
  first_ex_df_long$y = as.numeric(sub(".*;", "", first_ex_df_long$value))
  first_ex_df_long$my_ID =  NULL
  first_ex_df_long$variable = NULL
  first_ex_df_long$value = NULL
  return(first_ex_df_long)
}

get_first_last_ex_df = function(plotting_table) {
  first_last_ex_df = subset(plotting_table, ExPosition=="first;last")
  first_last_ex_df$first_point = paste0(first_last_ex_df$FakeStart, ";", 
                                        as.character(first_last_ex_df$Order+0.25))
  first_last_ex_df$second_point = paste0(as.character(first_last_ex_df$FakeStart + (abs(first_last_ex_df$FakeStop - first_last_ex_df$FakeStart))/2),
                                         ";", as.character(first_last_ex_df$Order+0.5))
  first_last_ex_df$third_point = paste0(as.character(first_last_ex_df$FakeStop), 
                                        ";", as.character(first_last_ex_df$Order+0.25))
  first_last_ex_df$fourth_point = paste0(as.character(first_last_ex_df$FakeStart + (abs(first_last_ex_df$FakeStop - first_last_ex_df$FakeStart))/2),
                                         ";", as.character(first_last_ex_df$Order))
  first_last_ex_df$my_ID = paste0(first_last_ex_df$ExonID,  ";", first_last_ex_df$State, ";", first_last_ex_df$Filling_status, ";", first_last_ex_df$AnnotStatus, ";", first_last_ex_df$IsoformExs)
  first_last_ex_df = first_last_ex_df[,c("my_ID", "first_point", "second_point", "third_point", "fourth_point")]
  first_last_ex_df_long = melt(first_last_ex_df, id="my_ID")
  first_last_ex_df_long$ExonID = sapply(strsplit(first_last_ex_df_long $my_ID,";"), `[`, 1)
  first_last_ex_df_long$State = sapply(strsplit(first_last_ex_df_long $my_ID,";"), `[`, 2)
  first_last_ex_df_long$Filling_status = sapply(strsplit(first_last_ex_df_long $my_ID,";"), `[`, 3)
  first_last_ex_df_long$AnnotStatus = sapply(strsplit(first_last_ex_df_long $my_ID,";"), `[`, 4)
  first_last_ex_df_long$IsoformExs = sapply(strsplit(first_last_ex_df_long $my_ID,";"), `[`, 5)
  #first_last_ex_df_long$ExonID = sub(";.*", "", first_last_ex_df_long$my_ID)
  #first_last_ex_df_long$State  = sub(".*;", "", first_last_ex_df_long$my_ID)
  first_last_ex_df_long$x = as.numeric(sub(";.*", "", first_last_ex_df_long$value))
  first_last_ex_df_long$y = as.numeric(sub(".*;", "", first_last_ex_df_long$value))
  first_last_ex_df_long$my_ID =  NULL
  first_last_ex_df_long$variable = NULL
  first_last_ex_df_long$value = NULL
  return(first_last_ex_df_long)
}