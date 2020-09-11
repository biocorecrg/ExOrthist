#######################################
########## EXINT PLOTTER  #############
#######################################

#The script will take as argument:
#1. The geneID
#2. The output folder

### The reference species can be derived from the geneID.
#This is just for now, for the debugging. This part will have to be modified accordingly.
#.libPaths(c("/Users/federica/mnt/software/R/3.5", "/Users/federica/mnt/software/R", .libPaths()))
#.libPaths(c("/software/mi/Rlib3.6/", "/software/as/el7.2/EasyBuild/CRG/software/R/3.6.0-foss-2019a/lib64/R/library", .libPaths()))

#upload libraries
suppressWarnings(library("ggplot2")) 
suppressWarnings(library("ggpubr"))
suppressWarnings(library("cowplot"))
suppressWarnings(library("reshape2"))
suppressWarnings(library("hashmap"))


#How it will be when reading the arguments
args <- commandArgs(TRUE)
my_gene = args[1]
my_query_species = args[2]
my_folder = args[3]
my_script_folder = args[4]
gene_clusters = args[5]
my_ordered_species_raw = args[6] #The order of the species is defined within the nextflow.
my_isoform_exons_raw = args[7]
my_isorform_id = as.vector(unlist(strsplit(my_isoform_exons_raw, ",")))[1]
if (nchar(my_isorform_id) == 0) {my_isorform_id = "None"}
my_isoform_exons = as.vector(unlist(strsplit(my_isoform_exons_raw, ",")))[2:length(as.vector(unlist(strsplit(my_isoform_exons_raw, ","))))]

if (length(args)==8) {
  my_interesting_exons = args[8]
  interesting_exons = as.vector(read.table(my_interesting_exons, header=FALSE)$V1)
  } else {
    interesting_exons = ""
    }


######## Set up
my_input_folder = paste0(my_folder, "/processed_tables/")
ordered_target_species = unlist(strsplit(my_ordered_species_raw, split=","))
source(paste0(my_script_folder, "/exint_plotter_functions.R"))

###### for testing
#my_gene="ENSG00000006788"
# my_gene="ENSG00000152601"
# my_query_species="Hs2"
# my_folder = paste0("/Users/federica/mnt/projects/regulatory_ancenstry/src/nextflow/EXORTHIST/ExOrthist/exint_plotter/output_exint/", my_gene, "/")
# #gene_clusters = "/Users/federica/mnt/projects/regulatory_ancenstry/src/nextflow/EXORTHIST/ExOrthist/exint_plotter/data/clusters/gene_clusters"
# gene_clusters = "/Users/federica/mnt/projects/regulatory_ancenstry/src/nextflow/EXORTHIST/ExOrthist/output_prot_sim/gene_cluster_file.gz"
# my_input_folder = paste0(my_folder, "/processed_tables/")
# #my_ordered_species_raw="Hs2,Mm2,Bta"
# my_ordered_species_raw="Hs2,Mm2,Dme"
# ordered_target_species = unlist(strsplit(my_ordered_species_raw, split=","))
# interesting_exons = ""
# my_isoform_exons = "None"

# interesting_exons = as.vector(read.table("/Users/federica/mnt/projects/regulatory_ancenstry/src/nextflow/EXORTHIST/ExOrthist/exint_plotter/data/clusters/interesting_exons.txt", header=FALSE)$V1)
# my_isoform_exons = c("chr17:10362360-10362503", "chr17:10362118-10362274", "chr17:10360161-10360188")




######## Main code ##########
########### QUERY SPECIES SECTION ################
#I need to find a way to define the species_query and the order of the other_species starting from the geneID.
#provide order of the species that you want to plot.
#In principle they should be ordered by evolutionary distance.

########
species_query_table = read.delim(paste0(my_input_folder, my_query_species, "_exons_cluster_info-fakecoords.tab"), sep="\t", header=TRUE, row=1)
species_query_table = species_query_table[species_query_table$State!="",] #removing some weird spurious entries

my_gene_table = subset(species_query_table, GeneID == my_gene) #select the gene of interest
my_gene_table$FakeCoords = as.character(paste0(my_gene_table$FakeStart, ";", my_gene_table$FakeStop)) #generate coordinates ID
query_coords_dict = hashmap(as.character(my_gene_table$ClusterID), my_gene_table$FakeCoords) #dictionary with key=cluster ID, value =  pair of fake coordinates.
raw_clusterID = unique(sub("\\..*", "", as.vector(my_gene_table$ClusterID))) #Remove the ".???" suffix (relative to the exon cluster) from the gene clusterID.
my_gene_clusterID = raw_clusterID[grep("^GF", raw_clusterID)] #this I will use to upload the best-hits table by gene_cluster. It will be much faster.

my_pos_table = read.delim(paste0(my_input_folder, my_query_species, "_all_exons_positions.tab"),
                          header=FALSE,
                          col.names=c("GeneID", "GeneClusterID", "ExonID", "position"), 
                          stringsAsFactors = FALSE) #upload table with all exons index by cluster.
my_exon_position_dict = hashmap(my_pos_table$ExonID, my_pos_table$position) #dictionary with exonID, exon index.
my_gene_table$ExPosition = my_exon_position_dict$find(my_gene_table$ExonID) #add the exon position the the query_species table.
my_gene_table$ExPosition[is.na(my_gene_table$ExPosition)] = "Internal" #whatever exon without label we keep as Internal.

####### ORTHOLOGS #########
#Get orthologs either from the general clusters or from the files provided by the user.
species_orthologs_table = read.delim(gene_clusters,
                                     sep="\t", header=FALSE)
colnames(species_orthologs_table)[1:4] = c("ClusterID", "Species", "GeneID", "GeneName")

my_clusterID = unique(as.vector(subset(species_orthologs_table, GeneID == my_gene)$ClusterID)) #isolate clusterID
my_gene_orthologs = as.vector(subset(species_orthologs_table, ClusterID == my_clusterID & Species != my_query_species)$GeneID) #generate named vector with all the one2one orthologs of the gene of interest.
names(my_gene_orthologs) = as.vector(subset(species_orthologs_table, ClusterID == my_clusterID & Species != my_query_species)$Species) #use species as vector names

###########################
#Generate final table
final_table = subset(my_gene_table, State != "Intron") #This could be simplified in the input. I do not really need the intron information so far.
my_gene_table = final_table #This is weird
final_table$Levels = rep(1, nrow(final_table))
final_table$Species = rep(my_query_species, nrow(final_table))
#Get exon coordinates from ExonID. 
my_stop_coords = sub(".*-", "", final_table$ExonID)
my_start_coords = sub(".*:", "", sub("-.*", "", final_table$ExonID))
final_table$ExonLength = as.numeric(my_stop_coords) - as.numeric(my_start_coords)

########### TARGET SPECIES SECTION #################
considered_species = ordered_target_species[(ordered_target_species  %in% names(my_gene_orthologs))] #select only the species where there is an ortholog

#problem: when there are only exons which are not in the exon clusters.
for (my_target_species in considered_species) {
  #upload target species table wihout fake coordinates
  target_species_table = read.delim(paste0(my_input_folder, my_target_species, "_exons_cluster_info-fakecoords.tab"), 
                                    sep="\t", header=TRUE, row=1, stringsAsFactors=FALSE)
  target_species_table = target_species_table[target_species_table$State!="",] #removing some weird spurious entries
  
  #cycle on all the orthologoues genes in target species
  for (target_gene in my_gene_orthologs[names(my_gene_orthologs) == my_target_species]) {
    my_gene_target_table_complete = subset(target_species_table, GeneID==target_gene)
    my_gene_target_table_complete$FakeCoords = query_coords_dict$find(my_gene_target_table_complete$ClusterID) #translate GeneClusterID with the fake coords in query species.
    
    #temporarily remove the Introns from the gene_target_table
    my_gene_target_table =  subset(my_gene_target_table_complete, State != "Intron") #I should not add the Introns in the first place.
    #Get exon coordinates starting from the ExonID
    my_stop_coords = sub(".*-", "", my_gene_target_table$ExonID)
    my_start_coords = sub(".*:", "", sub("-.*", "", my_gene_target_table$ExonID))
    my_gene_target_table$ExonLength = as.numeric(my_stop_coords) - as.numeric(my_start_coords)
    
    #upload best hits table
    best_hits_table = read.delim(paste0(my_input_folder, my_query_species, "_", my_target_species, "-best_scores_hits_exons.translated_coords"), header=TRUE)
    
    ### remove useless columns
    best_hits_table$ExonID_query = rep(NA, nrow(best_hits_table))
    best_hits_table$ExonID_target = rep(NA, nrow(best_hits_table))
    #This in theory is just to remove columns
    best_hits_table = subset(best_hits_table, select=-c(Score_C1, Score_I1, Score_I2, Score_C2))
    best_hits_table = unique(best_hits_table) #Check this. IDK why we should have duplicated rows.
    
    
    ##################################################################################################################
    ####### Consider the case where more exons in query species belong to the same exon cluster ######################
    ##################################################################################################################
    
    #This part I actually still have to code. Within each cluster, I should assign the one with the best total score.
    
    ##################################################################################################################
    ####### Consider the case where in query species there is an exon not in the exon clusters #######################
    #Consider also the cases where the exon in species query is in the exon clusters but none of the species target exons are
    ##################################################################################################################
    
    #belonging to an exon cluster whose cluster is not represented in query species || not belonging to an exon cluster
    not_conserved_query = unique(as.vector(my_gene_table$ExonID[!(my_gene_table$ClusterID %in% my_gene_target_table$ClusterID)])) 
    all_best_hits_df = vector()
    id_to_keep = as.vector(subset(my_gene_target_table, !(grepl(";", FakeCoords)))$ExonID) #select only the target exons still not matched by exon cluster.
    
    #all of this makes sense only if there are still unmatched exons both in query and unmatched exons in target.
    if (length(not_conserved_query) > 0 & length(id_to_keep) > 0) {
      for (unmatched_exon in not_conserved_query) { #for each unmatched exon
        sequence_hits_df = subset(best_hits_table, OverlappingID_query==unmatched_exon & GeneID_query==my_gene & GeneID_target==target_gene) #subset the table for the query unmatched exon
        sequence_hits_df = subset(sequence_hits_df, OverlappingID_target %in% id_to_keep) #select only the target exons which are still unmatched.
        
        if (nrow(sequence_hits_df) > 1) {#if still there are more than one hits, select the best based on Total score.
          max_score = max(as.vector(sequence_hits_df$Total_exon_score)) #select the max Total score between the IDs present.
          sequence_hits_df = subset(sequence_hits_df, Total_exon_score == max_score) #subset by the max_score.
          my_entry = c(as.vector(sequence_hits_df$OverlappingID_target), unmatched_exon, as.numeric(sequence_hits_df$Total_exon_score)) #get the entry for the final df.
          all_best_hits_df = rbind(all_best_hits_df, my_entry) #this will be in the format:ExonID_target, ExonID_query, Total_score.
        } else if (nrow(sequence_hits_df) == 1) {
          my_entry = c(as.vector(sequence_hits_df$OverlappingID_target), unmatched_exon, as.numeric(sequence_hits_df$Total_exon_score))
          all_best_hits_df = rbind(all_best_hits_df, my_entry)
        }
      }
      
      all_best_hits_df = as.data.frame(all_best_hits_df) #make sure it is a dataframe.
      
      
      if (nrow(all_best_hits_df) > 0) { 
        colnames(all_best_hits_df) = c("target_ID", "query_ID", "Total_score")
        #In case of duplicated best hits (same target exon being the best hit for more query exons), keep the one with the highest total score. 
        #This means that only the best query exon will have a match.
        duplicated_target_exons = unique(as.vector(all_best_hits_df$target_ID[duplicated(as.vector(all_best_hits_df$target_ID))]))
        if (length(duplicated_target_exons) > 0 ) {
          to_process_df = all_best_hits_df
          all_best_hits_df = subset(all_best_hits_df, !(target_ID %in% duplicated_target_exons)) #remove the duplicated target exons from the all_best_hits_df.
          for (target_exon in duplicated_target_exons) {
            max_score = max(as.vector(subset(to_process_df, target_ID  == target_exon)$Total_score)) #compute the max Total Score between all query exons matching the same target.
            to_keep_df = subset(to_process_df, target_ID == target_exon & Total_score == max_score) #select only the entries where the target_ID is different from target_exon AND the Total_score is the max.
            all_best_hits_df = rbind(all_best_hits_df, to_keep_df)
          }
        }
        #cycle on query exons and translate the relative target coordinates.
        #This was outside the loop, but I think it should be here
        for (element1 in all_best_hits_df$query_ID) {
          target_ID = as.character(all_best_hits_df$target_ID[all_best_hits_df$query_ID == element1])
          my_gene_target_table$FakeCoords[my_gene_target_table$ExonID == target_ID] = as.vector(my_gene_table$FakeCoords[my_gene_table$ExonID == element1])
          my_gene_target_table$State[my_gene_target_table$ExonID == target_ID] = "Exon_added"
        }
      }
    }
    #reorder gene target table so that the coordinates are ordered. This will be necessary for the exons in the next section to be inserted in the right position.
    my_gene_target_table$FakeStart = as.numeric(sub(";.*", "", my_gene_target_table$FakeCoords)) 
    my_gene_target_table$FakeStop = as.numeric(sub(".*;", "", my_gene_target_table$FakeCoords))
    #I think this line is useless, I am commenting it out. It will have to be ordered by actual coords.
    #my_gene_target_table = my_gene_target_table[order(my_gene_target_table$FakeStart),] #order by FakeStart coordinates.
    
    
    ##################################################################################################################
    ####### Consider the case where in target species still remain exons which are not in the exon clusters. #########
    ##################################################################################################################
    
    #consider also the case where in the target species an exon belong to the exon clusters but to none of the clusters in query species
    my_gene_target_table$ClusterID = as.vector(my_gene_target_table$ClusterID) #isolate all the exon cluster IDs.
    my_gene_target_table$Start = sub(".*:", "", sub("-.*", "", my_gene_target_table$ExonID)) #add the start to ordered the table (From the ExonID)
    my_gene_target_table = my_gene_target_table[order(my_gene_target_table$Start),] #order the table by start coordinate.
    
    not_conserved_target = as.vector(subset(my_gene_target_table, !(grepl(";", FakeCoords)))$ExonID) #This are the ones that after all still remain unmatched.
    all_matched_exons = as.vector(subset(my_gene_target_table, grepl(";", FakeCoords))$ExonID) #select all matched exons.
    all_matched_exons_indexes = match(all_matched_exons, as.vector(my_gene_target_table$ExonID)) #get the indexes of all matched exons.
    
    if (length(all_matched_exons) == 0) {
      for (unmatched_target in not_conserved_target) {
        my_gene_target_table = generate_NA_coordinates(unmatched_target,  my_gene_target_table) #Just give NAs, because there's no exon with fakeCoords to take as reference.
      }
    }
    
    if (length(all_matched_exons) > 0 & length(not_conserved_target > 0)) { #if there are still unique exons in target.
      for (unmatched_target in not_conserved_target) {
        #The unmatched exon gets in the right position looking at upstream and downstream matched exons.
        my_gene_target_table = position_unmatched_exons(unmatched_target, my_gene_target_table, all_matched_exons_indexes) #function in the functions section.
      }
    }
    
    ##### Final modifications
    #dictionary with key=exon cluster ID, values=number of exons belonging to that cluster (same fake coords)
    my_levels_dict = hashmap(names(table(as.vector(my_gene_target_table$ClusterID))), unname(table(as.vector(my_gene_target_table$ClusterID))))
    my_gene_target_table$Levels = my_levels_dict$find(my_gene_target_table$ClusterID)
    my_gene_target_table$Species = rep(my_target_species, nrow(my_gene_target_table)) #add species.
    my_gene_target_table$Start = NULL #remove start
    
    ###### work on the  exon position
    #position means first, last, internal
    my_target_pos_table = read.delim(paste0(my_input_folder, my_target_species, "_all_exons_positions.tab"),
                                     header=FALSE,
                                     col.names=c("GeneID", "GeneClusterID", "ExonID", "position"), 
                                     stringsAsFactors = FALSE)
    my_target_exon_position_dict = hashmap(my_target_pos_table$ExonID, my_target_pos_table$position)
    my_gene_target_table$ExPosition = my_target_exon_position_dict$find(my_gene_target_table$ExonID) #Add position using ExonID as key
    my_gene_target_table$ExPosition[is.na(my_gene_target_table$ExPosition)] = "Internal" #whatever is not classified we will plot as a normal rectangle (internal).
    
    #Append to final table.
    final_table = rbind(final_table, my_gene_target_table)
  }
}


###########################################################
########### Prepare plotting input ########################
###########################################################
plotting_table = final_table
plotting_table$GeneID = as.vector(plotting_table$GeneID)
plotting_table$Species = as.vector(plotting_table$Species)

#########
#in case there are more genes per species
my_ordered_genes = unique(as.vector(plotting_table$GeneID))
my_genes_order_dict = hashmap(rev(my_ordered_genes), seq(1, length(my_ordered_genes))) #create ordered index for plotting.
plotting_table$Order = my_genes_order_dict$find(plotting_table$GeneID)
#########
gene_exon_number_dict = hashmap(names(table(plotting_table$GeneID)), unname(table(plotting_table$GeneID)))
plotting_table$ExonNumber = gene_exon_number_dict$find(plotting_table$GeneID) #total number of exons

###### Setting up the phases
#upstream phases
plotting_table$FinalPhaseUp = plotting_table$UpPhase #insert phase to print
plotting_table$FinalPhaseUp[plotting_table$Strand  ==  "-"] = plotting_table$DownPhase[plotting_table$Strand  ==  "-"] #invert the phases in case the strand is negative.
plotting_table$FinalPhaseUp[plotting_table$FinalPhaseUp == "last"] = NA
plotting_table$FinalPhaseUp = gsub("\\..*", "", as.character(plotting_table$FinalPhaseUp))
plotting_table$FinalPhaseUp[plotting_table$Levels > 1] = NA #set phase to NA when there are more exons in the same intron.

#downstream phases
plotting_table$FinalPhaseDown = plotting_table$DownPhase #insert phase to print
plotting_table$FinalPhaseDown[plotting_table$Strand  ==  "-"] = as.character(plotting_table$UpPhase)[plotting_table$Strand  ==  "-"]
plotting_table$FinalPhaseDown[plotting_table$FinalPhaseDown == "last"] = NA
plotting_table$FinalPhaseDown = gsub("\\..*", "", as.character(plotting_table$FinalPhaseDown))
plotting_table$FinalPhaseDown[plotting_table$Levels > 1] = NA #set phase to NA when there are more exons in the same intron. (or more exons in the same exon cluster).
plotting_table$Levels[plotting_table$Levels == 1] = NA #If there is only one exon, set Levels to NA (we do not need to plot the number in that case).

plotting_table$ExonNumberPlot = paste0(plotting_table$GeneID, ";", plotting_table$ExonNumber)
plotting_table$ExonNumberPlot = factor(plotting_table$ExonNumberPlot, levels=unique(plotting_table$ExonNumberPlot))

####### Setting up the colors.
#Adding special colors for the exons to highlight and their orthologs.
all_colors_vector = terrain.colors(length(interesting_exons)+1)[1:length(interesting_exons)] #generate as many colors as interesting exons
plotting_table$Filling_status = rep("default", nrow(plotting_table)) #set default filling for all exons
#group_colors_vector = c("default"="lightsteelblue3")
group_colors_vector = c("default"="gray66")
index=1
for (my_exon in interesting_exons) {
  fake_coords = plotting_table$FakeCoords[plotting_table$ExonID == my_exon] #isolate the coords of each relevant ex.
  groupID = paste0("group", index) #create a groupID
  plotting_table$Filling_status[plotting_table$FakeCoords == fake_coords] = groupID #replace "default" by groupID for all the exons in the same coordinates
  group_colors_vector = c(group_colors_vector, all_colors_vector[index]); names(group_colors_vector)[length(group_colors_vector)] = groupID #update color vector for plotting.
  index=index+1
}

#Adding colors to highlight the exons in the query gene which belong to a given isoform
plotting_table$IsoformExs = rep("black", nrow(plotting_table))
plotting_table$IsoformExs[plotting_table$ExonID %in% my_isoform_exons] = "brown2"

#All these functions are in the functions script
internal_ex_df = subset(plotting_table, ExPosition=="Internal")
if ("first" %in% unique(as.vector(plotting_table$ExPosition))) {
  first_ex_df = get_first_ex_df(plotting_table) } else {
    first_ex_df = data.frame(matrix(ncol = 7, nrow = 0)); colnames(first_ex_df) = c("x", "y", "ExonID", "State", "Filling_status", "AnnotStatus", "IsoformExs")} #create empty dataframe.
if("last" %in% unique(as.vector(plotting_table$ExPosition))) {
  last_ex_df = get_last_ex_df(plotting_table)} else {
    last_ex_df = data.frame(matrix(ncol = 7, nrow = 0)); colnames(last_ex_df) = c("x", "y", "ExonID", "State", "Filling_status", "AnnotStatus", "IsoformExs")} #create empty dataframe.
if ("first;last" %in% unique(as.vector(plotting_table$ExPosition))) {
  first_last_ex_df = get_first_last_ex_df(plotting_table)} else {
  first_last_ex_df = data.frame(matrix(ncol = 7, nrow = 0)); colnames(first_last_ex_df) = c("x", "y", "ExonID", "State", "Filling_status", "AnnotStatus", "IsoformExs")} #create empty dataframe.

#Only plot the exon length for the reference gene.
plotting_table$ExonLength[plotting_table$GeneID != my_gene] = NA

######### Get unique table for names
unique_table_for_names = unique(plotting_table[c("Species", "ExonNumberPlot", "GeneID", "Order")])
unique_table_for_names = unique_table_for_names[rev(order(unique_table_for_names$Order)),]


##################################################
########### Make the plot ########################
##################################################

my_plot = ggplot()  +
  geom_rect(data=internal_ex_df, aes(xmin=FakeStart, xmax=FakeStop, ymin=Order, ymax=Order+0.5, alpha=AnnotStatus, fill=Filling_status, linetype=State, size=IsoformExs), color=internal_ex_df$IsoformExs) + #internal exons.
  geom_polygon(data=first_ex_df, aes(x=x, y=y, alpha=AnnotStatus, group=ExonID, fill=Filling_status, linetype=State, size=IsoformExs), color=first_ex_df$IsoformExs) + #first exons.
  geom_polygon(data=last_ex_df, aes(x=x, y=y, alpha=AnnotStatus, group=ExonID, fill=Filling_status, linetype=State, size=IsoformExs), color=last_ex_df$IsoformExs) + #last exons.
  geom_polygon(data=first_last_ex_df, aes(x=x, y=y, alpha=AnnotStatus, group=ExonID, fill=Filling_status, linetype=State, size=IsoformExs), color=first_last_ex_df$IsoformExs) + #exons which are both first and last.
  scale_fill_manual(values=group_colors_vector, name = "Exs", labels=c("default", paste0(interesting_exons, " (", my_query_species,")"))) + #the order of the labels should be the same as in group_colors_vector. 
  #geom_rect(data=internal_ex_df, aes(xmin=FakeStart, xmax=FakeStop, ymin=Order, ymax=Order+0.5, alpha=State), fill="lightsteelblue3", color="dimgray") + #internal exons.
  #geom_polygon(data=first_ex_df, aes(x=x, y=y, alpha=State, group=ExonID), fill="lightsteelblue3", color="dimgray") + #first exons.
  #geom_polygon(data=last_ex_df, aes(x=x, y=y, alpha=State, group=ExonID), fill="lightsteelblue3", color="dimgray") + #last exons.
  #geom_polygon(data=first_last_ex_df, aes(x=x, y=y, alpha=State, group=ExonID), fill="lightsteelblue3", color="dimgray") + #exons which are both first and last.
  
  #scale_alpha_manual(values=c("Exon"=1, "Exon_added"=0.4)) + #different transparencies to exons actually matching or not the reference gene.
  scale_alpha_manual(values=c("annotated"=1, "not_annotated"=0)) + #color depending on the annotation status.
  scale_linetype_manual(values=c("Exon"="solid", "Exon_added"="dashed")) +
  scale_size_manual(values=c("brown2"=2, "black"=0.5)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank() ,
        legend.position = "bottom",
        legend.title = element_text(color="black", size=18),
        legend.text = element_text(color="black", size=15),
        plot.title = element_text(color="black", hjust=0, size=20)
  )  +
  xlim(-60,max(plotting_table$FakeStop)+5) + #limit axis +
  ggtitle(paste0(unique(subset(plotting_table, Species==my_query_species)$GeneID), " and Orthologs:\nHighlighted isoform: ", my_isorform_id)) +
  guides(alpha=FALSE, size=FALSE) +
  geom_text(data=(unique_table_for_names), aes(x=-50, y=unique(unique_table_for_names$Order)+0.25, 
                                               label=paste0(unique_table_for_names$Species, ", ",
                                                            sub(".*;", "", unique(unique_table_for_names$ExonNumberPlot)), " ex.\n",
                                                            unique_table_for_names$GeneID)), size=7) +
  geom_text(aes(x=plotting_table$FakeStart+(plotting_table$FakeStop-plotting_table$FakeStart)/2, y=plotting_table$Order+0.25, label=plotting_table$Levels), size=7) + #plot number of matching exons
  geom_text(aes(x=plotting_table$FakeStart+(plotting_table$FakeStop-plotting_table$FakeStart)/2, y=plotting_table$Order+0.75, label=plotting_table$ExonLength+1), size=7) + #plot the exon length
  geom_point(data=plotting_table, aes(x=FakeStart-1.5, y=Order+0.25, color=FinalPhaseUp), shape=8) +
  geom_point(data=plotting_table, aes(x=FakeStop+1.5, y=Order+0.25, color=FinalPhaseDown), shape=8) +
  scale_color_manual(values=c("0"="coral3","1"="forestgreen","2"="mediumblue", "shit"="shit"), name = "Intron Phases",  labels=c("0", "1", "2"), breaks=c("0", "1", "2"))

#Save pdf to output file (in a folder called exint plot)
#This is generate the right proportions in the plot.
my_width = as.numeric(nrow(subset(plotting_table, GeneID==my_gene)))+10 #number of exons 
my_height = length(unique(as.vector(plotting_table$GeneID))) #Number of orthologs
pdf(paste0(my_folder, "exint_plot.pdf"), width=my_width, height=my_height)
my_plot
dev.off()
