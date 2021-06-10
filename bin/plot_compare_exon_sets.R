### R script to generate a graphical output for the compare_exon_sets.pl ExOrsthist module
#author: Federica Mantica (federica.mantica93@gmail.com)

#The script takes four arguments:
#1: stat input
#2: query species
#3: target species
#4: path to output file

#stat input is structured as follows (probably temporary structure)
#NB: first line = header
# my_class     category          freq  percentage  conservation
# gene_status  gene_orth         453   0.56        conserved
# gene_status  no_gene_orth      374   0.45        not_conserved
# exon_status  exon_orth         43    0.09        conserved
# exon_status  no_exon_orth      410   0.91        not_conserved
# reg_status   reg_exon_orth     4     0.09        conserved
# reg_status   reg_exon_no_orth  10    0.23        conserved
# reg_status   no_reg_exon       29    0.68        not_conserved
# all_regs     orth              4     0.05        conserved
# all_regs     best_hit          8     0.10        conserved
# all_regs     unclear           7     0.09        not_conserved
# all_regs     no_orth           61    0.76        not_conserved

##### upload libraries #####
suppressMessages(library(ggplot2))

##### read arguments #####
args <- commandArgs(TRUE)
input_file = args[1] #stat input
sp1 = args[2] #species query ID
sp2 = args[3] #species target ID
output_file = args[4] #Path to output file

##### generate variables #####
class_levels = c("gene_status", "exon_status", "reg_status", "all_regs")
category_levels = c("gene_orth", "no_gene_orth", 
                    "exon_orth", "no_exon_orth", 
                    "reg_exon_orth", "reg_exon_no_orth", "no_reg_exon", 
                    "orth", "best_hit", "unclear", "no_orth")

label_levels = c("gene_orth_conserved", "no_gene_orth_not_conserved", 
                 "exon_orth_conserved", "no_exon_orth_not_conserved",
                 "reg_exon_orth_conserved", "reg_exon_no_orth_conserved", "no_reg_exon_not_conserved",
                 "orth_conserved", "best_hit_conserved", "unclear_not_conserved", "no_orth_not_conserved")



fill_color_vector = c("gene_orth_conserved"="mediumblue",
                      "no_gene_orth_not_conserved"="mediumblue",
                      "exon_orth_conserved"="goldenrod2",
                      "no_exon_orth_not_conserved"="goldenrod2",
                      "reg_exon_orth_conserved"="firebrick1",
                      "reg_exon_no_orth_conserved"="springgreen4",
                      "no_reg_exon_not_conserved"="springgreen4",
                      "orth_conserved"="darkorchid4",
                      "best_hit_conserved"="darkorchid4",
                      "unclear_not_conserved"="darkorchid3",
                      "no_orth_not_conserved"="darkorchid3")

alpha_vector = c("gene_orth_conserved"=0.8,
                 "no_gene_orth_not_conserved"=0.3,
                 "exon_orth_conserved"=0.8,
                 "no_exon_orth_not_conserved"=0.3,
                 "reg_exon_orth_conserved"=0.8,
                 "reg_exon_no_orth_conserved"=0.8,
                 "no_reg_exon_not_conserved"=0.3,
                 "orth_conserved"=0.8, 
                 "best_hit_conserved"=0.6,
                 "unclear_not_conserved"=0.5,
                "no_orth_not_conserved"=0.2)
 
legend_labels = c("gene_orth_conserved"= paste0("conserved gene in ", sp2),
                  "no_gene_orth_not_conserved"= paste0("not conserved gene in ", sp2),
                  "exon_orth_conserved"= paste0("ortholog in ", sp2),
                  "no_exon_orth_not_conserved"= paste0("no ortholog in ", sp2),
                  "reg_exon_orth_conserved" = paste0("regulated ortholog in ", sp2),
                  "reg_exon_no_orth_conserved" = paste0("regulated not-ortholog in ", sp2),
                  "no_reg_exon_not_conserved"= paste0("no regulated exon in ", sp2),
                  "orth_conserved" = paste0("regulated orthologous exon in ", sp2),
                  "best_hit_conserved" = paste0("regulated best-hit exon in ", sp2),
                  "unclear_not_conserved" = paste0("regulated exon (unclear orthology) in ", sp2),
                  "no_orth_not_conserved" = paste0("regulated not-orthologous exon in ", sp2))


legend_title = paste0(sp1, " regulated exons with: ")


##### generate plot input #####
my_input = read.delim(input_file, header=TRUE)
my_input$my_class = factor(my_input$my_class, levels=class_levels)
my_input$category = factor(my_input$category, levels=category_levels)

#add label to join the legend
my_input$label = paste0(my_input$category, "_", my_input$conservation)
my_input$label = factor(my_input$label, levels=label_levels)

#coordinates for the first shaded area
A_first_x_coord = 1.25
A_second_x_coord = 1.75
A_first_y_coord = subset(my_input, my_class=="gene_status" & category=="no_gene_orth")$percentage
#coordinates for the second shaded area
B_first_x_coord = 2.25
B_second_x_coord = 2.75
B_first_y_coord = subset(my_input, my_class=="exon_status" & category=="no_exon_orth")$percentage
#coordinates for the third shaded area
C_first_x_coord = 3.25
C_second_x_coord = 3.75
C_first_y_coord = 1 - subset(my_input, my_class=="reg_status" & category=="reg_exon_orth")$percentage
C_second_y_coord = 1 - subset(my_input, my_class=="all_regs" & category=="orth")$percentage

#coordinates for the fourth shaded area
#compute the total number of all_regs
tot_regs = sum(subset(my_input, my_class=="all_regs")$freq)
percentage_tot_subtract = subset(my_input, my_class=="reg_status" & category=="reg_exon_no_orth")$freq/tot_regs

D_first_x_coord = 3.25
D_second_x_coord = 3.75
D_first_y_coord_1 = subset(my_input, my_class=="reg_status" & category=="no_reg_exon")$percentage
D_first_y_coord_2 = 1 - subset(my_input, my_class=="reg_status" & category=="reg_exon_orth")$percentage
D_second_y_coord_1 = subset(my_input, my_class=="all_regs" & category=="no_orth")$percentage
D_second_y_coord_2 = subset(my_input, my_class=="all_regs" & category=="no_orth")$percentage - percentage_tot_subtract

#build dataframe to draw the polygon
ids = c(rep("A",4), rep("B",4), rep("C",4), rep("D",4), rep("A",4), rep("B",4), rep("C",4), rep("D",4))
polygon_dataframe = data.frame(id = ids,
                               x = c(A_first_x_coord, A_first_x_coord, A_second_x_coord, A_second_x_coord,
                                     B_first_x_coord, B_first_x_coord, B_second_x_coord, B_second_x_coord,
                                     C_first_x_coord, C_first_x_coord, C_second_x_coord, C_second_x_coord,
                                     D_first_x_coord, D_first_x_coord, D_second_x_coord, D_second_x_coord),
                               y = c(A_first_y_coord, 1, 1, 0,
                                     B_first_y_coord, 1, 1, 0,
                                     C_first_y_coord, 1, 1, C_second_y_coord,
                                     D_first_y_coord_1, D_first_y_coord_2, D_second_y_coord_1, D_second_y_coord_2))

##### plot ######
compare_plot = ggplot(data=my_input, aes(fill=label)) +
  geom_bar(data=my_input, aes(x=my_class, y=freq, fill=label, alpha=label),
           stat="identity", position="fill", width = 0.5, color="black") +
  geom_polygon(data=polygon_dataframe, aes(x=x, y=y, group=id), 
               fill="dimgray", alpha=0.3) +
  geom_text(data=my_input, aes(x=my_class, y=freq, label=freq), 
            position=position_fill(vjust=0.5), size=3.5, color="black") +
  scale_fill_manual(values=fill_color_vector, labels=legend_labels) +
  scale_alpha_manual(values=alpha_vector, labels=legend_labels) +
  theme_minimal() +
  theme(axis.title = element_text(color="black", size=12),
        axis.text = element_text(color="black", size=10),
        axis.ticks.y = element_line(),
        plot.title = element_text(color="black", size=14, hjust=0.5),
        plot.subtitle = element_text(color="black", size=10, hjust=0.5),
        legend.key.size = unit(5,"mm")) +
  xlab("Conservation level") +
  ylab("Percentage") +
  labs(fill=legend_title, alpha=legend_title) +
  ggtitle(paste0(sp1, " vs ", sp2), "compare_exons_sets.pl stats") +
  scale_x_discrete(labels=c("gene_status"="Gene\n(genomic)", "exon_status"="Exon\n(genomic)", 
                            "reg_status"="Exon\n(genomic and\nregulatory)", "all_regs"="Exon\n(regulatory)")) +
  guides(color=FALSE)


##### save to output ######
pdf(output_file, height=5, width = 8)
suppressWarnings(print(compare_plot))
#dev.off()
