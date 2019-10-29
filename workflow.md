step 1
generate_annotations.pl -sp sp1 -GTF sp1.gtf -G sp1.fasta 
generate_annotations.pl -sp sp2 -GTF sp2.gtf -G sp2.fasta

step 2
splitcluster.pl cluster_file.tab // cluster of genes) 
- sp1_sp2_01.tab
- sp1_sp2_02.tab
- sp1_sp2_03.tab
- sp1_sp3.tab
...

step 3
align.pl sp1.gtf sp1.fasta sp2.gtf sp2.fasta sp1_sp2_01.tab
align.pl sp1.gtf sp1.fasta sp2.gtf sp2.fasta sp1_sp2_02.tab
align.pl sp1.gtf sp1.fasta sp2.gtf sp2.fasta sp1_sp2_03.tab

step 4
refinement on splitted alignment

step 5
collecting evertything
output cluster of exons.

