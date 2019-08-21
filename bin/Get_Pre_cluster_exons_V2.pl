#!/usr/bin/perl
#use warnings;
#use strict;

my $i1=$ARGV[0]; ##Gene cluster file
my $i2=$ARGV[1]; ##Exons scores by species pair
my $nc=$ARGV[2]; ##number of clusters in each part
my $out=$ARGV[3]; ##Outfile name
$fs=`ls cls_*`;
$dir=`pwd`;
chomp($dir);
@files=split(/\n/,$fs);
$f=0;
$num=scalar(@files);
for ($l=0; $l<scalar(@files); $l++){
	$in=$files[$l];
	$out=$dir."/EX$in";
	$in=$dir."/".$files[$l];
	$f++;
	$rf="cls$f.R";
	open (RSC, ">$rf");
	print RSC "#!/usr/bin/env Rscript\n\n";
	print RSC "library(igraph)\n";
	print RSC "e <- read.table (\"$in\")\n";
	print RSC "m <- as.matrix (e)\n";
	print RSC "g2 <- graph_from_edgelist(m, directed = FALSE)\n";
	print RSC "memb <- cluster_edge_betweenness(g2, weights =NULL)\n";
	print RSC "e <- membership (memb)\n";
	print RSC "write.table(as.matrix(e), file =\"$out\", sep=\"\\t\")\n";
	close (RSC);
}
$inf=$dir."/cls\${SGE_TASK_ID}.R";
`mkdir tmp`;
$jout="$dir/tmp";
open (ARR, ">array_job.sh");
print ARR "#!/bin/bash\n";
print ARR "#\$ -N array_test\n";
print ARR "#\$ -t 1-$num\n";
print ARR "#\$ -l virtual_free=5G,h_rt=6:00:00\n";
print ARR "#\$ -o $jout/\$TASK_ID.out\n";
print ARR "#\$ -e $jout/\$TASK_ID.err\n";
print ARR "\nRscript --vanilla $inf\n";
