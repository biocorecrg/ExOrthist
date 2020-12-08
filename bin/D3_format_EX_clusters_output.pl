#!/usr/bin/env perl
use warnings;
use strict;

##Joining the parts
open (T1, ">EX_clusters_info.tab");
print T1 "ExonID\tClusterID\tOut_degree\tIn_degree\tSPECIES_exons_in_cluster\tTOT_exons_in_cluster\tN_reciprocals\tMembership_score\n";
close T1;
system "cat EXcls_* >> EX_clusters_info.tab";
system "cat unclustered_EXs_* >> unclustered_EXs.txt"; #Fede: joining the unclustered exons (exons with hits but ending up in single-exon clusters)

open (TMP, ">tmp_excls.tab");
open (INFILE, "EX_clusters_info.tab");
# Format: GF_000006|ENSG00000103266|chr16:680526-680684:+|Hs2 1 2 2 1 3 2 1
while (<INFILE>){
    chomp($_);
    if ($_=~/\|/){
	$_=~s/\"//g; # no harm
	$_=~s/\|/\t/g;
	my @l=split(/\t/,$_);
	my $id=$l[0].".".sprintf("%.3d",$l[4]);
	$l[10] = sprintf("%.3f",$l[10]);
	print TMP  "$id\t$l[1]\t$l[2]\t$l[3]\t$l[10]\n"; # adds Memb score
    }   
}
close (TMP);

open (O, ">EX_clusters.tab") || die "Cannot open the output file\n";
print O "ExCluster_ID\tGeneID\tCoordinate\tSpecies\tMembership_score\n";
close O;

system ("cat tmp_excls.tab | sort -k1 >> EX_clusters.tab");
system "gzip EX_clusters_info.tab";
system "rm tmp_excls.tab";
