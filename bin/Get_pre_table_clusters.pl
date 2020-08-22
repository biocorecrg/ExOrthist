#!/usr/bin/env perl
use warnings;
use strict;

##Joining the parts
`cat EXcls_* > pre_cluster_exons.txt`; 

open (TMP, ">tmp_excls.tab");
open (INFILE,"pre_cluster_exons.txt");
# Format: "GF000033|WHL22.236185|Scaffold196:625640-625797:+|Spu"11
while (<INFILE>){
    chomp($_);
    if ($_=~/\|/){
	$_=~s/\"//g;
	$_=~s/\|/\t/g;
	my @l=split(/\t/,$_);
	my $id=$l[0].".".sprintf("%.3d",$l[4]);
	print TMP  "$id\t$l[1]\t$l[2]\t$l[3]\n";
    }   
}
close (TMP);

system ("cat tmp_excls.tab | sort -k1 > Exon_Clusters.tab");
system "rm pre_cluster_exons.txt tmp_excls.tab";
