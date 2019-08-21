#!/usr/bin/perl
#use warnings;
#use strict;

my $f=$ARGV[0]; ##path to the folder with the score file parts


##Joining the parts
`cat $f/PART_*/EXcls_* > pre_cluster_exons.txt`; 

#"GF000033|WHL22.236185|Scaffold196:625640-625797:+|Spu"11
my ($id1, $sp1, $sp2, $l, $id2, $n);
open (TMP, ">tmp_excls.tab");
open (INFILE,"pre_cluster_exons.txt");
while (<INFILE>){
    chomp($_);
    if ($_=~/\|/){
	$_=~s/\"//g;
	$_=~s/\|/\t/g;
	@l=split(/\t/,$_);
	$id=$l[0].".".sprintf("%.3d",$l[4]);
	print TMP  "$id\t$l[1]\t$l[2]\t$l[3]\n";
    }   
}
close (TMP);
system ("cat tmp_excls.tab | sort -k1 > Exon_Clusters.tab");

