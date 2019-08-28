#!/usr/bin/env perl
#use warnings;
use strict;
my $i1=$ARGV[0]; #
my $out=$ARGV[1];
my ($id1, $sp1, $sp2, $l, $id2, $n);
my %prot;
my %ccds;
my %refseq;
my %genecode;
my %tr;
my @line;
my %ref;
my %longest;
my %size;
my %lsize;

#ENSP00000379387	ENSG00000185219	1031
open (INONE,"$i1") || die "Cannot open $i1\n";
while (<INONE>){
    chomp($_); 
    @line=split(/\t/,$_);
    if (!$longest{$line[1]}){
	$prot{$line[1]}=1;
	$lsize{$line[1]}=$line[2];
	$longest{$line[1]}=$line[0];
    }
    elsif ($line[2]>$lsize{$line[1]}){
	$lsize{$line[1]}=$line[2];
	$longest{$line[1]}=$line[0];
    }
}
close (INFILE);

my @keys=keys(%prot);
my $el;
my $t;
my $sc=scalar(@keys);
open (OUT,">$out");
foreach $el (@keys){
    if ($longest{$el})  {
	$t=$longest{$el};
	print OUT "$el\t$longest{$el}\n";
    }
}














