#!/usr/bin/env perl
use warnings;
use strict;

##SCRIPT FOR SCORING EXONS BTW PAIR OF SPECIES##
##
my $i1=$ARGV[0]; ## score all exons: Score_all_exons.txt
my $out=$ARGV[1]; ## Best_scores_pair_exons.txt
my (%print);
my ($pair1, $pair2);
my @l;
open (OUT, ">$out");
#GF000004.03ENSP00000264360|ENSG00000138650exon_11-878chr4133150141-133152771+1ENSMUSP00000131073|ENSMUSG000000491001-87897.6198.1800.00ENSMUSP00000131073|ENSMUSG00000049100exon_11-878chr345379253-45381883+Hs2Mm2
open (IN,"$i1") || die "Missing exon file ($i1)\n";  
while (<IN>){
    chomp($_);
    @l=split(/\t/,$_);
    $pair1=$_;
    $pair2=$l[0]."\t".$l[14]."\t".$l[15]."\t".$l[16]."\t".$l[17]."\t".$l[18]."\t".$l[19]."\t".$l[7]."\t".$l[1]."\t".$l[3]."\t".$l[10]."\t".$l[11]."\t".$l[12]."\t".$l[13]."\t".$l[1]."\t".$l[2]."\t".$l[3]."\t".$l[4]."\t".$l[5]."\t".$l[6]."\t".$l[21]."\t".$l[20];
    if (!$print{$pair1}){
	$print{$pair1}=1;
	print OUT "$pair1\n";
    }
    if (!$print{$pair2}){
	$print{$pair2}=1;
	print OUT "$pair2\n";
    }   
}
close OUT;
close IN;
