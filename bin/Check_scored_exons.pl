#!/usr/bin/env perl
use warnings;
#use strict;


my $dir=$ARGV[0];
my $gcl=$ARGV[1];
my $sp1=$ARGV[2];
my $sp2=$ARGV[3];
my $f1=$dir."/$gcl-part_*";
my $l1=`ls $f1`;
my $f2=$dir."/score_exons_*_part_*.txt";
my $l2=`ls $f2`;
my $f3=$dir."/score_introns_*part_*.txt";
my $l3=`ls $f3`;
my $f4=$dir."/score_proteins_*part_*.txt";
my $l4=`ls $f4`;

my @s1=split(/\n/,$l1);
my @s2=split(/\n/,$l2);
my @s3=split(/\n/,$l3);
my @s4=split(/\n/,$l4);
$b=0;
$o1=$dir."/missing_scores_".$sp1."_".$sp2.".log";
open (OUT, ">$o1");
for ($w=0; $w<scalar(@s1); $w++){
	$p=$w+1;
	$i1=$s1[$w];
	$i2=$s2[$w];
	$i3=$s3[$w];
	$i4=$s4[$w];	
	my (%gn, %miss);
	print "$p\n";
	open (FILE, "$i4");
	#GF000020	Protein	ENSMODP00000009305|ENSMODG00000007499	ENSP00000462903|ENSG00000221926	46.57	72.13	42.49	0.35	Mdo	Hs2
	while (<FILE>){
		chomp($_);
		@l=split(/\t/,$_);
		if ($l[1] ne "CID"){
			if ($l[2]=~/\|/){
				@k1=split(/\|/,$l[2]);
				$gn{$k1[1]}=1;
			}
			if ($l[3]=~/\|/){
				@k2=split(/\|/,$l[3]);				
				$gn{$k2[1]}=1;
			}
		}
	}
	open (FILE, "$i2");
	#GF000001	ENSMUSP00000122881|ENSMUSG00000041439	exon_3	548-601	chr1	52663285-52663446	-	1	ENSP00000281416|ENSG00000151690	545-598	100.00
	while (<FILE>){
		chomp($_);
		@l=split(/\t/,$_);
		if ($l[1] ne "CID"){
			if ($l[1]=~/\|/){
				@k1=split(/\|/,$l[1]);
				$gn{$k1[1]}=1;
			}
			if ($l[8]=~/\|/){
				@k2=split(/\|/,$l[8]);				
				$gn{$k2[1]}=1;
			}
		}
	}
	open (FILE, "$i2");
	#GF000001	ENSMUSP00000122881|ENSMUSG00000041439	exon_3	548-601	chr1	52663285-52663446	-	1	ENSP00000281416|ENSG00000151690	545-598	100.00
	while (<FILE>){
		chomp($_);
		@l=split(/\t/,$_);
		if ($l[1] ne "CID"){
			if ($l[1]=~/\|/){
				@k1=split(/\|/,$l[1]);
				$gn{$k1[1]}=1;
			}
			if ($l[8]=~/\|/){
				@k2=split(/\|/,$l[8]);				
				$gn{$k2[1]}=1;
			}
		}
	}
	open (FILE, "$i3");
	#GF000001	ENSP00000406837|ENSG00000151690	intron_1	6	2	ENSMUSP00000084991|ENSMUSG00000041439	intron_2	547	2	0	u,u	10	Hs2	Mm2
	while (<FILE>){
		chomp($_);
		@l=split(/\t/,$_);
		if ($l[1] ne "CID"){
			if ($l[1]=~/\|/){
				@k1=split(/\|/,$l[1]);
				$gn{$k1[1]}=1;
			}
			if ($l[5]=~/\|/){
				@k2=split(/\|/,$l[5]);				
				$gn{$k2[1]}=1;
			}
		}
	}
	open (FILE, "$i1");
	#GF000001	Hs2	ENSG00000151690	MFSD6	OG0000376	O=0,R=1,F>R
	while (<FILE>){
		chomp($_);
		@l=split(/\t/,$_);
		if ($l[1] eq $sp1 || $l[1] eq $sp2){
			if (!$gn{$l[2]}){
				print OUT "$_\tpart\t$p\n";
				$b=1;
			}
		}
	}
}
if ($b==1) { print STDERR "Missing genes to be scored... Check file $o1 !!!\n"; }


