#!/usr/bin/perl
#use warnings;
#use strict;

my $i1=$ARGV[0]; ##reclustered genes
my $i2=$ARGV[1]; ##exon clusters from Module II
my $out=$ARGV[2]; ##name output file
my ($id1, $as, $ts, $l, $id2, $n);

#GF000001.01	Mdo	ENSMODG00000027707	NA
open (INONE, $i1);
while (<INONE>){
    	chomp($_); 
    	@line=split(/\t/,$_);
    	$cid{$line[2]}=$line[0];
	$s=$line[1];
	$sp{$s}=1;
}

#GF000001.004	Xtr	ENSXETG00000032975	no_gene_name	GL173719:90495-90617	123	no_vast_id
open (INONE, $i2);
while (<INONE>){
    	chomp($_); 
   	@line=split(/\t/,$_);
 	if ($_=~/GF/ && $sp{$line[3]}){ 
		@tmp=split(/\./,$line[0]);
    		if ($cid{$line[1]}){
			$nid=$cid{$line[1]}.".".$tmp[1];
			if (!$spcl{$nid}{$line[3]}){
			    $spcl{$nid}{$line[3]}=1;
			    $sps{$nid}++;				
			}
    		}
	}
}

open (OUT, ">$out");
$c=0;
open (INONE, $i2);
while (<INONE>){
    	chomp($_); 
   	@line=split(/\t/,$_);
	if ($_=~/GF/ && $sp{$line[3]}){
		@tmp=split(/\./,$line[0]);
    		if ($cid{$line[1]}){ 
			$nid=$cid{$line[1]}.".".$tmp[1];
			if ($sps{$nid}>1){
				print OUT "$nid\t$line[1]\t$line[2]\t$line[3]\n"; 			
			}

    		}
	}
}
close (OUT);

##sorting##
$tmp="tmp_".$out;

`cat $out | sort -k1 > $tmp`;
`mv $tmp $out`;
