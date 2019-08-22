#!/usr/bin/perl
#use warnings;
#use strict;

##
my $i1=$ARGV[0]; ##Gene cluster file
my $i2=$ARGV[1]; ##Exons scores by species pair
my $nc=$ARGV[2]; ##number of clusters in each part
my $out=$ARGV[3]; ##Outfile name

#GF000009	Hsa	ENSG00000257008	GPR142	OG0000061O=0,R=1,F>R
my ($id1, $sp1, $sp2, $l, $id2, $n);
open (INFILE,"$i1");
while (<INFILE>){
    chomp($_);
    @l=split(/\t/,$_);
    $cid{$l[2]}=$l[0];
}
#BL00000	Sc0000095:751518-751620:-	CD14496	clodip2_s0015:466798-466932:+	Bla	Cdi
my $tmp="tmp.out";
open (TMP, ">$tmp");
open (INFILE,"$i2");
while (<INFILE>){
    chomp($_);
    @l=split(/\t/,$_);
    if ($cid{$l[0]}){  
	print TMP "$cid{$l[0]}\t$_\n";
	
    }elsif ($cid{$l[2]}){
	print TMP "$cid{$l[2]}\t$_\n";
    }
}
close (TMP);
system ("cat $tmp | sort -k1 > sorted_scores.out");
$cid="NA";
open (INFILE,"sorted_scores.out");
while (<INFILE>){
    chomp($_);
    @l=split(/\t/,$_);
    if ($cid eq "NA"){
	$cid2{$l[0]}=1;
	$cid=$l[0];
	$part="PART_01";
	$p=1;
	$cls=1;
	system ("mkdir PART_01");
	$out="$part/cls_$cid.tab";
	open (SC, ">>$out");
	print SC "$l[0]|$l[1]|$l[2]|$l[5]\t$l[0]|$l[3]|$l[4]|$l[6]\n";
    }
    elsif (!$cid2{$l[0]}){
	if ($cls<$nc){
		$cls++;
	   	$cid2{$l[0]}=1;
		$cid=$l[0];
		$out="$part/cls_$cid.tab";
		open (SC, ">>$out");
		print SC "$l[0]|$l[1]|$l[2]|$l[5]\t$l[0]|$l[3]|$l[4]|$l[6]\n";
	}
	else {
	   	$cls=1;
		$p++;
	   	$cid2{$l[0]}=1;
		$cid=$l[0];
		if ($p<=9) { $part="PART_0$p"; }
		else { $part="PART_$p";    }	
		system ("mkdir $part");
		$out="$part/cls_$cid.tab";
		open (SC, ">>$out");
		print SC "$l[0]|$l[1]|$l[2]|$l[5]\t$l[0]|$l[3]|$l[4]|$l[6]\n";
	}

    }
    elsif ($cid2{$l[0]}){
	print SC "$l[0]|$l[1]|$l[2]|$l[5]\t$l[0]|$l[3]|$l[4]|$l[6]\n";
    }
    
}
print STDERR "Number of parts:\t$p\n\n";

