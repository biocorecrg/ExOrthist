#!/usr/bin/env perl
use warnings;
#use strict;

##SCRIPT FOR GET RIG OF REDUNDANT PAIR OF EXONS##
##
my $i1=$ARGV[0]; ##overlapping exons
my $i2=$ARGV[1]; ##best score hit exons
my $out=$ARGV[2]; ##liftover exons
my $i3=$ARGV[3]; ##outfile name

# OV_EX_Bla_36	BL00001	Sc0000095:474219-474451:-	Bla	1015
    open (IN,"$i1") || die "Missing overlapping exons"; 
while (<IN>){
    chomp($_);
    @l=split(/\t/,$_);
    $exid{$l[1]."\t".$l[2]}=$l[0];
    if (!$max{$l[0]}){
	$max{$l[0]}=$l[4];
	$bh{$l[0]}=$l[1]."\t".$l[2];
    }
    elsif ($l[4]>$max{$l[0]}) {
	$max{$l[0]}=$l[4];
	$bh{$l[0]}=$l[1]."\t".$l[2];
    }
}

open (OUT, ">$out");
##BL00000	Sc0000095:736652-736717:-	NEMVEDRAFT_v1g239017	NEMVEscaffold_11:253907-253972:-	Bla	Nve
open (IN,"$i2") || die "Missing best score exons"; 
while (<IN>){
    chomp($_);
    @l=split(/\t/,$_);
    $id1=$l[0]."\t".$l[1];
    $ovid1=$exid{$id1};
    $bh1=$bh{$ovid1};
    $id2=$l[2]."\t".$l[3];
    $ovid2=$exid{$id2};
    $bh2=$bh{$ovid2};
    
    if (!$pair{$bh1."\t".$bh2} && !$pair{$bh2."\t".$bh1}){
	print OUT "$bh1\t$bh2\t$l[4]\t$l[5]\n";
	$pair{$bh1."\t".$bh2}=1;
    }
}
if ($i3){
open (IN,"$i3") || die "Missing lift over file"; 
while (<IN>){
    chomp($_);
    @l=split(/\t/,$_);
    $id1=$l[0]."\t".$l[1];
    $ovid1=$exid{$id1};
    $bh1=$bh{$ovid1};
    $id2=$l[2]."\t".$l[3];
    $ovid2=$exid{$id2};
    $bh2=$bh{$ovid2};
    
    if (!$pair{$bh1."\t".$bh2} && !$pair{$bh2."\t".$bh1}){
	print OUT "$bh1\t$bh2\t$l[4]\t$l[5]\n";
	$pair{$bh1."\t".$bh2}=1;
    }
}
}







