#!/usr/bin/env perl
use warnings;
use strict;

##SCRIPT TO GET RIGD OF REDUNDANT PAIR OF EXONS##
my $infile1=$ARGV[0]; ## overlapping exon info (Overlap_exons_by_sp.tab)
my $infile2=$ARGV[1]; ## best score hit exons (Best_score_exon_hits_filtered_$m-$int-$id.tab)
my $outfile=$ARGV[2]; ## outfile name
my $infile3=$ARGV[3]; ## liftover pairs, if provided

my %exid; my %max; my %bh;
my %pair;

open (IN, $infile1) || die "It cannot open overlapping exons info ($infile1)\n"; 
# Format: OV_EX_Bla_36	BL00001	Sc0000095:474219-474451:-	Bla	1015
while (<IN>){
    chomp($_);
    my @l=split(/\t/,$_);
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
close IN;

open (OUT, ">$outfile") || die "It cannot open the output file $outfile\n";
open (IN2, $infile2) || die "It cannot open fil with best score exons ($infile2)\n"; 
## Format: BL00000	Sc0000095:736652-736717:-	NEMVEDRAFT_v1g239017	NEMVEscaffold_11:253907-253972:-	Bla	Nve SCORE?
<IN2>; # discards header
while (<IN2>){
    chomp($_);
    my @l=split(/\t/,$_);
    my $id1=$l[0]."\t".$l[1];
    my $ovid1=$exid{$id1};
    my $bh1=$bh{$ovid1};
    my $id2=$l[2]."\t".$l[3];
    my $ovid2=$exid{$id2};
    my $bh2=$bh{$ovid2};
    
    # it originally excluded no-symetric hits
    if (!$pair{$bh1."\t".$bh2}){# && !$pair{$bh2."\t".$bh1}){ # now it allows reciprocal matches   
	print OUT "$bh1\t$bh2\t$l[4]\t$l[5]\n";
	$pair{$bh1."\t".$bh2}=1; # simply to avoid repeated lines
    }
}
close IN2;

if ($infile3){
    open (IN3,"$infile3") || die "It cannot open the provided liftover file ($infile3)"; 
    ## Format: BL00000	Sc0000095:736652-736717:-	NEMVEDRAFT_v1g239017	NEMVEscaffold_11:253907-253972:-	Bla	Nve
    <IN3>;
    while (<IN3>){
	chomp($_);
	my @l=split(/\t/,$_);
	my $id1=$l[0]."\t".$l[1];
	my $ovid1=$exid{$id1};
	my $bh1=$bh{$ovid1};
	my $id2=$l[2]."\t".$l[3];
	my $ovid2=$exid{$id2};
	my $bh2=$bh{$ovid2};
	
	if (!$pair{$bh1."\t".$bh2}){# && !$pair{$bh2."\t".$bh1}){ # allows reciprocal matches
	    print OUT "$bh1\t$bh2\t$l[4]\t$l[5]\n";
	    $pair{$bh1."\t".$bh2}=1; # to avoid repeated lies
	}
    }
    close IN3;
}
