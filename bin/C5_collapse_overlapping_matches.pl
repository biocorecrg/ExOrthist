#!/usr/bin/env perl
use warnings;
use strict;

##SCRIPT TO GET RIGD OF REDUNDANT PAIR OF EXONS##
my $infile1=$ARGV[0]; ## overlapping exon info (overlapping_EXs_by_species.tab)
my $infile2=$ARGV[1]; ## best score hit exons (filtered_best_scored_EX_matches_by_targetgene.tab);
my $outfile=$ARGV[2]; ## outfile name => filtered_best_scored_EX_matches_by_targetgene-NoOverlap.tab
my $infile3=$ARGV[3]; ## liftover pairs, if provided

my %exid; my %max; my %bh;
my %pair;

open (IN, $infile1) || die "It cannot open overlapping exons info ($infile1)\n"; 
# Format: OV_EX_Bla_36	BL00001	Sc0000095:474219-474451:-	Bla	1015
while (<IN>){
    chomp($_);
    my @l=split(/\t/,$_);
    $exid{$l[1]."\t".$l[2]}=$l[0]; # individual_exon => Exon_group_ID
    if (!$max{$l[0]}){
	$max{$l[0]}=$l[4];
	$bh{$l[0]}=$l[1]."\t".$l[2]; # exon_group_ID => individual_exon
    }
    elsif ($l[4]>$max{$l[0]}) {
	$max{$l[0]}=$l[4];
	$bh{$l[0]}=$l[1]."\t".$l[2];
    }
}
close IN;

open (OUT, ">$outfile") || die "It cannot open the output file $outfile\n";

### First is to pre-load the data, to know which exons have a match with given genes
my %exon_to_gene_match; # keeps info of whether an actual exon has a hit in the target gene
# liftover
if ($infile3){
    open (IN3_A,"$infile3") || die "It cannot open the provided liftover file ($infile3)"; 
    ## Format: BL00000	Sc0000095:736652-736717:-	NEMVEDRAFT_v1g239017	NEMVEscaffold_11:253907-253972:-	Bla	Nve
    <IN3_A>;
    while (<IN3_A>){
	chomp($_);
	my @l=split(/\t/,$_);
	my $id1=$l[0]."\t".$l[1];
	my $g2=$l[2];
	$exon_to_gene_match{$id1}{$g2}=1;
    }
    close IN3_A;
}
# Exortist
open (IN2_A, $infile2) || die "It cannot open fil with best score exons ($infile2)\n"; 
## Format: BL00000	Sc0000095:736652-736717:-	NEMVEDRAFT_v1g239017	NEMVEscaffold_11:253907-253972:-	Bla	Nve
<IN2_A>; # discards header
while (<IN2_A>){
    chomp($_);
    my @l=split(/\t/,$_);
    my $id1=$l[0]."\t".$l[1];
    my $g2 = $l[2];
    $exon_to_gene_match{$id1}{$g2}=1;
}
close IN2_A;


### Then, it uses liftOver info, if available
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
	my $g2=$l[2];
	
	if ($id1 eq $bh1){ # checks if it's the best rep of the exon group
	    if (!$pair{$ovid1."\t".$g2}){ # to avoid redundancy
		print OUT "$bh1\t$bh2\t$l[4]\t$l[5]\n";
		$pair{$ovid1."\t".$g2}=1; # to avoid repeated lines
	    }
	}
	elsif (!$exon_to_gene_match{$bh1}{$g2}){ # it's not the best rep, but the best rep does not hit this gene
	    if (!$pair{$ovid1."\t".$g2}){ # to avoid redundancy; only one non-rep per exon group
		print OUT "$bh1\t$bh2\t$l[4]\t$l[5]\n";
		$pair{$ovid1."\t".$g2}=1; # to avoid repeated lines
	    }	    
	}
    }
    close IN3;
}

open (IN2, $infile2) || die "It cannot open fil with best score exons ($infile2)\n"; 
## Format: BL00000	Sc0000095:736652-736717:-	NEMVEDRAFT_v1g239017	NEMVEscaffold_11:253907-253972:-	Bla	Nve
<IN2>; # discards header
while (<IN2>){
    chomp($_);
    my @l=split(/\t/,$_);
    my $id1=$l[0]."\t".$l[1];
    my $ovid1=$exid{$id1}; # Gets exon_group_ID
    my $bh1=$bh{$ovid1};  # best representative of the exon_group_ID
    my $id2=$l[2]."\t".$l[3];
    my $ovid2=$exid{$id2};
    my $bh2=$bh{$ovid2};
    my $g2 = $l[2];
    
    # it originally excluded no-symetric hits
    if ($id1 eq $bh1){ # checks if it's the best rep of the exon group
	if (!$pair{$ovid1."\t".$g2}){ # to avoid redundancy
	    print OUT "$bh1\t$bh2\t$l[4]\t$l[5]\n";
	    $pair{$ovid1."\t".$g2}=1; # to avoid repeated lines
	}
    }
    elsif (!$exon_to_gene_match{$bh1}{$g2}){ # it's not the best rep, but the best rep does not hit this gene
	if (!$pair{$ovid1."\t".$g2}){ # to avoid redundancy; only one non-rep per exon group
	    print OUT "$bh1\t$bh2\t$l[4]\t$l[5]\n";
	    $pair{$ovid1."\t".$g2}=1; # to avoid repeated lines
	}	    
    }
}
close IN2;

