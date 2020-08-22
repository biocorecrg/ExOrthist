#!/usr/bin/env perl
use warnings;
use strict;

my $infile1=$ARGV[0]; ## Gene cluster file
my $infile2=$ARGV[1]; ## Exons scores by species pair => Best_score_exon_hits_pairs.txt
my $nc=$ARGV[2]; ## Number of clusters in each part => default 500
my $outfile=$ARGV[3]; ## Outfile name per folder

my %cid;
open (INFILE, $infile1) || die "It cannot open the input file ($infile1)\n";
# Format: GF000009	Hsa	ENSG00000257008	GPR142	OG0000061O=0,R=1,F>R
while (<INFILE>){
    chomp($_);
    my @l=split(/\t/,$_);
    $cid{$l[2]}=$l[0];
}
close INFILE;

my $tmp="tmp.out";
open (TMP, ">$tmp") || die "Tt cannot open the temporary tmp.out";
open (INFILE2, $infile2) || die "It cannot open the Exons scores by species pair ($infile2)\n";
# Format: BL00000	Sc0000095:751518-751620:-	CD14496	clodip2_s0015:466798-466932:+	Bla	Cdi
while (<INFILE2>){
    chomp($_);
    my @l=split(/\t/,$_);
    if (defined $cid{$l[0]}){  
	print TMP "$cid{$l[0]}\t$_\n";
    } 
    elsif (defined $cid{$l[2]}){
	print TMP "$cid{$l[2]}\t$_\n";
    }
}
close INFILE2;
close (TMP);

system ("cat $tmp | sort -k1 > sorted_scores.out");
my $cid="NA"; 
my ($p, $part, $cls);
my %cid2;

system "rm -r PART_*/" if (-e "PART_01");

open (TEMP_INFILE, "sorted_scores.out") || die "It cannot open the temporary input sorted_scores.out\n";
while (<TEMP_INFILE>){
    chomp($_);
    my @l=split(/\t/,$_);

    if ($cid eq "NA"){ # first case
	$cid2{$l[0]}=1;
	$cid=$l[0];
	$part="PART_01";
	$p=1;
	$cls=1;
	system ("mkdir PART_01");
	$outfile="$part/cls_$cid.tab";
	open (SC, ">>$outfile");
	print SC "$l[0]|$l[1]|$l[2]|$l[5]\t$l[0]|$l[3]|$l[4]|$l[6]\n";
    }
    elsif (!$cid2{$l[0]}){
	if ($cls < $nc){
	    $cls++;
	    $cid2{$l[0]}=1;
	    $cid=$l[0];
	    $outfile="$part/cls_$cid.tab";
	    open (SC, ">>$outfile");
	    print SC "$l[0]|$l[1]|$l[2]|$l[5]\t$l[0]|$l[3]|$l[4]|$l[6]\n";
	}
	else { # re-starts file
	    $cls=1;
	    $p++;
	    $cid2{$l[0]}=1;
	    $cid=$l[0];
	    if ($p<=9) { $part="PART_0$p";}
	    else { $part="PART_$p";}	
	    close SC;
	    system ("mkdir $part");
	    $outfile="$part/cls_$cid.tab";
	    open (SC, ">>$outfile");
	    print SC "$l[0]|$l[1]|$l[2]|$l[5]\t$l[0]|$l[3]|$l[4]|$l[6]\n";
	}
    }
    elsif ($cid2{$l[0]}){
	print SC "$l[0]|$l[1]|$l[2]|$l[5]\t$l[0]|$l[3]|$l[4]|$l[6]\n";
    }    
}
close TEMP_INFILE;
close SC;
print STDERR "Number of parts:\t$p\n";

system "rm tmp.out sorted_scores.out";
