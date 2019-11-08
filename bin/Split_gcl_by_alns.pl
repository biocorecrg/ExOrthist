#!/usr/bin/perl
#use warnings;
#use strict;
use Getopt::Long;
use Data::Dumper;

$N=10000;
$w=0;

### Get options:
Getopt::Long::Configure("no_auto_abbrev");
GetOptions( "gene_cluster=s" => \$cl, ##gene clusters for sp1 and sp2
	         "s1=s" => \$sp1,
	         "s2=s" => \$sp2,
	         "clean" => \$clean,
	         "split" => \$split,
	         "ex1=s" => \$ex1,
	         "ex2=s" => \$ex2,
	         "bin=s" => \$bin,
	         "N_split=i" => \$N,
	         "help" => \$help
    );

### Help
if (!defined ($cl) || !defined($sp1) || !defined($sp2) || !defined($ex1) || !defined($ex2) || defined ($help)){
    die "
Usage: Split_gcl_by_alns.pl -s1 sp1 -s2 sp2 --gene_cluster FILE --ex1 exint_file_sp1 --ex2 exint_file_sp2 [options]

OPTIONS
    --clean              Force removal of previous data [def OFF]
    --split              Split in parts the original cluster file  [def OFF]
    --N_split int        Number of alignments in subfile [def 10000] 

";
}

#>GB42790-PA|GB42790  33.1 67.2 135.2
open (ONE, "$ex1");
while (<ONE>){
    chomp($_);
    if ($_=~/\>/){
	@l=split(/\s+/,$_);
	@t=split(/\|/,$l[0]);
	$trs{$t[1]}++;
    }
}
open (TW0, "$ex2");
while (<TW0>){
    chomp($_);
    if ($_=~/\>/){
	@l=split(/\s+/,$_);
	@t=split(/\|/,$l[0]);
	$trs{$t[1]}++;
    }
}

#GF000001 Mmu ENSMUSG00000090353Gm17555OG0000001O=0,R=0
open (CL, "$cl");
while (<CL>){
    chomp($_);
    @l=split(/\t/,$_);
    $nt{$l[0]}{$l[1]}+=$trs{$l[2]};
    $cid{$l[0]}=1;
}
my @k=sort(keys(%cid));
my $el;
open (NALN, ">Num_alns_by_cl_$sp1-$sp2.log");
foreach $el (@k){
    $naln=$nt{$el}{$sp1}*$nt{$el}{$sp2};
    $taln+=$naln;
    if ($naln>$N){  
	print NALN "$el\t$naln\tThis cluster exceed $N alignments: Reduce number of total transcripts in exint file for this cluster or increase number of alignments!!!\n";
	$warn{$el}=1;
	$w++;
    }
    else { print NALN "$el\t$naln\n"; }
    $nalns{$el}=$naln;
    if (!$max){ $max=$naln; }
    elsif ($naln>$max){$max=$naln; }
}
print NALN "\nTotal alignments\t$taln\t$max\nNumber of clusters with warnings\t$w\n";
if ($w){ print STDERR "\n\nWARNING!!! $w clusters exceed the maximum number of alignments: $N  !!!!\n\n"; }
if ($split){
    $tmp="GCL_".$sp1."-".$sp2."_part_";
    if ($clean) { `rm $tmp*`; } 
    if ($w) { print STDERR "WARNING!!! Skipping $w clusters in the splitted files !!!\n"; }
    $sum=0;
    $part=1;
    open (CL, "$cl");
    while (<CL>){
	chomp($_);
	@l=split(/\t/,$_);
	if (!$warn{$l[0]}){
	    $out="GCL_".$sp1."-".$sp2."_part_".$part;
	    open (OUT, ">>$out"); # first part already defined
	    if ($sum<$N){
		if (!exists($prcid{$l[0]})){
		    $sum+=$nalns{$l[0]};
		    $prcid{$l[0]}=1;
		}
		print OUT "$_\n";
	    }
	    else {       
		if ($prcid{$l[0]}){
		    print OUT "$_\n";
		}
		elsif (!exists($prcid{$l[0]})){
		    $sum=$nalns{$l[0]};
		    $part++;
		    $out="GCL_".$sp1."-".$sp2."_part_".$part;
		    open (OUT, ">>$out");
		    print OUT "$_\n";
		    $prcid{$l[0]}=1;
		}
	    }
	}
    }
}


