#!/usr/bin/env perl
use warnings;
use strict;

my $input_file1 = $ARGV[0]; # filtered_best_scored_EX_matches_by_targetgene.tab
my $output_file = $ARGV[1]; # EX_matches_count_by_species.tab
my $input_file2 = $ARGV[2]; # Liftover exons (optional)
 
my %ex; # count of exons

open (OUT, ">$output_file") || die "It cannot open the output file ($output_file)\n";

### Processes the first input file
# Format: BL00000	Sc0000095:745259-745367:-	SMAR011687	JH432088:105973-106084:-	Bla	Sma
open (INONE, $input_file1) || die "It cannot open the first input file\n";
while (<INONE>){
    chomp($_); 
    my @l=split(/\t/,$_);
    if ($l[0]!~/GeneID_sp/){ # ignores header
	my $e1=$l[0]."\t".$l[1]."\t".$l[4];
	my $e2=$l[2]."\t".$l[3]."\t".$l[5]; 
	$ex{$e1}++;
	$ex{$e2}++;
    }
}
close INONE;

### Processes the liftover pairs if provided
# Format: ENSG00000171055  chr2:36552056-36552268:-   ENSMUSG00000056121   chr17:78377717-78377890:-  Sp1  Sp2
if (defined $input_file2){
    open (INTWO, $input_file2) || die "It cannot open the Liftover file\n";
    while (<INTWO>){
	chomp($_); 
	my @l=split(/\t/,$_);
	if ($l[0]!~/GeneID_sp/){ # ignores header
	    my $e1=$l[0]."\t".$l[1]."\t".$l[4];
	    my $e2=$l[2]."\t".$l[3]."\t".$l[5];	
	    $ex{$e1}++;
	    $ex{$e2}++;
	    $ex{$e1}=$ex{$e1}+1000; # to give more weight to bonafide exons
	    $ex{$e2}=$ex{$e2}+1000; # to give more weight to bonafide exons
	}
    }
    close INTWO;
}

my @keys=keys(%ex);
@keys=sort (@keys);
foreach my $el (@keys){
    print OUT "$el\t$ex{$el}\n";
}
close OUT;








