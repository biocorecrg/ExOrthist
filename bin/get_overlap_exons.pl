#!/usr/bin/env perl
use strict;
use warnings;

#Declaration of variables and arguments
my ($i, $infile, $outfile);
for ($i=0; $i<=$#ARGV; $i++){
    if ($ARGV[$i] eq "-i"){
	$i++;
	$infile = $ARGV[$i]; 		
    }
    if ($ARGV[$i] eq "-o"){	
	$i++;
	$outfile = $ARGV[$i]; 		
    }
}
open (OUT, ">$outfile") || die "It cannot open the output file ($outfile)\n";

my %junctions;
my %ov_juncs;
my ($pos5p,$pos3p)=0;
my $count=0;
my $id;

open (INFILE, $infile) || die "It cannot open the input file (Exon_count_hits_by_sp.tab)\n";
# Format: BL00113	Sc0000002:5027516-5027714:+	Bla	18
while (<INFILE>){ 
    if ($_){
	chomp($_);
	my @line=split(/\t/,$_);
	$line[0]=~s/\s+//; # clean up of geneID
	$line[1]=~s/\s+//; # clean up of coord
	my @crs=split(/\:/,$line[1]); # chr, pos, strand
	my @pos=split(/\-/,$crs[1]);
	my $sp=$line[2];

	if (!$junctions{$line[0]}){ # if there is nothing on junctions for the geneID
	    if ($pos5p != 0) {
		if (defined $id){
		    print OUT "$ov_juncs{$id}\n" if (defined $ov_juncs{$id});
		}
	    }
	    $count++;
	    $id="OV_EX_".$sp."_".$count; # defined above since it's kept defind for the next loop
	    $ov_juncs{$id}=$id."\t".$_;
	    $junctions{$line[0]} = $line[1];
	    $pos5p=$pos[0];
	    $pos3p=$pos[1]; 
	}	
	else {
	    # this needs the input to be sorted
	    if ($pos[0]>=$pos5p && $pos[0]<=$pos3p){ ##the junction overlaps
		if ($pos[1]>$pos3p){
		    $pos3p=$pos[1];
		}
		$ov_juncs{$id}.="\n".$id."\t".$_;
	    }
	    else { ### the junction does not overlap, print all the information from previous junctions
		$pos5p=$pos[0];
		$pos3p=$pos[1];
		print OUT "$ov_juncs{$id}\n";
		$count++;
		$id="OV_EX_".$sp."_".$count;
		$ov_juncs{$id}=$id."\t".$_;
	    }
	}
    }
}
close (INFILE);

###printing the last junction
print OUT "$ov_juncs{$id}\n";
close OUT;



