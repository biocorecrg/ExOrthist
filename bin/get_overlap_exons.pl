#!/usr/bin/env perl
#use strict;
#Declaration of variables
#Arguments
my ($i, $infile, $sp);
for ($i=0; $i<=$#ARGV; $i++){
	if($ARGV[$i] eq "-i"){
		$i++;
		$infile=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-o"){	
		$i++;
		$out=$ARGV[$i]; 		
	}
}
open (OUT, ">$out");
my @line;
my %junctions;
my @pos;
my $sum_junc=0;
my $flag=0;
my ($pos5p,$pos3p)=0;
my $gene;
my $bandera=0;
my $count=0;
my $id;
my %ov_juncs;
open (INFILE, $infile);
while (<INFILE>){ #BL00113	Sc0000002:5027516-5027714:+	Bla	18
	if($_){
		chomp($_);
		@line=split(/\t/,$_);
		$line[0]=~s/\s+//; ###saving_gene_id
		$line[1]=~s/\s+//;
		@crs=split(/\:/,$line[1]);
		@pos=split(/\-/,$crs[1]);
		$sp=$line[2];
		if(!$junctions{$line[0]}){
			if ($pos5p!=0) {
				print OUT "$ov_juncs{$id}\n";
			}
			$count++;
			$id="OV_EX_".$sp."_".$count;
			$ov_juncs{$id}=$id."\t".$_;
			$gene=$line[0];
			$junctions{$line[0]} = $line[1];
			$sum_junc=1;
			$flag=1;
			$pos5p=$pos[0];
			$pos3p=$pos[1]; 
		}	
		else {
			if ($pos[0]>=$pos5p && $pos[0]<=$pos3p){ ##the junction overlaps
				if ($pos[1]>$pos3p){
					$pos3p=$pos[1];
				}
				$sum_junc++;
				$flag=0;
				$ov_juncs{$id}.="\n".$id."\t".$_;
			}
			else{ ### the junction does not overlap, print all the information from previous junctions
				$sum_junc=0;
				$pos5p=$pos[0];
				$pos3p=$pos[1];
				$sum_junc=1;				
				$flag=0;
				$bandera=0;
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




