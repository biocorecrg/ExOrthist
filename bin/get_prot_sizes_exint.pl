#!/usr/bin/env perl -w

use strict;

#Declaration of variables
#Arguments

my ($i, $gtf, $exint, $sfile, $sp, $p, $f, $l, $n, $s);


for ($i=0; $i<=$#ARGV; $i++){

	if($ARGV[$i] eq "-e"){ ##folder with GTFs
		$i++;
		$exint=$ARGV[$i]; 
		
	}

	if($ARGV[$i] eq "-o"){  ##for generating exint file -F 1, for reference file -F 2
	
		$i++;
		$sfile=$ARGV[$i]; 
		
	}



}

open (EXINT, "$exint");

open (SIZES, ">$sfile");

my @l;

while (<EXINT>){

		if($_=~/\>/){ 

			@l=split(/\s+/,$_);

			$n=$l[0];

			$n=~s/\>//; $n=~s/\|/\t/;

		}

		else { chomp($_); $s=length($_); print SIZES "$n\t$s\n";  }
				
}


