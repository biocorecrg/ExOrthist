#!/usr/bin/perl -w
use strict;
#Declaration of variables
#Arguments
my ($i, $gtf, $genome, $bin, $ref, $sp, $p, $f);
for ($i=0; $i<=$#ARGV; $i++){
	if($ARGV[$i] eq "-GTF"){ ##folder with GTFs
		$i++;
		$gtf=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-G"){  ##folder with genomes
		$i++;
		$genome=$ARGV[$i]; 		
	}	
	if($ARGV[$i] eq "-REF"){ ##folder with reference proteins	
		$i++;
		$ref=$ARGV[$i]; 		
	}	
	if($ARGV[$i] eq "-S"){  ##list of sp, separated by comma	
		$i++;
		$sp=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-F"){  ##for generating exint file -F 1, for reference file -F 	
		$i++;
		$f=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-bin"){  ##bin directory                                                                            
	    $i++;
	    $bin=$ARGV[$i];
        }
	if($ARGV[$i] eq "-h"){	
		print "For running the program: perl get_ref_proteins.pl -GTF <folder GTF files> -G <folder with Genome sequences> -S <species list separated by comma> -F <1 or 2> -bin <bin_directory>
\n";
		print "-F 1 for generating the exint file; -F 2 for generating the reference file\n"
	}
}
my @species=split(/\,/,$sp);
my (@l);
my ($fgtf, $fref, $fgen, $out1, $out2, $rp, $j, $sfile, $n, $s, $hq);
	for ($j=0; $j<scalar(@species); $j++){
		$fgtf=$gtf.$species[$j]."_annot.gtf";	
		$fgen=$genome.$species[$j]."_gDNA.fasta";
		$rp=$species[$j]."_ref_proteins.txt";
		$out1=$species[$j].".exint";
		$sfile=$species[$j]."_prot_sizes.txt";
		$out2=$species[$j]."_ref_proteins.exint";
		$hq=$species[$j]."_HQ_prots.txt";
		#1) Getting exint file
		if ($f==1){
			`perl $bin/get_exint_file.pl -GTF $fgtf -G $fgen -out $out1`;
		}
		else { 
		#2) Getting sizes of annotated proteins
		open (EXINT, "$out1");
		open (SIZES, ">$sfile");
		while (<EXINT>){
			if($_=~/\>/){ 
				@l=split(/\s+/,$_);
				$n=$l[0];
				$n=~s/\>//; $n=~s/\|/\t/;
			}
			else { chomp($_); $s=length($_); print SIZES "$n\t$s\n";  }				
		}
		#3) Getting reference proteins
		`perl $bin/get_ref_prot.pl $sfile $rp`;
		`perl $bin/get_ref_prot_exint_file.pl $rp $out1 $out2`;
	}
}








