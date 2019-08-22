#!/usr/bin/perl 
use strict;
#use warnings;
#Arguments
my ($i, $sp, $e1, $e2, $ex1, $ex2, $p1, $p2, $p3, $ie1, $ie2, $pi1, $pi2, $pe1, $pe2, $s1, $s2);
for ($i=0; $i<=$#ARGV; $i++){
	if($ARGV[$i] eq "-sp"){	
		$i++;
		$sp=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-e1"){	
		$i++;
		$e1=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-e2"){	
		$i++;
		$e2=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-ex1"){
		$i++;
		$ex1=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-ex2"){
		$i++;
		$ex2=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-p1"){
		$i++;
		$p1=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-p2"){
		$i++;
		$p2=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-p3"){
		$i++;
		$p3=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-pi1"){
		$i++;
		$pi1=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-pi2"){
		$i++;
		$pi2=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-pe1"){
		$i++;
		$pe1=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-pe2"){
		$i++;
		$pe2=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-s1"){
		$i++;
		$s1=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-s2"){
		$i++;
		$s2=$ARGV[$i]; 
	}
}

my $tmp=$sp."_annot_fake.exint";
`cat $e1 $e2 > $tmp`;
print "Joined_exint_files: $tmp\n";
$tmp="Final_exons_annot_fake_".$sp.".txt";
`cat $ex1 $ex2 > $tmp`;
print "Joined_exons_files: $tmp\n";
$tmp=$sp."_annot_fake_protein_ids_exons.txt";
`cat $p1 $p2 $p3 > $tmp`;
print "Joined_prot_ids_exon_files: $tmp\n";
$tmp=$sp."_annot_fake_protein_ids_intron_pos_CDS.txt";
`cat $pi1 $pi2 > $tmp`;
print "Joined_prot_ids_intron_files: $tmp\n";
$tmp=$sp."_annot_fake_protein_ids_exons_pos.txt";
`cat $pe1 $pe2 > $tmp`;
print "Joined_exon_pos_files: $tmp\n";
$tmp=$sp."_annot_fake_exon_sizes.txt";
`cat $s1 $s2 > $tmp`;
print "Joined_exon_sizes_files: $tmp\n";
