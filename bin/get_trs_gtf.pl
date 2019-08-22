#!/usr/bin/env perl 
use strict;
#Declaration of variables
#Arguments
my ($i, $infile1, $gtf, $out1, $out2, $out3, $out4);
for ($i=0; $i<=$#ARGV; $i++){
	if($ARGV[$i] eq "-gtf"){	
		$i++;
		$infile1=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-o1"){
		$i++;
		$out1=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-o2"){
		$i++;
		$out2=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-o3"){
		$i++;
		$out3=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-o4"){
		$i++;
		$out4=$ARGV[$i]; 
	}	
	if($ARGV[$i] eq "-h"){
		print "OPTIONS: -gtf <gtf file> -o1 <outfile 1> -o2 <outfile 2> -o3 <outfile 3>-h help\n";
	}
}
my (@line, @coords1, @coords2, @gene, @l1, @l2, @l3);
my ($l, $grep, $gid, $int, $r, $prot, $tmp, $trid, $n);
my (%tr1, %tr2, %tpid);
my %prot;
my (%coords1, %coords2);
my (%strand1, %strand2);
my %intron;
my (%chr1, %chr2);
open (INFILEONE, $infile1); 
#1	protein_coding	exon	860260	860328	.	+	.	gene_id "ENSG00000187634"; transcript_id "ENST00000420190"; exon_number "1"; gene_name "SAMD11"; gene_biotype "protein_coding"; transcript_name "SAMD11-011"; exon_id "ENSE00001637883";
while (<INFILEONE>){ 
	if($_){ 
		chomp($_);
		@line=split(/\t/,$_);
		if ($line[2] eq "exon"){
			@l1=split(/\"/,$line[8]);
			$gid=$l1[1];
			$tmp="transcript_id";
			@l2=split(/$tmp/,$line[8]);
			@l3=split(/\"/,$l2[1]);
			$trid=$l3[1];
			$strand1{$trid."|".$gid}=$line[6];
			$chr1{$trid."|".$gid}=$line[0];
			if (!$tr1{$trid."|".$gid}){
				$tr1{$trid."|".$gid}=$line[3].",".$line[4];
			}
			else {
				$tr1{$trid."|".$gid}.="-".$line[3].",".$line[4];
			}	
		}
		elsif ($line[2] eq "CDS"){
			@l1=split(/\"/,$line[8]);
			if ($line[8]=~/protein_id/) { $tmp="protein_id"; }
			else { $tmp="transcript_id"; }
			@l2=split(/$tmp/,$line[8]);
			@l3=split(/\"/,$l2[1]);
			$prot=$l3[1];
			$tmp="transcript_id";
			@l2=split(/$tmp/,$line[8]);
			@l3=split(/\"/,$l2[1]);
			$trid=$l3[1];
			$tpid{$trid."\t".$prot}=1;
			$strand2{$prot."|".$gid}=$line[6];
			$chr2{$prot."|".$gid}=$line[0];
			if (!$tr2{$prot."|".$gid}){ ##modifying when the first exon of the CDS is not in phase 0
				if ($line[7] == 1){
					if ($line[6] eq "+"){ 
						$n=$line[3]+1;
						$tr2{$prot."|".$gid}=$n.",".$line[4];
					} else { 
						$n=$line[4]-1;
						$tr2{$prot."|".$gid}=$line[3].",".$n;
					}
				}
				elsif ($line[7] == 2) {
					if ($line[6] eq "+"){ 
						$n=$line[3]+2;
						$tr2{$prot."|".$gid}=$n.",".$line[4];
					} else { 
						$n=$line[4]-2;
						$tr2{$prot."|".$gid}=$line[3].",".$n;
					}
				}
				else {
					$tr2{$prot."|".$gid}=$line[3].",".$line[4];
				}
			}
			else {
				$tr2{$prot."|".$gid}.="-".$line[3].",".$line[4];
			}	
		}
	}
}
close (INFILEONE);

my @keys=keys(%tr1);
my $element;
open (OUTFILE, ">$out1");
foreach $element (@keys){
	if ($strand1{$element} eq "+"){
		print OUTFILE "$element\t$chr1{$element}\t$strand1{$element}\t$tr1{$element}\n";
	}
	else {
		@coords1=split(/\-/,$tr1{$element});	
		print OUTFILE "$element\t$chr1{$element}\t$strand1{$element}\t";
		for ($l=scalar(@coords1)-1; $l>=1; $l--){
			print OUTFILE  "$coords1[$l]-";
		}
		print OUTFILE "$coords1[0]\n";
	}
}

my @keys2=keys(%tr2);
my $element2;
my %print;
open (OUTFILE, ">$out2");
open (IDS, ">$out4");
foreach $element2 (@keys2){
	if ($strand2{$element2} eq "+"){
		print OUTFILE "$element2\t$chr2{$element2}\t$strand2{$element2}\t$tr2{$element2}\n";
		if (!$print{$element2}){
			print IDS "$element2\n";
			$print{$element2}=1;
		}
	}
	else {
		@coords2=split(/\-/,$tr2{$element2});	
		print OUTFILE "$element2\t$chr2{$element2}\t$strand2{$element2}\t";
		if (!$print{$element2}){
			print IDS "$element2\n";
			$print{$element2}=1;
		}
		for ($l=scalar(@coords2)-1; $l>=1; $l--){
			print OUTFILE  "$coords2[$l]-";
		}
		print OUTFILE "$coords2[0]\n";
	}
}
my @keys3=keys(%tpid);
@keys=sort(@keys3);
my $tp;
open (OUTFILE, ">$out3");
foreach $tp (@keys3){	
		print OUTFILE "$tp\n";
}













