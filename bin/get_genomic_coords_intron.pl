#!/usr/bin/perl -w
use strict;
#Declaration of variables
#Arguments
my ($i1, $i2, $out1, $out2, $c1, $c2, $p, $t, $v, $i);
for ($i=0; $i<=$#ARGV; $i++){
	if($ARGV[$i] eq "-c1"){ 
		$i++;
		$i1=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-c2"){ 
		$i++;
		$i2=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-p"){ 
		$i++;
		$p=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-v"){ 
		$i++;
		$v=$ARGV[$i]; 
	}	
	if($ARGV[$i] eq "-t"){ 	
		$i++;
		$t=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-o1"){ 	
		$i++;
		$out1=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-o2"){ 	
		$i++;
		$out2=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-h"){	
		print "For running the program: perl get_genomic_coords_intron.pl -\n";
	}
}
##Opening GTF file##
my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
my (%pid, %tid, %ntseq);
my $b=0;
my (@l, @l1, @line, @l2, @l3, @coords);
my ($l, $size, $tmpseq, $c);
open (PROTS,"$p");
while (<PROTS>){
	chomp($_); 
	$pid{$_}=1;
}
close (PROTS);
if ($v){
	open (VAST,"$v");
	while (<VAST>){
		chomp($_); 
		$pid{$_}=1;
	}
	close (VAST);
}
#ENSBTAT00000000005	ENSBTAP00000000005
open (TRS,"$t");
while (<TRS>){
	chomp($_); 
	@line=split(/\t/,$_);
	$tid{$line[0]}=$line[1];
}
close (TRS);
open (ONE, ">$out1");
open (FILE,"$i1");
while (<FILE>){
	chomp($_); 
	@line=split(/\t/,$_);
	@l=split(/\|/,$line[0]);
	if ($tid{$l[0]}){ 
		$tmp=$tid{$l[0]}."|".$l[1];
		if ($pid{$tmp}){
			if ($line[3]=~/\-/){
				@coords=split(/\,/,$line[3]);
				if ($line[2] eq "+"){
					$c=0;
					for ($l=1; $l<scalar(@coords)-1; $l++){
						$c++;
						print ONE "$tmp\tintron\t$c\t$line[1]\t$coords[$l]\t$line[2]\n";				
					}
				}
				elsif ($line[2] eq "-"){
					$c=0;
					for ($l=scalar(@coords)-2; $l>=1; $l--){
						$c++;
						print ONE "$tmp\tintron\t$c\t$line[1]\t$coords[$l]\t$line[2]\n";				
					}			
				}
			}
			else {
				print ONE "$tmp\tintronless\tNA\t$line[1]\tNA\t$line[2]\n";
			}
		}
	}
}
close (FILE);
close (ONE);
open (TWO, ">$out2");
open (FILE,"$i2");
while (<FILE>){
	chomp($_); 
	@line=split(/\t/,$_);
	$tmp=$line[0];
	if ($pid{$tmp}){
		if ($line[3]=~/\-/){
			@coords=split(/\,/,$line[3]);
			if ($line[2] eq "+"){
				$c=0;
				for ($l=1; $l<scalar(@coords)-1; $l++){
					$c++;
					print TWO "$tmp\t$line[1]\t$line[2]\tIntron_$c\t$coords[$l]\n";				
				}
			}
			elsif ($line[2] eq "-"){
				$c=0;
				for ($l=scalar(@coords)-2; $l>=1; $l--){
					$c++;
					print TWO "$tmp\t$line[1]\t$line[2]\tIntron_$c\t$coords[$l]\n";				
				}			
			}
		}
		else {
			print TWO "$tmp\t$line[1]\t$line[2]\tIntronless\n";
		}
	}
}
close (FILE);
close (TWO);
