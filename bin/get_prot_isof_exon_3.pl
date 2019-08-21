#!/usr/bin/perl 
use strict;
#use warnings;
#Arguments
my ($i, $i1, $i2, $i3, $i4, $i5, $out, $out2, $out3, $out4, $sp);
for ($i=0; $i<=$#ARGV; $i++){
	if($ARGV[$i] eq "-sp"){	
		$i++;
		$sp=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-i1"){	
		$i++;
		$i1=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-i2"){
		$i++;
		$i2=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-i3"){
		$i++;
		$i3=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-i4"){
		$i++;
		$i4=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-i5"){
		$i++;
		$i5=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-o1"){
		$i++;
		$out=$ARGV[$i]; 
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
		print "OPTIONS: -i1 <infile1> -i2 <infile 2> -i3 <infile 3> -i4 <infile 4> -i5 <infile 5> -o1 <outfile 1> -o2 <outfile 2> -o3 <outfile 3> -o4 <outfile 4> -h help\n";
	}
}
print "##$out4\n";
#Declaration of variables
my (@line);
my %pid;
my %lprot;
my ($ex, $l, $intron);
my (@t1,@t2, @c);
my %gene;
my %psize;
my $pr;
#BL21561_evm0	BL21561	476
open (INONE, $i1); 
while (<INONE>){ 
	chomp($_);		
	@line=split(/\t/,$_);
	$pr=$line[0];
	$psize{$pr}=$line[2];
	$gene{$pr}=$line[1];	
}
close (INONE);

open (INTWO, $i2); 
my %trid;
while (<INTWO>){ 
	chomp($_);
	@line=split(/\t/,$_);
	$trid{$line[0]}=$line[1];	
}
close (INTWO);

#BL20899_cuf1|BL20899	Sc0000317	-	255783,255866-258240,258374-259202,259294
my ($tmp, $k, $r, $s);
my @n;
my %coords;
open (INTHR, $i3); 
while (<INTHR>){ 
	chomp($_);		
	@line=split(/\t/,$_);
	$s=$line[3]; $s=~s/\,/\|/g; $s=~s/\-/\,/g; $s=~s/\|/\-/g;
	@n=split(/\|/,$line[0]);
	$coords{$trid{$n[0]}}=$line[1].":".$line[3].":".$line[2];
	if ($line[3]=~/\-/){
		@c=split(/\-/,$_);
		for ($l=1; $l<scalar(@c)-1; $l++){
			my $r=$l-1; my $k=$l+1;
			@t1=split(/\,/,$c[$r]); @t2=split(/\,/,$c[$k]);
			$tmp=$c[$l];
			$tmp=~s/\,/\-/;
			if ($line[2] eq "+"){
				$ex=$line[1].":".$t1[1].",".$tmp.",".$t2[0].":".$line[2];
			}
			else {  $ex=$line[1].":".$t2[0].",".$tmp.",".$t1[1].":".$line[2]; }
			if (!$pid{$ex}){
				$pid{$ex}=$trid{$n[0]};
				$lprot{$ex}=$psize{$trid{$n[0]}};
			}
			else {  
				if ($psize{$trid{$n[0]}} > $lprot{$ex}){
					$pid{$ex}=$trid{$n[0]};
					$lprot{$ex}=$psize{$trid{$n[0]}};
				}
			}
		}
		@c=split(/\,/,$_);
		for ($l=1; $l<scalar(@c)-1; $l++){
			@t1=split(/\-/,$c[$l]);
			if ($line[2] eq "+"){
				$intron=$line[1].":".$t1[0].",".$t1[1].":".$line[2];
			}
			else {  $intron=$line[1].":".$t1[1].",".$t1[0].":".$line[2]; }
			if (!$pid{$intron}){
				$pid{$intron}=$trid{$n[0]};
				$lprot{$intron}=$psize{$trid{$n[0]}};
			}
			else {  
				if ($psize{$trid{$n[0]}} > $lprot{$intron}){
					$pid{$intron}=$trid{$n[0]};
					$lprot{$intron}=$psize{$trid{$n[0]}};
				}
			}
		}
	}	
}
close (INTHR);

#FAM13A	BlaEX0015343	Sc0000095:756974-757072	99	Sc0000095:757428,756974-757072,756243	S	Sc0000095:757428,756974-757072,756243:-	99	Sc0000095:757428-757560	Sc0000095:756974-757072	Sc0000095:756125-756243
my ($id,$g, $b, $trcoords);
my @ex;
open (OUT,">$out");
open (MISS,">$out2");
open (INFIVE, $i5); 
while (<INFIVE>){ 
	chomp($_);		
	@line=split(/\t/,$_);
	if ($pid{$line[6]}){
		$g=$gene{$pid{$line[6]}};
		$trcoords=$coords{$pid{$line[6]}};
		print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$pid{$line[6]}|$g\tAnnotated_exon_C1,A,C2*\t$trcoords\n";
	}
	else {
		print MISS "$_\n";
	}
}
close (INFIVE);
##Getting final file of exons
system ("cat $out $i4 > $out3");
##Getting file with vast exons and protein ids
my %print;
open (OUTF,">$out4");
my $mex=$sp."_not_annotated_exons.txt";
open (MEX, ">$mex");
open (IN, $out3); 
while (<IN>){ 
	chomp($_);		
	@line=split(/\t/,$_);
	if (!$print{$line[6]}){
		print OUTF "$line[6]\n";
		$print{$line[6]}=1;
	}
	if ($line[7]=~/Not_annotated_exon/){
		print MEX "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
	}	
}
close (MEX);
close (OUTF);