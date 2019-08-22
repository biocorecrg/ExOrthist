#!/usr/bin/perl 
use strict;
#Arguments
my ($i, $i1, $i2, $i3, $out, $out2);
for ($i=0; $i<=$#ARGV; $i++){
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
	if($ARGV[$i] eq "-o1"){
		$i++;
		$out=$ARGV[$i]; 
	}
	if($ARGV[$i] eq "-o2"){
		$i++;
		$out2=$ARGV[$i]; 
	}	
	if($ARGV[$i] eq "-h"){
		print "OPTIONS: -i1 <infile1> -i2 <infile 2> -i3 <infile 3> -o1 <outfile 1> -o2 <outfile 2>, -h help\n";
	}
}

#Declaration of variables
my (@line);
my %pid;
my %lprot;
my ($ex, $l, $intron, $ex2);
my (@t1,@t2, @c);
my %psize;
my $pr;
#BL21561_evm0	BL21561	476
open (INONE, $i1); 
while (<INONE>){ 
	chomp($_);
	@line=split(/\t/,$_);
	$pr=$line[0]."|".$line[1];
	$psize{$pr}=$line[2];	
}
close (INONE);

#BL20899_cuf1|BL20899	Sc0000317	-	255783,255866-258240,258374-259202,259294
my ($tmp, $k, $r, $s);
my %coords;
open (INTWO, $i2); 
while (<INTWO>){ 
	chomp($_);		
	@line=split(/\t/,$_);
	$s=$line[3]; $s=~s/\,/\|/g; $s=~s/\-/\,/g; $s=~s/\|/\-/g;
	$coords{$line[0]}=$line[1].":".$line[3].":".$line[2];
	if ($line[3]=~/\-/){
		@c=split(/\-/,$_);
		for ($l=1; $l<scalar(@c)-1; $l++){
			my $r=$l-1; my $k=$l+1;
			@t1=split(/\,/,$c[$r]); @t2=split(/\,/,$c[$k]);
			$tmp=$c[$l];
			$tmp=~s/\,/\-/;
			$ex2=$line[1].":".$tmp;
			if ($line[2] eq "+"){ 
				$ex=$line[1].":".$t1[1].",".$tmp.",".$t2[0].":".$line[2];
			}
			else {  $ex=$line[1].":".$t2[0].",".$tmp.",".$t1[1].":".$line[2]; }
			if (!$pid{$ex}){
				$pid{$ex}=$line[0];				
				$pid{$ex2}=$line[0];
				$lprot{$ex}=$psize{$line[0]};
				$lprot{$ex2}=$psize{$line[0]};
			}
			else {  
				if ($psize{$line[0]} > $lprot{$ex}){
					$pid{$ex}=$line[0];
					$pid{$ex2}=$line[0];
					$lprot{$ex}=$psize{$line[0]};
					$lprot{$ex2}=$psize{$line[0]};
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
				$pid{$intron}=$line[0];
				$lprot{$intron}=$psize{$line[0]};
			}
			else {  
				if ($psize{$line[0]} > $lprot{$intron}){
					$pid{$intron}=$line[0];
					$lprot{$intron}=$psize{$line[0]};
				}
			}
		}
	}	
}
close (INTWO);

#FAM13A	BlaEX0015343	Sc0000095:756974-757072	99	Sc0000095:757428,756974-757072,756243	S	Sc0000095:757428,756974-757072,756243:-	99	Sc0000095:757428-757560	Sc0000095:756974-757072	Sc0000095:756125-756243
my ($id,$id4, $b, $trcoords);
my @ex;
open (OUT,">$out");
open (MISS,">$out2");
open (INTHREE, $i3); 
while (<INTHREE>){ 
	chomp($_);		
	@line=split(/\t/,$_);
	if ($line[1]=~/EX/){ ##processing only exons
		if ($pid{$line[6]}){
			$trcoords=$coords{$pid{$line[6]}};
			print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$pid{$line[6]}\tAnnotated_exon_C1,A,C2\t$trcoords\n";
		}
		elsif ($pid{$line[2]}) {
			$trcoords=$coords{$pid{$line[2]}};
			print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$pid{$line[2]}\tAnnotated_exon_A\t$trcoords\n";
		}
		else {
			@t1=split(/\,/,$line[6]);	
			$id=$t1[0].",".$t1[2]; ##searching then only the intron, the exon might be not annotated
			if ($pid{$id}){
				$trcoords=$coords{$pid{$id}};
				print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$pid{$id}\tNot_annotated_exon|intron:$id\t$trcoords\n"; 
			}	
			else { print MISS "$_\n"; }
		}
	}
}
close (INTHREE);





