#!/usr/bin/perl -w
use strict;
#Declaration of variables
#Arguments
my ($i, $gtf, $p, $v,$outfile);
for ($i=0; $i<=$#ARGV; $i++){
	if($ARGV[$i] eq "-GTF"){ 
		$i++;
		$gtf=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-i"){ ##proteins ids annotated exons
		$i++;
		$p=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-v"){ ##vastdb protein ids exons
		$i++;
		$v=$ARGV[$i]; 		
	}	
	if($ARGV[$i] eq "-out"){ 	
		$i++;
		$outfile=$ARGV[$i]; 		
	}
}
##Opening proteins file
my (%pid, %protein);
my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m, $ex);
my (%header, %ntseq);
my $b=0;
my (@l, @l1, @line, @l2, @l3, @s);
open (PROTS, "$p");
while (<PROTS>){
	chomp($_);
	$protein{$_}=1;
}
close (PROTS);
if ($v){
	open (VAST, "$v");
	while (<VAST>){
		chomp($_);
		$protein{$_}=1;
	}
	close (VAST);
}

##Opening GTF file##
#Scaffold1348	protein_coding	CDS	153968	154192	.	+	.	gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "4"; protein_id "WHL22.100348.0";
#Scaffold1348	protein_coding	exon	162899	162997	.	+	.	gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "5";
#Scaffold1348	protein_coding	CDS	162899	162997	.	+	.	gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "5"; protein_id "WHL22.100348.0";
#Scaffold1348	protein_coding	exon	163794	163930	.	+	.	gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "6";
my %rseq;
my ($res, $size, $new, $id, $e);
my %fex; ##saving the phase of the first exon
open (OUT, ">$outfile");
open (GTF,"$gtf");
while (<GTF>){
	chomp($_); 
	@line=split(/\t/,$_);
	if ($line[2] eq "CDS"){
		@l1=split(/\"/,$line[8]);
		$gid=$l1[1]; ##change here
		if ($line[8]=~/protein_id/){ $tmp="protein_id"; }
		else { $tmp="transcript_id";  }
		@l2=split(/$tmp/,$line[8]);
		@l3=split(/\"/,$l2[1]);
		$prot=$l3[1];
		$id=$prot."|".$gid;
		$tmp="exon_number";
		@l2=split(/$tmp/,$line[8]);
		@l3=split(/\"/,$l2[1]);
		$ex=$l3[1];
		if ($protein{$id}){
			if (!$pid{$prot}){
				$e=1;
				$pid{$prot}=1;
				$m=0;
				if ($line[7] ne "."){  $phase=$line[7]; $fex{$prot}=$phase; $m=$phase; } ##
				else { 
					$phase=0; 
					$fex{$prot}=$phase;
				}
				$nuc=0;
				$size=$line[4]-$line[3]+1;
				$size=$size-$m;
				$res=$size%3;
				$nuc=$nuc+$size;
				$aa=int(1+$nuc/3); 
				print OUT "$prot|$gid\texon_$e\t1-$aa\t$line[0]\t$line[3]-$line[4]\t$line[6]\n"; ###exon numbering according to CDS not to transcript, otherwise $ex variable shoud be use
			}
			else {
				$e++;
				$new=$aa+1;
				if ($line[7] ne "."){  $phase=$line[7]; }
				else {  
					if ($res==0){ $phase=0; } elsif ($res==1 ){ $phase=2; } elsif ($res==2) { $phase=1; } 
				}
				$size=($line[4]-$line[3])+1;
				$nuc=$nuc+$size;
				$size=$size+$res;		
				$res=$size%3;
				$aa=int(1+$nuc/3);
				print OUT "$prot|$gid\texon_$e\t$new-$aa\t$line[0]\t$line[3]-$line[4]\t$line[6]\n"; 						
			}
		}
	}
}
close (GTF);

