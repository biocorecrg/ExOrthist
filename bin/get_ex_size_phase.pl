#!/usr/bin/perl -w
use strict;
#Arguments
my ($i, $gtf, $genome, $outfile);
for ($i=0; $i<=$#ARGV; $i++){
	if($ARGV[$i] eq "-GTF"){ 
		$i++;
		$gtf=$ARGV[$i]; 
	}	
	if($ARGV[$i] eq "-out"){ 	
		$i++;
		$outfile=$ARGV[$i]; 		
	}	
	if($ARGV[$i] eq "-h"){
		print "For running the program: perl get_ex_size_phase.pl -GTF <gtf_file> -out <name_outfile>\n";
	}
}

##Opening GTF file##
my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
my (%pid, %header, %ntseq);
my $b=0;
my (@l, @l1, @line, @l2, @l3, @s);
my %rseq;
my ($res, $size, $tmpseq);
my %fex; ##saving the phase of the first exon
open (GTF,"$gtf");
while (<GTF>){
	chomp($_); 
	@line=split(/\t/,$_);
	if ($line[2] eq "CDS"){
		@l1=split(/\"/,$line[8]);
		$gid=$l1[1]; ##change here 
		if ($line[8]=~/protein_id/){
			$tmp="protein_id";
		} else { $tmp="transcript_id"; }
		@l2=split(/$tmp/,$line[8]);
		@l3=split(/\"/,$l2[1]);
		$prot=$l3[1];
		if (!$pid{$prot}){
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
			$pid{$prot}=1;
			$header{$prot}="$prot|$gid $size.$phase";
		}
		else {
			if ($line[7] ne "."){  $phase=$line[7]; }
			else {  
				if ($res==0){ $phase=0; } elsif ($res==1 ){ $phase=2; } elsif ($res==2) { $phase=1; } 
			}
			$size=($line[4]-$line[3])+1;
			$nuc=$nuc+$size;
			$header{$prot}.=" $size.$phase";
		}
	}
}
close (GTF);

my @keys=keys(%header);
my $el;
my $p;
my ($l, $triplet, $protein, $nseq, $name);
my $c=0;
open (OUT, ">$outfile");
foreach $el(@keys){ 
	$name=$header{$el};
	print OUT "$name\n";
}
close (OUT);


