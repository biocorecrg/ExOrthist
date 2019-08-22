#!/usr/bin/env perl -w
use strict;
#Declaration of variables
#Arguments

my ($i, $gtf, $genome, $outfile);
for ($i=0; $i<=$#ARGV; $i++){
	if($ARGV[$i] eq "-GTF"){ 
		$i++;
		$gtf=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-G"){ 
		$i++;
		$genome=$ARGV[$i]; 		
	}	
	if($ARGV[$i] eq "-out"){ 	
		$i++;
		$outfile=$ARGV[$i]; 		
	}	
	if($ARGV[$i] eq "-h"){	
		print "For running the program: perl get_exint_file.pl -GTF <gtf_file> -G <genome.fa> -out <name_outfile>\n";
	}
}
##Opening genome file
my (%seq);
my $chr;
my @g;
open (GENOME, "$genome");
while (<GENOME>){
	chomp($_);
	if ($_=~/\>/){ $_=~s/\>//; @g=split(/\s+/,$_); $chr=$g[0]; }
	else { $seq{$chr}.=$_;  }    
}
close (GENOME);

##Opening GTF file##
my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
my (%pid, %header, %ntseq);
my $b=0;
my (@l, @l1, @line, @l2, @l3, @s);
#Scaffold1348	protein_coding	CDS	153968	154192	.	+	.	gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "4"; protein_id "WHL22.100348.0";
#Scaffold1348	protein_coding	exon	162899	162997	.	+	.	gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "5";
#Scaffold1348	protein_coding	CDS	162899	162997	.	+	.	gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "5"; protein_id "WHL22.100348.0";
#Scaffold1348	protein_coding	exon	163794	163930	.	+	.	gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "6";
my %rseq;
my ($res, $size, $tmpseq);
my %fex; ##saving the phase of the first exon
open (GTF,"$gtf");
while (<GTF>){
	chomp($_); 
	@line=split(/\t/,$_);
	if ($seq{$line[0]}){
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
				$tmpseq=substr($seq{$line[0]},($line[3]-1),($size));
				$size=$size-$m;			
				if ($line[6] eq "-"){  $tmpseq=~tr/ATCG/TAGC/; $tmpseq=reverse($tmpseq); }
				$ntseq{$prot}=$tmpseq;
				$rseq{$prot}=$tmpseq;
				$res=$size%3;			
				$nuc=$nuc+$size;			
				$aa=int(1+$nuc/3); 
				$pid{$prot}=1;						
				$header{$prot}="$prot|$gid ";
			}
			else {
				$header{$prot}.=" $aa.";
				if ($line[7] ne "."){  $phase=$line[7]; }
				else {  
					if ($res==0){ $phase=0; } elsif ($res==1 ){ $phase=2; } elsif ($res==2) { $phase=1; } 
				}
				$size=($line[4]-$line[3])+1;
				$tmpseq=substr($seq{$line[0]},($line[3]-1),($size));
				if ($line[6] eq "-"){  $tmpseq=~tr/ATCG/TAGC/; $tmpseq=reverse($tmpseq); }
				$ntseq{$prot}.=$tmpseq;
				$rseq{$prot}.="|".$tmpseq;			
				$nuc=$nuc+$size;
				$size=$size+$res;		
				$res=$size%3;			
				$aa=int(1+$nuc/3); ##adding the residues that were left
				$header{$prot}.="$phase";
			}
		}
	}	
}
close (GTF);

###Getting protein sequence####
my %codons;
$codons{"TTT"}="F";
$codons{"TTC"}="F";
$codons{"TTA"}="L";
$codons{"TTG"}="L";
$codons{"CTT"}="L";
$codons{"CTC"}="L";
$codons{"CTA"}="L";
$codons{"CTG"}="L";
$codons{"TCT"}="S";
$codons{"TCC"}="S";
$codons{"TCA"}="S";
$codons{"TCG"}="S";
$codons{"AGT"}="S";
$codons{"AGC"}="S";
$codons{"TAT"}="Y";
$codons{"TAC"}="Y";
$codons{"TGT"}="C";
$codons{"TGC"}="C";
$codons{"TGG"}="W";
$codons{"CCT"}="P";
$codons{"CCC"}="P";
$codons{"CCA"}="P";
$codons{"CCG"}="P";
$codons{"CAT"}="H";
$codons{"CAC"}="H";
$codons{"CAA"}="Q";
$codons{"CAG"}="Q";
$codons{"CGT"}="R";
$codons{"CGC"}="R";
$codons{"CGA"}="R";
$codons{"CGG"}="R";
$codons{"AGA"}="R";
$codons{"AGG"}="R";
$codons{"ATT"}="I";
$codons{"ATC"}="I";
$codons{"ATA"}="I";
$codons{"ATG"}="M";
$codons{"ACT"}="T";
$codons{"ACC"}="T";
$codons{"ACA"}="T";
$codons{"ACG"}="T";
$codons{"AAT"}="N";
$codons{"AAC"}="N";
$codons{"AAA"}="K";
$codons{"AAG"}="K";
$codons{"GTT"}="V";
$codons{"GTC"}="V";
$codons{"GTA"}="V";
$codons{"GTG"}="V";
$codons{"GCT"}="A";
$codons{"GCC"}="A";
$codons{"GCA"}="A";
$codons{"GCG"}="A";
$codons{"GAT"}="D";
$codons{"GAC"}="D";
$codons{"GAA"}="E";
$codons{"GAG"}="E";
$codons{"GGT"}="G";
$codons{"GGC"}="G";
$codons{"GGA"}="G";
$codons{"GGG"}="G";
$codons{"TAG"}="stop";
$codons{"TAA"}="stop";
$codons{"TGA"}="stop";
my @keys=keys(%ntseq);
my $el;
my $p;
my ($l, $triplet, $protein, $nseq, $name);
my $c=0;

open (OUT, ">$outfile");
foreach $el(@keys){ 
	$nseq=$ntseq{$el};
	$name=$header{$el};
	$p=$fex{$el}; ##checking phase of the first exon
	$protein="";
	$c=0;
	for($l=$p; $l<length($nseq); $l=$l+3){
		$triplet=substr($nseq,$l,3);
		if ($c==0){
			if ($codons{$triplet}){
				if ($l>=(length($nseq)-3) && $codons{$triplet} eq "stop"){
						$c=1;	
				}
				elsif ($codons{$triplet} eq "stop"){ $c=1; } #$protein.=$codons{$triplet}; } #$c=1; }						
				else {  $protein.=$codons{$triplet};   }
			}
			else {
				$protein.="X";
			}
		}
	}
	print OUT ">$name\n$protein\n";
}
close (OUT);
