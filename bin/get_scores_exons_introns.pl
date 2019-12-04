#!/usr/bin/env perl
use warnings;
use strict;


##SCRIPT FOR SCORING EXONS BTW PAIR OF SPECIES##

##
my $s1=$ARGV[0]; ##species 1 
my $s2=$ARGV[1]; ##species 2
my $i1=$ARGV[2]; ##aligned protein sp1 sp2
my $i2=$ARGV[3]; ##exon scores
my $i3=$ARGV[4]; ##intron scores
my $i4=$ARGV[5]; ##exint file sp1
my $i5=$ARGV[6]; ##exint file sp1
my $i6=$ARGV[7]; ##intron coords sp1
my $i7=$ARGV[8]; ##intron coords sp1
my $outf=$ARGV[9];

##protein ids exons species 1
open (IN,"$i1") || die "Missing aligned proteins file";  
#13	Protein	ENSDARP00000115151|ENSDARG00000056001	ENSP00000351113|ENSG00000075539	69.33	67.82	57.95	0.57	Dre	Hsa
my (%pid1);
my $id;
my (%p1, %p2, %aln);
my (@t1, @t2, @l, @line);
my ($pint1,$pint2);
while (<IN>){
	chomp($_);
	@l=split(/\t/,$_);
	$pid1{$l[2]}=$l[0];
	$pid1{$l[3]}=$l[0];
	$aln{$l[2]."#".$l[3]}=$_;
}
##PROCESSING EXINT FILES
my $f1=$i4;
my $f2=$i5;
my $m=0; 
my @tmp;
my ($e,$i);
my %introns;
my %exons;
my %intpos;
#>ENSP00000379387|ENSG00000185219  144.0 200.2 232.0 274.2 311.2
open (EXONE, "$f1") || die "Missing exint file species 1";
while (<EXONE>){
	chomp($_);
	if ($_=~/\>/){
		@l=split(/\s+/,$_);
		$l[0]=~s/\>//;
		if ($pid1{$l[0]}){
			if (scalar(@l)==1){
				$exons{$l[0]}=1;
			}
			else {  
 				$e=scalar(@l);
				$exons{$l[0]}=$e;
			}	
		}
	}			
}
close (EXONE); 

open (EXTWO, "$f2")|| die "Missing exint file species 2";
while (<EXTWO>){
	chomp($_);
	if ($_=~/\>/){
		@l=split(/\s+/,$_);
		$l[0]=~s/\>//;
		if ($pid1{$l[0]}){
			if (scalar(@l)==1){
				$exons{$l[0]}=1;
			}
			else {  
				$e=scalar(@l);
				$exons{$l[0]}=$e;
			}	
		}
	}			
}
close (EXTWO);
##GETTING GENOMIC COORDS OF INTRONS
#BL22196_evm0|BL22196	Sc0000025	-	Intron_1	366299-366745
my @k;
my ($posint1, $posint2);
open (IN,"$i6")|| die "Missing intron pos species 1"; ###exon intron scores
while (<IN>){
	chomp($_);
    	@line=split(/\t/,$_);
    	$line[3]=~s/I/i/;
    	if ($line[3] ne "Intronless" ){
	    @k=split(/\-/,$line[4] && $line[4]=~/\-/);
	    $posint1=$k[0]+1; $posint2=$k[1]-1;
	    $intpos{$line[0]."_".$line[3]}=$line[1].":".$posint1."-".$posint2.":".$line[2];
    	}
}
open (IN,"$i7")|| die "Missing intron pos species 2"; ###exon intron scores
while (<IN>){
	chomp($_);
	@line=split(/\t/,$_);
    	$line[3]=~s/I/i/;
    	if ($line[3] ne "Intronless" && $line[4]=~/\-/){
	    @k=split(/\-/,$line[4]);
            $posint1=$k[0]+1; $posint2=$k[1]-1;
            $intpos{$line[0]."_".$line[3]}=$line[1].":".$posint1."-".$posint2.":".$line[2];
    	}

}
##PROCESSING EXON POSITION FILES
my %exsc;
my %insc;
my $rs1;
my %pos;
my @n; 
my %coords;
my %sps;
open (IN,"$i2")|| die "Missing exon intron scores file"; ###exon intron scores
while (<IN>){
	chomp($_);	
	@line=split(/\t+/,$_);##change of format
	#Exon	ENSMUSP00000000028|ENSMUSG00000000028	exon_10	236-275	chr16	18795918-18796037	-	ENSP00000397434|ENSG00000093009	NO_EXON_ALN
	#13	ENSDARP00000044293|ENSDARG00000056001	exon_28	1152-1188	chr10	33997813-33997923	+	1	ENSP00000261575|ENSG00000073910	1-3	2.70	5.41	34	91.89	ENSP00000261575|ENSG00000073910	exon_1	1-54	chr13	32768274-32768434	+	Dre	Hsa
	if ($line[2]=~/exon/ && $line[7]==1){
		$coords{$line[1]."#".$line[2]}=$line[4].":".$line[5].":".$line[6];	
		$coords{$line[14]."#".$line[15]}=$line[17].":".$line[18].":".$line[19];
	        $sps{$line[1]}=$line[20];
	        $sps{$line[8]}=$line[21];
			$id=$line[1]."#".$line[2]."#".$line[8];		
		if ($line[11]>0){ $line[11]=sprintf("%.2f",($line[11]/100));  
			$exsc{$id}=$line[0]."\t".$line[1]."\t".$line[2]."\t".$line[4].":".$line[5].":".$line[6]."\t".$line[8]."\t".$line[15]."\t".$line[17].":".$line[18].":".$line[19]."\t".$line[11]."\t".$line[20]."\t".$line[21];
		}
		else { }
	}
}
#GF000006	ENSP00000268043|ENSG00000140451	intron_4	324	1	ENSMUSP00000049046|ENSMUSG00000041064	intron_4	333	1	0	z,z	25	Hs2	Mm2
open (IN,"$i3")|| die "Missing exon intron scores file";
while (<IN>){
	chomp($_);	
	@line=split(/\t+/,$_);
	if ($line[2]=~/intron/){
		$id=$line[1]."#".$line[2]."#".$line[5];
		$pint1=$line[1]."_".$line[2];
		#$pint2=$line[5]."_".$line[6];
		#print "**$id\n";
		if ($intpos{$line[5]."_".$line[6]}){ $pint2=$intpos{$line[5]."_".$line[6]}; }
		else { $pint2="NA";  }
		$insc{$id}=$line[0]."\t".$line[1]."\t".$line[2]."\t$intpos{$pint1}\t".$line[5]."\t".$line[6]."\t$pint2\t".$line[11]."\t".$line[12]."\t".$line[13];
		#$id=$line[5]."#".$line[6]."#".$line[1];
		#$insc{$id}=$line[0]."\t".$line[5]."\t".$line[6]."\t".$line[7]."\t".$line[8]."\t".$line[1]."\t".$line[2]."\t".$line[3]."\t".$line[4]."\t".$line[9];
		#print "**$id\n";
	}
}
my @keys;
@keys=keys(%aln);
my $el;
my $prots;
my $nex;
my $id2;
my $cr;
my $cl;
my $ne1;
my $ne2;
#$outf="Final_aln_scores_".$s1."_".$s2.".txt";
##HEADER##
open (OUT, ">$outf");
print OUT "CID\tProt_query\tFeat_query\tCoords_query\tProt_subject\tFeat_subject\tCoords_subject\tScore_feature\tSp_query\tSp_subject\n";
foreach $el (@keys){
	@tmp=split(/\#/,$el);
	$prots=$aln{$el};
	print "$prots\n";
	my $ex1=0; my $ex2=0;
	if ($exons{$tmp[0]}){ $cl=$pid1{$tmp[0]}; 
	if ($exons{$tmp[0]}==1){
		$id=$tmp[0]."#exon_1#".$tmp[1];
		if ($exsc{$id}){		
			print OUT "$exsc{$id}\n";
		}
		else {
		    $cr=$coords{$tmp[0]."#"."exon_1"};
		    print OUT "$cl\t$tmp[0]\texon_1\t$cr\t$tmp[1]\tNO_EXON_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
		}
	}
	else {
		$nex=$exons{$tmp[0]};
		for ($m=1; $m<$nex; $m++){
			$id=$tmp[0]."#exon_".$m."#".$tmp[1];
			$cr=$coords{$tmp[0]."#"."exon_".$m};
			if ($exsc{$id}){
				print OUT "$exsc{$id}\n";
			}		
			else {
				print OUT "$cl\t$tmp[0]\texon_$m\t$cr\t$tmp[1]\tNO_EXON_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
			}
			$id2=$tmp[0]."#intron_".$m."#".$tmp[1];
			#print "####$id2###\n";
			if ($insc{$id2}){
				print OUT "$insc{$id2}\n";
			}		
			else {
				$ne1=$tmp[0]."#exon_".$m."#".$tmp[1]; 
				$ne2=$tmp[0]."#exon_".($m+1)."#".$tmp[1];
				$pint1=$tmp[0]."_intron_".$m;
				if (!$exsc{$ne1} && !$exsc{$ne2}){					
					print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
				}
				elsif ($exsc{$ne1}<0.2 && !$exsc{$ne2}){
					print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
				}
				elsif ($exsc{$ne1}<0.2 && $exsc{$ne2}<0.2){
					print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
				}
				elsif (!$exsc{$ne1} && $exsc{$ne2}<0.2){
					print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
				}
				else {
					print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_INTRON\tNA\t0\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
				}
			}
		}
		$id=$tmp[0]."#exon_".$nex."#".$tmp[1];
		$cr=$coords{$tmp[0]."#exon_".$nex};
		if ($exsc{$id}){
			print OUT "$exsc{$id}\n";
		}		
		else {
			print OUT "$cl\t$tmp[0]\texon_$nex\t$cr\t$tmp[1]\tNO_EXON_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
		}
	}
	}
}
