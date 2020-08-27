#!/usr/bin/env perl
use warnings;
use strict;

my $input_file=$ARGV[0]; # scores all introns btw pair of species (all merged Final_aln_scores_Sp1_Sp2.txt)
my $outf=$ARGV[1]; # output folder => Sp1_Sp2
my ($info, $id, $r, $k, $pid, $s1, $s2);
my $b=0;
my (@i1, @i2, @i1b, @i2b, @line, @t1, @t2, @t3, @t4, @name, @l1, @l2);
my (%sintron);
my (%psc, %esc, %isc, %cid, %le, %sps);

open (IN, $input_file) || die "It cannot open File with Aln scores (Final_aln_scores_Sp1_Sp2.txt)\n";
### Format (multi-line: various formats in a single file
# (1) EXONS   # Exon	PAC30764619|Ocbimv22009530m	exon_1	1-33	Scaffold1679	955644-955741	+	ENSP00000403502|ENSG00000178035	1-33	78.79%	78.79%	0	0.00%	ENSP00000403502|ENSG00000178035	exon_1	1-33	chr3	49066686-49066783	-
# (2) INTRONS # CID:6948	PAC30764619|Ocbimv22009530m	intron_1	33	phase_1	ENSP00000403502|ENSG00000178035	intron_1	33	phase_1	conserved;u,u
# (3) PROTEIN # CID:6948	Protein	Query:PAC30764619|Ocbimv22009530m	Subject:ENSP00000403502|ENSG00000178035	43.07%	81.36%	38.33%	0.29
while (<IN>){
    chomp($_); 
    @line=split(/\t/,$_);
    if ($line[2]=~/exon/){ # (1)
	$sps{$line[1]}=$line[8];
	$sps{$line[4]}=$line[9];
	if ($line[5] ne "NO_EXON_ALN"){
	    $id=$line[1]."#".$line[2]."#".$line[4];
	    $esc{$id}=$line[7];
	}
	$le{$line[1]."#".$line[4]}=$line[2]; ##saving last exon
    }
    elsif ($line[2]=~/intron/){ # (2)
	$id=$line[1]."#".$line[2]."#".$line[4];
	if ($line[5] eq "NO_ALN"){ }
	elsif ($line[5] eq "NO_INTRON"){ $isc{$id}=$line[5]; }
	else {
	    if ($line[7]<0){ $isc{$id}=-0.25; }
	    else {
		my $ps;
		for ($ps=1; $ps<=10; $ps++){ 	
		    if($line[7]==$ps) { $isc{$id}=($ps*0.025); } ##change here now is max score = 10
		}
	    }
	}
    } 
    ### (3) not used
}

open (IN, $input_file) || die "It cannot open again File with Aln scores (Final_aln_scores_Sp1_Sp2.txt)\n";
my $scorefile=$outf."/All_scores_by_exon.txt";
open (OUTONE, ">$scorefile");
print OUTONE "CID\tFeature\tExon_prot_loc\tProt_query\tExon_number_query\tExon_coords_query\tProt_Subject\tExon_number_subject\tExon_coords_subject";
print OUTONE "\tScore_C1\tScore_I1\tScore_A\tScore_I2\tScore_C2\tTotal_exon_score\tSp_query\tSp_subject\n";
my ($e1, $i1, $e2, $i2, $c);
my ($nid, $pscore, $sC1, $sC2, $sA, $sI1, $sI2); 
my $TSC;
my %totalSC;
my %hex;
my (@g1,@g2);
my $idex;
my $string;
while (<IN>){
    chomp($_); 
    @line=split(/\t/,$_);
    if ($line[2]=~/exon/){
	($sC1, $sC2, $sA, $sI1, $sI2, $TSC)=(0,0,0,0,0,0);
	$nid=$line[1]."#".$line[4];
	@g1=split(/\|/,$line[1]);
	@g2=split(/\|/,$line[4]);
	$idex=$g1[1]."\t".$line[3]."\t".$g2[1];	
	$c=$line[0];
	$string="";
	if ($line[2] eq "exon_1"){ ##N-terminal exon, we evaluate only I2 and C2;
	    $e2=$line[1]."#exon_2#".$line[4];
	    $i2=$line[1]."#intron_1#".$line[4];
	    if ($isc{$i2}) { if ($isc{$i2} eq "NO_INTRON"){ $sI2=0; } else { $sI2=$isc{$i2};}  } else { $sI2=-0.25 }
	    if ($esc{$e2}) { $sC2=($esc{$e2} * 0.15); } # prev = 0.16
	    if ($line[5] ne "NO_EXON_ALN") { 
		$sA=($line[7] * 0.20); # prev = 0.18
		$TSC=$sA+$sI2+$sC2; ##Total Score
		print OUTONE "$c\tExon\tN_terminal\t$line[1]\t$line[2]\t$line[3]\t";
		print OUTONE "$line[4]\t$line[5]\t$line[6]\tNA\tNA\t$sA\t$sI2\t$sC2\t$TSC\t$line[8]\t$line[9]\n";
		$string=$c."\tExon"."\t"."N_terminal\t".$line[1]."\t".$line[2]."\t".$line[3]."\t";
		$string.=$line[4]."\t".$line[5]."\t".$line[6]."\tNA\tNA\t".$sA."\t".$sI2."\t".$sC2."\t".$TSC."\t".$line[8]."\t".$line[9];
	    }
	    else { 
		$sA=-0.1; ##Nothing aligned
		$TSC=$sA+$sI2+$sC2+0.41; ##Total Score  (to sum up to 1, the missing exon and intron score are added)
		print OUTONE "$c\tExon\tN_terminal\t$line[1]\t$line[2]\t$line[3]\t";
		print OUTONE "$line[4]\t$line[5]\tNA\tNA\tNA\t$sA\t$sI2\t$sC2\t$TSC\t$line[8]\t$line[9]\n";
		$string=$c."\tExon"."\t"."N_terminal\t".$line[1]."\t".$line[2]."\t".$line[3]."\t";
		$string.=$line[4]."\t$line[5]\tNA\tNA\tNA\t".$sA."\t".$sI2."\t".$sC2."\t".$TSC."\t".$line[8]."\t".$line[9];
	    }
	}
	elsif ($line[2] eq $le{$nid}) { ##C-terminal exon, we evaluate only I1 and C1
	    @l1=split(/\_/,$line[2]);
	    $e1=$line[1]."#exon_".($l1[1]-1)."#".$line[4];
	    $i1=$line[1]."#intron_".($l1[1]-1)."#".$line[4];
	    if ($isc{$i1}) { if ($isc{$i1} eq "NO_INTRON"){ $sI1=0; } else { $sI1=$isc{$i1} } }  else{ $sI1=-0.25; }
	    if ($esc{$e1}) { $sC1=($esc{$e1} * 0.15); } # prev = 0.16
	    if ($line[5] ne "NO_EXON_ALN") { 
		$sA=($line[7] * 0.20); # prev = 0.18
		$TSC=$sC1+$sI1+$sA+0.41; ##Total Score  (to sum up to 1, the missing exon and intron score are added); ##Total Score
		print OUTONE "$c\tExon\tC_terminal\t$line[1]\t$line[2]\t$line[3]\t";
		print OUTONE "$line[4]\t$line[5]\t$line[6]\t$sC1\t$sI1\t$sA\tNA\tNA\t$TSC\t$line[8]\t$line[9]\n";
		$string=$c."\t"."Exon"."\t"."C_terminal\t".$line[1]."\t".$line[2]."\t".$line[3]."\t";
		$string.=$line[4]."\t".$line[5]."\t".$line[6]."\t".$sC1."\t".$sI1."\t".$sA."\tNA\tNA\t".$TSC."\t".$line[8]."\t".$line[9];
	    }
	    else { 
		$sA=-0.1; 
		$TSC=$sC1+$sI1+$sA;
		print OUTONE "$c\tExon\tC_terminal\t$line[1]\t$line[2]\t$line[3]\t";			
		print OUTONE "$line[4]\t$line[5]\tNA\t$sC1\t$sI1\t$sA\tNA\tNA\t$TSC\t$line[8]\t$line[9]\n";
		$string=$c."\t"."Exon"."\t"."C_terminal\t".$line[1]."\t".$line[2]."\t".$line[3]."\t";
		$string.=$line[4]."\t$line[5]\tNA\t".$sC1."\t".$sI1."\t".$sA."\tNA\tNA\t".$TSC."\t".$line[8]."\t".$line[9];
	    }
	}
	else { ##internal exon
	    @l1=split(/\_/,$line[2]);
	    $e1=$line[1]."#exon_".($l1[1]-1)."#".$line[4];
	    $i1=$line[1]."#intron_".($l1[1]-1)."#".$line[4];
	    $e2=$line[1]."#exon_".($l1[1]+1)."#".$line[4];
	    $i2=$line[1]."#intron_".($l1[1])."#".$line[4];
	    if ($isc{$i1}) { if ($isc{$i1} eq "NO_INTRON"){ $sI1=0; } else { $sI1=$isc{$i1}; }  } else { $sI1=-0.25; }
	    if ($esc{$e1}) { $sC1=($esc{$e1} * 0.15); } # prev = 0.16
	    if ($isc{$i2}) { if ($isc{$i2} eq "NO_INTRON") { $sI2=0; } else { $sI2=$isc{$i2}; } } else{ $sI2=-0.25; }
	    if ($esc{$e2}) { $sC2=($esc{$e2} * 0.15); } # prev = 0.16 
	    if ($line[5] ne "NO_EXON_ALN") { 
		$sA=$line[7] * 0.20; # prev = 0.18
		$TSC=$sC1+$sI1+$sA+$sI2+$sC2;			
		print OUTONE "$c\tExon\tInternal\t$line[1]\t$line[2]\t$line[3]\t";
		print OUTONE "$line[4]\t$line[5]\t$line[6]\t$sC1\t$sI1\t$sA\t$sI2\t$sC2\t$TSC\t$line[8]\t$line[9]\n";
		$string=$c."\t"."Exon"."\t"."Internal\t".$line[1]."\t".$line[2]."\t".$line[3]."\t";
		$string.=$line[4]."\t".$line[5]."\t".$line[6]."\t".$sC1."\t".$sI1."\t".$sA."\t".$sI2."\t".$sC2."\t".$TSC."\t".$line[8]."\t".$line[9];
	    }
	    else { 
		$sA=-0.1; 
		$TSC=$sC1+$sI1+$sA+$sI2+$sC2;
		print OUTONE "$c\tExon\tInternal\t$line[1]\t$line[2]\t$line[3]\t";			
		print OUTONE "$line[4]\t$line[5]\tNA\t$sC1\t$sI1\t$sA\t$sI2\t$sC2\t$TSC\t$line[8]\t$line[9]\n";
		$string=$c."\t"."Exon"."\t"."Internal\t".$line[1]."\t".$line[2]."\t".$line[3]."\t";
		$string.=$line[4]."\t$line[5]\tNA\t".$sC1."\t".$sI1."\t".$sA."\t".$sI2."\t".$sC2."\t".$TSC."\t".$line[8]."\t".$line[9];			
	    }
	}
	### Here is where the BEST hit per GENE is selected
	if ($TSC==0) { $TSC=0.000001; }
	if (!$totalSC{$idex}){
	    $totalSC{$idex}=$TSC;			
	    $hex{$idex}=$string;
	}
	elsif ($TSC>$totalSC{$idex}){
	    $totalSC{$idex}=$TSC;			
	    $hex{$idex}=$string;
	}
    }
}
close (OUTONE);

### Creates second output file
my $bhfile=$outf."/Best_score_hits_exons.txt";
open (OUT, ">$bhfile");
my @keys=keys(%hex);
@keys=sort(@keys);
my $el;
##PRINT HEADER
print OUT "CID\tFeature\tExon_prot_loc\tProt_query\tExon_number_query\tExon_coords_query\tProt_Subject\tExon_number_subject\tExon_coords_subject";
print OUT "\tScore_C1\tScore_I1\tScore_A\tScore_I2\tScore_C2\tTotal_exon_score\tSp_query\tSp_subject\n";
foreach $el (@keys){
    print OUT "$hex{$el}\n";
}
close (OUT);













