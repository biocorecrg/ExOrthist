#!/usr/bin/env perl
use warnings;
use strict;

my $input_file=$ARGV[0]; # scores all introns btw pair of species (all_PROT_EX_INT_aln_features_Sp1-Sp2.txt)
my $outf=$ARGV[1]; # output folder => Sp1-Sp2
my ($info, $id, $r, $k, $pid, $s1, $s2);
my $b=0;
my (@i1, @i2, @i1b, @i2b, @line, @t1, @t2, @t3, @t4, @name, @l1, @l2);
my (%sintron);
my (%psc, %esc, %isc, %cid, %le, %sps);
my %intron_N_hit; # hash with the intron number of the intron hit in Sp2

open (IN, $input_file) || die "It cannot open File with Aln scores (all_PROT_EX_INT_aln_features_Sp1-Sp2.txt)\n";
### Format (multi-line: various formats in a single file
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
		for ($ps=1; $ps<=10; $ps++){ # why not directly: $isc{$id}=($line[7]*0.025);
		    if($line[7]==$ps) { 
			$isc{$id}=($ps*0.025); ##change here now is max score = 10
			# if there are gaps, longer distance is allowed in B1 (14/11/20). This is incorporated in the score
			my ($temp_int_N) = $line[5]=~/intron_(\d+)/;
			$intron_N_hit{$id} = $temp_int_N; # stores intron N in Sp2
		    }
		}
	    }
	}
    } 
    ### (3) not used
}
close IN;

# open all_PROT_EX_INT_aln_features_Sp1-Sp2.txt again
open (IN, $input_file) || die "It cannot open again File with Aln scores (Final_aln_scores_Sp1_Sp2.txt)\n";
#my $scorefile=$outf."/All_scores_by_exon.txt"; # old name
my $scorefile=$outf."/all_scored_EX_matches.txt";
open (OUTONE, ">$scorefile");
print OUTONE "CID\tFeature\tExon_prot_loc\tProt_query\tExon_number_query\tExon_coords_query\tProt_Subject\tExon_number_subject\tExon_coords_subject";
print OUTONE "\tScore_C1\tScore_I1\tScore_A\tScore_I2\tScore_C2\tTotal_exon_score\tSp_query\tSp_subject\n";
my ($e1, $i1, $e2, $i2, $c);
my ($nid, $pscore, $sC1, $sC2, $sA, $sI1, $sI2); 
my $TSC;
my (@g1,@g2);
my $idex;
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
	my $string; # it can be deprecated
	if ($line[2] eq "exon_1"){ ##N-terminal exon, we evaluate only I2 and C2;
	    $e2=$line[1]."#exon_2#".$line[4];
	    $i2=$line[1]."#intron_1#".$line[4];
	    if ($isc{$i2}) { if ($isc{$i2} eq "NO_INTRON"){ $sI2=0; } else { $sI2=$isc{$i2};}  } else { $sI2=-0.25 }
	    if ($esc{$e2}) { $sC2=($esc{$e2} * 0.15); } # prev = 0.16
	    if ($line[5] ne "NO_EXON_ALN") { 
		$sA=($line[7] * 0.20); # prev = 0.18
		$TSC=$sA+$sI2+$sC2+0.4; ##Total Score
		print OUTONE "$c\tExon\tN_terminal\t$line[1]\t$line[2]\t$line[3]\t";
		print OUTONE "$line[4]\t$line[5]\t$line[6]\tNA\tNA\t$sA\t$sI2\t$sC2\t$TSC\t$line[8]\t$line[9]\n";
		$string=$c."\tExon"."\t"."N_terminal\t".$line[1]."\t".$line[2]."\t".$line[3]."\t";
		$string.=$line[4]."\t".$line[5]."\t".$line[6]."\tNA\tNA\t".$sA."\t".$sI2."\t".$sC2."\t".$TSC."\t".$line[8]."\t".$line[9];
	    }
	    else { 
		$sA = -1; ##Nothing aligned => previously -0.1 (14/11/20)
		$TSC=$sA+$sI2+$sC2+0.4; ##Total Score  (to sum up to 1, the missing exon and intron score are added)
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
		$TSC=$sC1+$sI1+$sA+0.40; ##Total Score  (to sum up to 1, the missing exon and intron score are added); ##Total Score
		print OUTONE "$c\tExon\tC_terminal\t$line[1]\t$line[2]\t$line[3]\t";
		print OUTONE "$line[4]\t$line[5]\t$line[6]\t$sC1\t$sI1\t$sA\tNA\tNA\t$TSC\t$line[8]\t$line[9]\n";
		$string=$c."\t"."Exon"."\t"."C_terminal\t".$line[1]."\t".$line[2]."\t".$line[3]."\t";
		$string.=$line[4]."\t".$line[5]."\t".$line[6]."\t".$sC1."\t".$sI1."\t".$sA."\tNA\tNA\t".$TSC."\t".$line[8]."\t".$line[9];
	    }
	    else { 
		$sA = -1; ##Nothing aligned => previously -0.1 (14/11/20) 
		$TSC=$sC1+$sI1+$sA+0.4; # the 0.4 is for "perfect" I2 C2
		print OUTONE "$c\tExon\tC_terminal\t$line[1]\t$line[2]\t$line[3]\t";			
		print OUTONE "$line[4]\t$line[5]\tNA\t$sC1\t$sI1\t$sA\tNA\tNA\t$TSC\t$line[8]\t$line[9]\n";
		$string=$c."\t"."Exon"."\t"."C_terminal\t".$line[1]."\t".$line[2]."\t".$line[3]."\t";
		$string.=$line[4]."\t$line[5]\tNA\t".$sC1."\t".$sI1."\t".$sA."\tNA\tNA\t".$TSC."\t".$line[8]."\t".$line[9];
	    }
	}
	else { ##internal exon
	    @l1=split(/\_/,$line[2]);
	    $e1=$line[1]."#exon_".($l1[1]-1)."#".$line[4];
	    $i1=$line[1]."#intron_".($l1[1]-1)."#".$line[4]; # intron -1 in Sp1
	    $e2=$line[1]."#exon_".($l1[1]+1)."#".$line[4];   
	    $i2=$line[1]."#intron_".($l1[1])."#".$line[4];   # intron +1 in Sp1

	    if ($isc{$i1}) { if ($isc{$i1} eq "NO_INTRON"){ $sI1=0; } else { $sI1=$isc{$i1}; }  } else { $sI1=-0.25; } 
	    if ($esc{$e1}) { $sC1=($esc{$e1} * 0.15); } # prev = 0.16
	    if ($isc{$i2}) { # the check for consecutive introns is done ARBITRARILY for intron 2
		if ($isc{$i2} eq "NO_INTRON") { $sI2=0; } 
		else { 
		    $sI2=$isc{$i2}; 
		    # adds the check for Sp2 introns -1 and +1 being consecutive
		    if ($sI2>0 && $sI1>0) { # if both introns seem OK
			$sI2=-0.01 if abs($intron_N_hit{$i1}-$intron_N_hit{$i2}) != 1; # also if it's the same intron (not sure it can happen)
		    }
		} 
	    } else{ $sI2=-0.25; }
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
		$sA = -1; ##Nothing aligned => previously -0.1 (14/11/20)  
		$TSC=$sC1+$sI1+$sA+$sI2+$sC2;
		print OUTONE "$c\tExon\tInternal\t$line[1]\t$line[2]\t$line[3]\t";			
		print OUTONE "$line[4]\t$line[5]\tNA\t$sC1\t$sI1\t$sA\t$sI2\t$sC2\t$TSC\t$line[8]\t$line[9]\n";
		$string=$c."\t"."Exon"."\t"."Internal\t".$line[1]."\t".$line[2]."\t".$line[3]."\t";
		$string.=$line[4]."\t$line[5]\tNA\t".$sC1."\t".$sI1."\t".$sA."\t".$sI2."\t".$sC2."\t".$TSC."\t".$line[8]."\t".$line[9];			
	    }
	}
    }
}
close (OUTONE);
