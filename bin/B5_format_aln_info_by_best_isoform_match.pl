#!/usr/bin/env perl
use warnings;
use strict;


##SCRIPT FOR SCORING EXONS BTW PAIR OF SPECIES##
my $s1=$ARGV[0]; ##species 1 
my $s2=$ARGV[1]; ##species 2
my $i1=$ARGV[2]; ##aligned protein sp1 sp2 (Aligned_proteins.txt)
my $i2=$ARGV[3]; ##exon scores (Best_scores_pair_exons.txt)
my $i3=$ARGV[4]; ##intron scores (Score_all_introns.txt)
my $i4=$ARGV[5]; ##exint file sp1 (Sp1.exint)
my $i5=$ARGV[6]; ##exint file sp1 (Sp2.exint)
my $i6=$ARGV[7]; ##intron coords sp1 (Sp1_protein_ids_intron_pos_CDS.txt)
my $i7=$ARGV[8]; ##intron coords sp2 (Sp2_protein_ids_intron_pos_CDS.txt)
my $outf=$ARGV[9];

##protein ids exons species 1
open (IN,"$i1") || die "Missing aligned proteins file";  
# Format: 13	Protein	ENSDARP00000115151|ENSDARG00000056001	ENSP00000351113|ENSG00000075539	69.33	67.82	57.95	0.57	Dre	Hsa
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
open (EXONE, "$f1") || die "Missing exint file species 1";
#Format: >ENSP00000379387|ENSG00000185219  144.0 200.2 232.0 274.2 311.2
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
my @k;
my ($posint1, $posint2);
open (IN,"$i6")|| die "It cannot open intron positions for species 1 ($i6)\n"; #  Intron coords sp1 (Sp1_protein_ids_intron_pos_CDS.txt)
# Format: BL22196_evm0|BL22196	Sc0000025	-	Intron_1	366299-366745
while (<IN>){
    chomp($_);
    @line=split(/\t/,$_);
    $line[3]=~s/I/i/;
    if ($line[3] ne "Intronless" && $line[3] ne "intronless"){ # it's small letters, but just in case both
	if ($line[4]=~/\-/){
	    @k=split(/\-/, $line[4]); # && $line[4]=~/\-/);
	    $posint1=$k[0]+1; $posint2=$k[1]-1;
	    $intpos{$line[0]."_".$line[3]}=$line[1].":".$posint1."-".$posint2.":".$line[2];
	}
    }
}
open (IN,"$i7")|| die "It cannot open intron positions for species 2 ($i7)\n"; #  Intron coords sp1 (Sp2_protein_ids_intron_pos_CDS.txt) 
while (<IN>){
    chomp($_);
    @line=split(/\t/,$_);
    $line[3]=~s/I/i/;
    if ($line[3] ne "Intronless" && $line[3] ne "intronless"){# && $line[4]=~/\-/){
	if ($line[4]=~/\-/){
	    @k=split(/\-/, $line[4]);
	    $posint1=$k[0]+1; $posint2=$k[1]-1;
	    $intpos{$line[0]."_".$line[3]}=$line[1].":".$posint1."-".$posint2.":".$line[2];
	}
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
open (IN,"$i2")|| die "It cannot open exon intron scores file (Best_scores_pair_exons.txt)\n"; ###exon intron scores
while (<IN>){
    chomp($_);	
    @line=split(/\t+/,$_); ##change of format
    #Exon	ENSMUSP00000000028|ENSMUSG00000000028	exon_10	236-275	chr16	18795918-18796037	-	ENSP00000397434|ENSG00000093009	NO_EXON_ALN
    #13	ENSDARP00000044293|ENSDARG00000056001	exon_28	1152-1188	chr10	33997813-33997923	+	1	ENSP00000261575|ENSG00000073910	1-3	2.70	5.41	34	91.89	ENSP00000261575|ENSG00000073910	exon_1	1-54	chr13	32768274-32768434	+	Dre	Hsa
    if ($line[2]=~/exon/ && $line[7]==1){
	$coords{$line[1]."#".$line[2]}=$line[4].":".$line[5].":".$line[6];	 # some exons are missing for some prots
	$coords{$line[14]."#".$line[15]}=$line[17].":".$line[18].":".$line[19];  # some exons are missing for some prots
	$sps{$line[1]}=$line[20];
	$sps{$line[8]}=$line[21];
	$id=$line[1]."#".$line[2]."#".$line[8];		
#	if ($line[11]>0){ 
	if ($line[11]=~/\d/){ # if sim = 0 will skip it (i.e. in the case of some microexons between flies and mammals)
	    $line[11]=sprintf("%.2f",($line[11]/100));  
	    $exsc{$id}=$line[0]."\t".$line[1]."\t".$line[2]."\t".$line[4].":".$line[5].":".$line[6]."\t".$line[8].
		"\t".$line[15]."\t".$line[17].":".$line[18].":".$line[19]."\t".$line[11]."\t".$line[20]."\t".$line[21];
	}
	else {}
    }
}
open (IN,"$i3")|| die "It cannot open intron scores (Score_all_introns.txt)\n";
# Format: GF000006	ENSP00000268043|ENSG00000140451	intron_4	324	1	ENSMUSP00000049046|ENSMUSG00000041064	intron_4	333	1	0	z,z	25	Hs2	Mm2
while (<IN>){
    chomp($_);	
    @line=split(/\t+/,$_);
    if ($line[2]=~/intron/){
	$id=$line[1]."#".$line[2]."#".$line[5];
	$pint1=$line[1]."_".$line[2];
	if ($intpos{$line[5]."_".$line[6]}){ $pint2=$intpos{$line[5]."_".$line[6]}; }
	else { $pint2="NA";  }
	$insc{$id}=$line[0]."\t".$line[1]."\t".$line[2]."\t$intpos{$pint1}\t".$line[5]."\t".$line[6]."\t$pint2\t".$line[11]."\t".$line[12]."\t".$line[13];

	### Reciprocal hit not enforced
	#$id=$line[5]."#".$line[6]."#".$line[1];
	#$insc{$id}=$line[0]."\t".$line[5]."\t".$line[6]."\t".$line[7]."\t".$line[8]."\t".$line[1]."\t".$line[2]."\t".$line[3]."\t".$line[4]."\t".$line[9];
    }
}

### Some exons without aln miss basic info (rescued 22/08/20)
### Tries to gather the info from the same EXONS_DB as intron annotations
my ($root_sp1)=$i6=~/(.+?)_protein_ids_intron_pos_CDS.txt/;
my ($root_sp2)=$i7=~/(.+?)_protein_ids_intron_pos_CDS.txt/;
open (EX_ANNOT_1, "$root_sp1"."_protein_ids_exons_pos.txt") || die "It cannot open the exon positions internally for Sp1 ($root_sp1"."_protein_ids_exons_pos.txt)\n";
# Format: GB42164-PA|GB42164   exon_3   94-166   chr1  10316-10534  -
while (<EX_ANNOT_1>){
    chomp($_);
    my @t=split(/\t/);
    my $ex_id = "$t[0]#$t[1]"; # ProtID|GeneID#exon_2  
    if (!defined $coords{$ex_id}){ # it only fills gaps
	$coords{$ex_id}="$t[3]:$t[4]:$t[5]";
    }
    if (!defined $sps{$t[0]}){
	$sps{$t[0]}=$s1;
    }
}
close EX_ANNOT_1;

open (EX_ANNOT_2, "$root_sp2"."_protein_ids_exons_pos.txt") || die "It cannot open the exon positions internally for Sp2 ($root_sp2"."_protein_ids_exons_pos.txt)\n";
# Format: GB42164-PA|GB42164   exon_3   94-166   chr1  10316-10534  -
while (<EX_ANNOT_2>){
    chomp($_);
    my @t=split(/\t/);
    my $ex_id = "$t[0]#$t[1]"; # ProtID|GeneID#exon_2  
    if (!defined $coords{$ex_id}){ # it only fills gaps
	$coords{$ex_id}="$t[3]:$t[4]:$t[5]";
    }
    if (!defined $sps{$t[0]}){
	$sps{$t[0]}=$s2;
    }
}
close EX_ANNOT_2;

##### Starts the actual processing and putting together
my @keys;
@keys=sort(keys(%aln));
my $el;
my $prots;
my $nex;
my $id2;
my $cr;
my $cl;
my $ne1;
my $ne2;
#### Opens output file
open (OUT, ">$outf") || die "It cannot open the output file ($outf)\n";
print OUT "CID\tProt_query\tFeat_query\tCoords_query\tProt_subject\tFeat_subject\tCoords_subject\tScore_feature\tSp_query\tSp_subject\n";
foreach $el (@keys){
    @tmp=split(/\#/,$el);
    $prots=$aln{$el};
#    print "$prots\n";
    my $ex1=0; my $ex2=0;
    if ($exons{$tmp[0]}){ # $tmp[0] => ProtID|GeneID,  $exons{$tmp[0]} => Number of exons
	$cl=$pid1{$tmp[0]}; 

	### Puts together the exon scores
	if ($exons{$tmp[0]}==1){ # only one exon
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
	    for ($m=1; $m<$nex; $m++){ # loops through all exons in the gene
		$id=$tmp[0]."#exon_".$m."#".$tmp[1];
		$cr=$coords{$tmp[0]."#"."exon_".$m};
		if ($exsc{$id}){ # if it has score, all good.
		    print OUT "$exsc{$id}\n";
		}		
		else {  # if it does NOT have score (some missed info; added extra code above to get annotations)
		    print OUT "$cl\t$tmp[0]\texon_$m\t$cr\t$tmp[1]\tNO_EXON_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
		}

		### Now it gets the intron info
		$id2=$tmp[0]."#intron_".$m."#".$tmp[1];
		if ($insc{$id2}){ # if there is score, it's easy
		    print OUT "$insc{$id2}\n";
		}		
		else {
		    $ne1=$tmp[0]."#exon_".$m."#".$tmp[1]; 
		    $ne2=$tmp[0]."#exon_".($m+1)."#".$tmp[1];
		    $pint1=$tmp[0]."_intron_".$m;
		    
		    if (!$exsc{$ne1} && !$exsc{$ne2}){ 
			print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
		    }
		    elsif ($exsc{$ne1} && !$exsc{$ne2}){
			my @temp1=split(/\t/,$exsc{$ne1});
			if ($temp1[7]<0.2){
			    print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
			}
			else { # if the Aln is decent, it's assigned as no intron
			    print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_INTRON\tNA\t0\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
			}
		    }
		    elsif ($exsc{$ne1} && $exsc{$ne2}){
			my @temp1=split(/\t/,$exsc{$ne1}); # the info for exon_X
			my @temp2=split(/\t/,$exsc{$ne2}); # the info for exon_X+1
			if ($temp1[7]<0.2 && $temp2[7]<0.2){
			    print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
			}
			else {
			    print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_INTRON\tNA\t0\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
			}
		    }
		    elsif (!$exsc{$ne1} && $exsc{$ne2}){
			my @temp2=split(/\t/,$exsc{$ne2}); # the info for exon_X+1  
			if ($temp2[7]<0.2){
			    print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_ALN\tNA\tNA\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
			}
			else {
			    print OUT "$cl\t$tmp[0]\tintron_$m\t$intpos{$pint1}\t$tmp[1]\tNO_INTRON\tNA\t0\t$sps{$tmp[0]}\t$sps{$tmp[1]}\n";
			}
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
close OUT;
