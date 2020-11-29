#!/usr/bin/env perl
use warnings;
use strict;

my ($i,$input_file,$sps,$out,$out2,$m,$int,$id,$C,$A);
for ($i=0; $i<=$#ARGV; $i++){
    if ($ARGV[$i] eq "-b"){ ## all_score_EX
	$i++;
	$input_file = $ARGV[$i]; 		
    }
    if($ARGV[$i] eq "-sps"){  ##species separate by comma
	$i++;
	$sps=$ARGV[$i]; 		
    }	
    if($ARGV[$i] eq "-out"){ ##folder with reference proteins (def: filtered_best_scored_EX_matches_by_targetgene.tab)
	$i++;
	$out=$ARGV[$i]; 		
    }	
    if($ARGV[$i] eq "-out2"){ ##folder with reference proteins (def: best_scored_EX_matches_by_targetgene.txt)
	$i++;
	$out2=$ARGV[$i]; 		
    }	
    if($ARGV[$i] eq "-max_size"){ ##maximum % in size difference between the two exons (for exons smaller than 30, this is fixed to 30%)
	$i++;
	$m=$ARGV[$i]; 		
    }
    if($ARGV[$i] eq "-int"){  ##number of conserved introns either 1 or 2 (def:2) 	
	$i++;
	$int=$ARGV[$i]; 		
    }
    if($ARGV[$i] eq "-id"){  ##% of protein similarity of the exon subject and upstream and downstream exon (def:20%)	
	$i++;
	$id=$ARGV[$i]; 		
    }
    if($ARGV[$i] eq "-h"){	
	print "For running the program: perl filt_exons_by_score.pl -b <all_score_EX> -sps <species list separated by comma> -out outfile -out2 best_scored
-max_size <max % dif size between scored exons> 
-id <min %sim between exons pairs and upstream and downstream exons>
-int <min number of neighboring introns with conserved positions (1 or 2)>
 
\n";
    }
}

if (!$m){ $m=0.65}
if ($int !~ /\d/){ die "ERROR: Intron number must be 0, 1 or 2\n";}
if (!$id){ $id = 0.2; $A = 0.20*0.2; $C = 0.15*0.2; } else { $A=0.20*$id; $C=0.15*$id; } # prev weights: A = 0.18 and C = 0.16
if (!$out){ $out="filtered_best_scored_EX_matches_by_targetgene.tab"; }
if (!$out2){ $out2="best_scored_EX_matches_by_targetgene.txt"; }
my @sps=sort(split(/\,/,$sps));
my %sid;
my $k;
for ($k=0; $k<scalar(@sps); $k++){ $sid{$sps[$k]}=$k+1;	} #assigning id to each species
my (%pre_exon_score, %best_hit_pre);               # store info pre-filtering
my (%post_exon_score, %best_hit_post, %best_pair); # store info post-filtering

### Opens and parses ALL scored exons
open (INONE, $input_file) || die "It cannot open input file ($input_file)\n";
# Format of input:
# 4282	Exon	Internal	ENSDARP00000111316|ENSDARG00000086505	exon_15	chr11:13246628-13246793:-	ENSP00000262811|ENSG00000099308	exon_16	chr19:18245627-18245792:+	0.1183	0.25	0.1372	0.25	0.1144	0.8699	Dre	Hsa
my $header = <INONE>; # keeps header for output
while (<INONE>){
    chomp($_); 
    my $rs=0; # whether the intron positions are conserved (def = NO)
    my @line=split(/\t/,$_);
    
    ### First it gets the best hit by target get by score alone (no filter). Equivalent to former C1 script.
    my @t_g1=split(/\|/,$line[3]); my @t_g2=split(/\|/,$line[6]);
    my $t_id1=$t_g1[1]."\t".$line[5];
    if (!defined $pre_exon_score{$t_id1}{$t_g2[1]}){
	$best_hit_pre{$t_id1}{$t_g2[1]}=$_;
	$pre_exon_score{$t_id1}{$t_g2[1]}=$line[14];
    } 
    else {
	if ($line[14] > $pre_exon_score{$t_id1}{$t_g2[1]}){
	    $best_hit_pre{$t_id1}{$t_g2[1]}=$_;
	    $pre_exon_score{$t_id1}{$t_g2[1]}=$line[14];
	}
    }
        
    ### Does the filtering and selection of best match POST-filtering
    if ($line[7] ne "NO_EXON_ALN"){
	# takes the coordinate and uses the two coordiantes to calculate the length of each exon
	my $ts1=$line[5]; $ts1=~s/\-/\:/;
	my @ts1=split(/\:/,$ts1); 
	my $sex1=$ts1[2]-$ts1[1]+1;
	my $ts2=$line[8]; $ts2=~s/\-/\:/;
	my @ts2=split(/\:/,$ts2);
	my $sex2=$ts2[2]-$ts2[1]+1;
	my ($max,$lim);
	# gets the maximum length (the largest of the two exons)
	if ($sex1>$sex2){ 
	    $max=$sex1;
	}else { $max=$sex2; }	    
	# setting a limit of the sizes of the two scored exons
	if ($max>30) { $lim=int($max*$m);}
	else { $lim=int($max*0.3);} # if it's a MIC of 18 and 6 => it's OK (5.4)
	
	# sets IDs, etc (before, it forced reciprocal fake mappings)
	my @g1=split(/\|/,$line[3]); my @g2=split(/\|/,$line[6]);		
	my $id1=$g1[1]."\t".$line[5];
	my $id2=$g2[1]."\t".$line[8];

	if ($line[2] eq "Internal") {
	    # tests if $int (1 or 2) intron positions are conserved
	    # add the explicit check for 1 or 0 introns.
	    if ($int==2){
		if ($line[10] > 0 && $line[12]>0){ $rs=1; }
	    } elsif ($int==1) {
		if ($line[10] > 0 || $line[12]>0){ $rs=1; }				
	    } elsif ($int==0) { 
		$rs=1; 
	    } else { 
		die "ERROR: The number of introns needs to be 0, 1 or 2\n";
	    }

	    # IF INTRON POSITIONS ARE CONSERVED:
	    if ($rs==1){
		# tests sequence conservation of the three exons
		# added a loophole by which tiny MICs get through without similarity
		if (($line[9]>=$C && $line[13]>=$C && $line[11]>=$A) || ($line[9]>=$C && $line[13]>=$C && $sex1<=9 && $sex2<=9 && $line[10]>0 && $line[12]>0)){ 
		    # filtering by sizes of the two exons that are being scored
		    if ($sex1 >= $lim && $sex2 >= $lim){ 
			### Get the best match post-filter and best pair
			if (!defined $post_exon_score{$id1}{$g2[1]}){
			    $best_hit_post{$id1}{$g2[1]}=$_;
			    $best_pair{$id1}{$g2[1]}="$id1\t$id2\t$line[15]\t$line[16]";
			    $post_exon_score{$id1}{$g2[1]}=$line[14];
			} 
			else {
			    if ($line[14] > $post_exon_score{$id1}{$g2[1]}){
				$best_hit_post{$id1}{$g2[1]}=$_;
				$best_pair{$id1}{$g2[1]}="$id1\t$id2\t$line[15]\t$line[16]";
				$post_exon_score{$id1}{$g2[1]}=$line[14];				
			    }
			}
		    }
		}
	    }
	}
	elsif ($line[2] eq "C_terminal"){
	    if (($line[10]>0 || $int == 0) && $line[9]>=$C && $line[11]>=$A){
		if($sex1>=$lim && $sex2>=$lim){ ##filtering by size diference 
		    ### Get the best match post-filter and best pair
		    if (!defined $post_exon_score{$id1}{$g2[1]}){
			$best_hit_post{$id1}{$g2[1]}=$_;
			$best_pair{$id1}{$g2[1]}="$id1\t$id2\t$line[15]\t$line[16]";
			$post_exon_score{$id1}{$g2[1]}=$line[14];
		    } 
		    else {
			if ($line[14] > $post_exon_score{$id1}{$g2[1]}){
			    $best_hit_post{$id1}{$g2[1]}=$_;
			    $best_pair{$id1}{$g2[1]}="$id1\t$id2\t$line[15]\t$line[16]";
			    $post_exon_score{$id1}{$g2[1]}=$line[14];				
			}
		    }		    
		}
	    }
	}
	elsif ($line[2] eq "N_terminal"){
	    if (($line[12] > 0 || $int == 0) && $line[11]>=$A && $line[13]>=$C){
		if ($sex1>=$lim && $sex2>=$lim){ ##filtering by size difference
		    ### Get the best match post-filter and best pair
		    if (!defined $post_exon_score{$id1}{$g2[1]}){
			$best_hit_post{$id1}{$g2[1]}=$_;
			$best_pair{$id1}{$g2[1]}="$id1\t$id2\t$line[15]\t$line[16]";
			$post_exon_score{$id1}{$g2[1]}=$line[14];
		    } 
		    else {
			if ($line[14] > $post_exon_score{$id1}{$g2[1]}){
			    $best_hit_post{$id1}{$g2[1]}=$_;
			    $best_pair{$id1}{$g2[1]}="$id1\t$id2\t$line[15]\t$line[16]";
			    $post_exon_score{$id1}{$g2[1]}=$line[14];				
			}
		    }
		}
	    }
	}
    } 
}
### It opens the output here, in case it crashes before
open (OUT, ">$out") || die "It cannot open the output file with filtered exons ($out)\n";
print OUT "GeneID_sp1\tExon_coords_sp1\tGeneID_sp2\tExon_coords_sp2\tSp1\tSp2\n";
foreach my $ex_sp1 (sort keys %best_pair){
    foreach my $g_sp2 (sort keys %{$best_pair{$ex_sp1}}){
	print OUT "$best_pair{$ex_sp1}{$g_sp2}\n";
    }
}
close OUT;

open (OUT2, ">$out2") || die "It cannot open the second output file ($out2)\n";
print OUT2 $header;
foreach my $ex_sp1 (sort keys %best_hit_pre){ # best_hit_pre should have hit for ALL exons (even without ALN)
    foreach my $g_sp2 (sort keys %{$best_hit_pre{$ex_sp1}}){
	if (defined $best_hit_post{$ex_sp1}{$g_sp2}){
	    print OUT2 "$best_hit_post{$ex_sp1}{$g_sp2}\n";
	}
	else {
	    print OUT2 "$best_hit_pre{$ex_sp1}{$g_sp2}\n";
	}
    }
}
close OUT2;
