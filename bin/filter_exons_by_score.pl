#!/usr/bin/env perl
#use warnings;
use strict;
my ($i,$i1,$sps,$out,$m,$int,$id,$C,$A);
for ($i=0; $i<=$#ARGV; $i++){
	if($ARGV[$i] eq "-b"){ ##best hit exons
		$i++;
		$i1=$ARGV[$i]; 		
	}
	if($ARGV[$i] eq "-sps"){  ##species separate by comma
		$i++;
		$sps=$ARGV[$i]; 		
	}	
	if($ARGV[$i] eq "-out"){ ##folder with reference proteins (def: Best_score_exon_hits_filtered_<parameters>.txt)	
		$i++;
		$out=$ARGV[$i]; 		
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
		print "For running the program: env perl filt_exons_by_score.pl -b <best_hit_exons> -s <species list separated by comma> -out outfile -max_size <max % dif size between scored exons> -int <introns to be scored>
-id <min %sim between exons pairs and upstream and downstream exons>\n";
	}
}

if (!$m){ $m=0.65}
if (!$int){ $int=2; }
if (!$id){ $id=0.2; $A=0.18*0.2; $C=0.16*0.2; } else { $A=0.18*$id; $C=0.16*$id; }
if (!$out){ $out="Best_score_exon_hits_filtered_$m-$int-$id".".tab"; }
#print "$i1\t$int\t$m\t$A\t$C\n";
my (@line);
my $c=0;
my @sps=split(/\,/,$sps);
@sps=sort(@sps);
my %sid;
my $k;
for ($k=0; $k<scalar(@sps); $k++){ $sid{$sps[$k]}=$k+1;	} #assigning id to each species
my ($s1, $s2, $id1, $id2, $l, $sex1, $sex2, $max, $lim, $ts1, $ts2, $rs);
my (@g1, @g2, @ts1, @ts2);
my %hits;
my %sp;
my (%sc,%pair,%ex, %pr);
#4282	Exon	Internal	ENSDARP00000111316|ENSDARG00000086505	exon_15	chr11:13246628-13246793:-	ENSP00000262811|ENSG00000099308	exon_16	chr19:18245627-18245792:+	0.1183	0.25	0.1372	0.25	0.1144	0.8699	Dre	Hsa
open (OUT, ">$out");
open (INONE,"$i1");
while (<INONE>){
	$c++;
	chomp($_); 
	$rs=0;
	if ($c==1){ } ##printing header
	@line=split(/\t/,$_);
	if ($line[7] ne "NO_EXON_ALN"){
	    $ts1=$line[5]; $ts1=~s/\-/\:/;
	    @ts1=split(/\:/,$ts1); 
	    $sex1=$ts1[2]-$ts1[1]+1;
	    $ts2=$line[8]; $ts2=~s/\-/\:/;
	    @ts2=split(/\:/,$ts2);
	    $sex2=$ts2[2]-$ts2[1]+1;
	    if ($sex1>$sex2){ 
		$max=$sex1;
	    }else { $max=$sex2; }
	    if ($max>30) { $lim=int($max*$m); }
	    else {  $lim=int($max*0.3);  } ##setting a limit of the sizes of the two scored exons
		@g1=split(/\|/,$line[3]); @g2=split(/\|/,$line[6]);		
		if ($sid{$line[15]}<$sid{$line[16]}) { 
			$id1=$g1[1]."\t".$line[5];
			$id2=$g2[1]."\t".$line[8];
			$sp{$id1."\t".$id2}=$line[15]."\t".$line[16];
			$sp{$id2."\t".$id1}=$line[16]."\t".$line[15];	
		} else { 
			$id2=$g1[1]."\t".$line[5];  
			$id1=$g2[1]."\t".$line[8]; 
			$sp{$id1."\t".$id2}=$line[16]."\t".$line[15];
			$sp{$id2."\t".$id1}=$line[15]."\t".$line[16];			
		}	   			
		if ($line[2] eq "Internal") {
		    if ($int==2){
					if ($line[10] > 0 && $line[12]>0){ $rs=1; }
				}else {
					if ($line[10] > 0 || $line[12]>0){ $rs=1; }				
				}
		    if ($rs==1){
    				if ($line[9]>=$C && $line[13]>=$C && $line[11]>=$A){
				       # print  "$_\n"; ##only if both intron positions are conserved					
			    		if ($sex1>=$lim && $sex2>=$lim){ ##adding a filtering by sizes of the two exons that are being scored
					    $hits{$id1."\t".$id2}=1;
							$sc{$id1."\t".$id2}+=$line[14]; 
							if (!$ex{$id1}){ $ex{$id1}=$id2; 
								$pair{$id1."\t".$id2}=1;				
							}
							else { 
								if (!$pair{$id1."\t".$id2}) {  
									$ex{$id1}.=",".$id2;  
									$pair{$id1."\t".$id2}=1;
								}
							} 
						}
			    	}
				}
		}
		elsif ($line[2] eq "C_terminal"){
			if ($line[10]>0 && $line[9]>=$C && $line[11]>=$A){
			    if($sex1>=$lim && $sex2>=$lim){ ##filtering by size diference 
				$hits{$id1."\t".$id2}=1;
				$sc{$id1."\t".$id2}+=$line[14]; 
				if (!$ex{$id1}){ $ex{$id1}=$id2; 
					$pair{$id1."\t".$id2}=1;				
				}
				else { 
					if (!$pair{$id1."\t".$id2}) {  
						$ex{$id1}.=",".$id2;  
						$pair{$id1."\t".$id2}=1;
					} 
				}
			    }
		       }
		}
		elsif ($line[2] eq "N_terminal"){
			if ($line[12] > 0 && $line[11]>=$A && $line[13]>=$C){
			    if ($sex1>=$lim && $sex2>=$lim){ ##filtering by size difference
				$hits{$id1."\t".$id2}=1;
				$sc{$id1."\t".$id2}+=$line[14]; 
				if (!$ex{$id1}){ $ex{$id1}=$id2; 
					$pair{$id1."\t".$id2}=1;				
				}
				else { 
					if (!$pair{$id1."\t".$id2}) {  
						$ex{$id1}.=",".$id2;  
						$pair{$id1."\t".$id2}=1;
					} 
				} 
			    }
			}
		}
	} 
}
print OUT "GeneID_sp1\tExon_coords_sp1\tGeneID_sp2\tExon_coords_sp2\tSp1\tSp2\n";
my @keys=keys(%ex);
@keys=sort(@keys);
my ($el, $nid, $score, $pr);
my (@e, @g); 
foreach $el (@keys){	
	if ($ex{$el}=~/\,/){
		@e=split(/\,/,$ex{$el});
		$score=0;
		for ($l=0; $l<scalar(@e);$l++){
			@g=split(/\t/,$e[$l]);
			#print "1)$el\t$g[0]\n";
			if ($sc{$el."\t".$e[$l]}> $sc{$el."\t".$g[0]}){
				$pr{$el."\t".$g[0]}=$el."\t".$e[$l];
				$score=$sc{$el."\t".$e[$l]};				
			}
		}
	}
	else {
		@g=split(/\t/,$ex{$el});
		$pr{$el."\t".$g[0]}=$el."\t".$ex{$el};
	}
}

@keys=keys(%pr);
@keys=sort(@keys);
foreach $el (@keys){
	$nid=$pr{$el};
	print OUT "$pr{$el}\t$sp{$nid}\n";
}




