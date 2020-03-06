#!/usr/bin/perl
#use warnings;
#use strict;
use Getopt::Long;

### AIMS OF SCRIPT
#Define a new set of clusters.

### Criteria to chose the best hit
#* Filter by expression and similarity patterns (sim_exp [ON|OFF]; def OFF, 
	#* if ON then two criteria are used: minimum expression values & max similarity and it requires a file with expression patterns and TS information
	#* Minimum expression values: to be chosen by the user min_glob (0 to 100%) or min_exp (cRPKM value) [def: min_glob 50 or min_exp 5]
	#* Highest similarity: Filter threshold of % of max similarity (0 to 100%, def: max_sim 75 )  


#### DECLARATION AND GENERATION OF VARIABLES
my $cl_file; ##Original cluster file
my $sps; #Species separated by comma y 3 letter code, to perform the reclustering
my $of_path=""; #orthofinder path 
my $verbose; # prints steps
my $help;
my $out;

### Get options:
Getopt::Long::Configure("no_auto_abbrev");
GetOptions(  "cl_file=s" => \$cl_file,	     
	     "sps=s" => \$sps,
	     "outfile=s" => \$out,
	     "of_path=s" => \$of_path,
	     "verbose" => \$verbose,
	     "help" => \$help
    );


if (!defined ($cl_file) || !defined($sps) || !defined($of_path)  || defined ($help)){
    die "
Usage: perl Recluster_genes_species.pl --cl_file Gene_cluster_file --sps Species --outfile name_outfile --of_path Orthofinder_path Similarity_file [options]

OPTIONS
    --cl_file	Gene cluster file
    --sps	Species separated by comma to perform the reclustering in 3 letter code; e.g. Danio rerio (Dre)
    --outfile	Outfile name
    --of_path	Path to orthologues folder by pair of species orthologues of Orthofinder output		
    --verbose   Prints commands in screen [def OFF]   
";
}

@sp=split(/\,/,$sps);
for ($m=0; $m<scalar(@sp); $m++){
	$sid{$sp[$m]}=1;
}

#GF000001	Mmu	ENSMUSG00000050459	NA	OG0000001	O=0,R=0
#1) Checking how many clusters have the two species and the number of genes in such clusters by species
open (CLS, $cl_file);
while (<CLS>){
	chomp($_);
	@l=split(/\t/,$_);
	if ($_=~/GF/){
		if ($sid{$l[1]}){
			$cid{$l[0]}{$l[1]}++;
			$csp{$l[0]}=1;
			$cgid{$l[2]}=$l[0];
			$sp{$l[2]}=$l[1];
			$nid{$l[2]}=$l[3];
			#print "0)$l[2]\t$l[0]\n";
		} 
	}	
}

my @keys=keys(%csp);
$r=scalar(@keys);
foreach $el (@keys){
	for ($m=0; $m<scalar(@sp); $m++){
		for ($n=0; $n<scalar(@sp); $n++){
			if ($m != $n && !$ch{$sp[$m]}{$sp[$n]} && !$ch{$sp[$n]}{$sp[$m]} ){
				$ch{$sp[$m]}{$sp[$n]}=1;
				if ($cid{$el}{$sp[$m]} && $cid{$el}{$sp[$n]}){
					$tc{$sp[$m]}{$sp[$n]}+=$cid{$el}{$sp[$m]};
					$tc{$sp[$n]}{$sp[$m]}+=$cid{$el}{$sp[$n]};
					$tc2sp{$sp[$m]}{$sp[$n]}++;
					$rcl{$el}{$sp[$m]}{$sp[$m]}=1;
					if ($cid{$el}{$sp[$m]}==1 || $cid{$el}{$sp[$n]}==1){
							$onehit{$el}=1;
					}
				}
			}	
		}
	}
	%ch=();
}
%ch=();
$ost="Stats_".$out;
open (STATS, ">$ost");
for ($m=0; $m<scalar(@sp); $m++){
	for ($n=0; $n<scalar(@sp); $n++){
		if ($m != $n && !$ch{$sp[$m]}{$sp[$n]} && !$ch{$sp[$n]}{$sp[$m]} ){
			$ch{$sp[$m]}{$sp[$n]}=1;
			print STATS "Total shared clusters $sp[$m],$sp[$n]:\t$tc2sp{$sp[$m]}{$sp[$n]}\n";
			print STATS "Max genes for $sp[$m] with $sp[$n]:\t$tc{$sp[$m]}{$sp[$n]}\n";
			print STATS"Max genes for $sp[$n] with $sp[$m]:\t$tc{$sp[$n]}{$sp[$m]}\n";
		}	
	}	print STATS "\n";
}

open (CLS, $cl_file);
while (<CLS>){
	chomp($_);
	@l=split(/\t/,$_);
	if ($_=~/GF/){
		if ($sid{$l[1]}){
			if ($onehit{$l[0]}){
				$nw=$l[0].".01";
				$ncid{$l[2]}=$nw;
				if (!$ncl{$nw}){
					$ncl{$nw}=$l[2];
					$uhit{$l[2]}=$nw;
				}else {  $ncl{$nw}.=",".$l[2]; }			
			}		
		} 
	
	}	
}


#2) Opening orthofinder files
$files=`ls $of_path/Orthologues_*`;
@fs=split(/\n/,$files);
for ($r=0; $r<scalar(@fs); $r++){
	@tn=split(/\_/,$fs[$r]);
	#print "1)$fs[$r]\n";
	$n1=substr($tn[0],0,1); 
	$n2=substr($tn[1],0,2);
	$nn=$n1.$n2;
	#print "#$nn\n";
	if ($sid{$nn}){
		$ns=$tn[0]."_".$tn[1];
		$f2=`ls $of_path/Orthologues_$ns/*.csv`;
		#print "2)$f2\n";
		@fs2=split(/\n/,$f2);
		for ($k=0; $k<scalar(@fs2); $k++){
			@tn2=split(/\_+/,$fs2[$k]);
			$m1=substr($tn2[scalar(@tn2)-2],0,1); 
			$m2=substr($tn2[scalar(@tn2)-1],0,2);
			$mm=$m1.$m2;
			if ($sid{$mm}  && !$file{$nn.",".$mm} && !$file{$mm.",".$nn}){
				$s1=$nn;
				$s2=$mm;
				$file{$s1.",".$s2}=1;
				#print "$s1 $s2\n";
				#OG0000000	ENSP00000366359|ENSG00000243729	ENSMUSP00000094933|ENSMUSG00000090675, ENSMUSP00000132208|ENSMUSG00000090894
				open (FILE,$fs2[$k]);
				while (<FILE>){
					chomp($_);
					@l=split(/\t/,$_);
					@p1=(); @p2=();
					$l[1]=~s/\s+//g;
					$l[2]=~s/\s+//g;
					#print "#5) $l[1]\t=$l[2]\n"; 
					if ($l[1]=~/\,/){	
						@p1=split(/\,/,$l[1]);
					}else { push (@p1,$l[1]); } 
					if ($l[2]=~/\,/){						
						@p2=split(/\,/,$l[2]);
					}else { push (@p2,$l[2]); } 
					$f=0;
					$nw="";
					for ($x=0; $x<scalar(@p1); $x++){
						@tmp1=split(/\|/,$p1[$x]);
						$cid1=$cgid{$tmp1[1]};
						#print "#1)$p1[$x]#\t*1a)$p2[$y]*###\n";
						if (!$uhit{$tmp1[1]}){ ##
							if ($cgid{$tmp1[1]}){
								if ($ncid{$tmp1[1]}){
									$nw=$ncid{$tmp1[1]};
								}
								elsif (!$nw) {
									if(!$cs1{$cid1}){
										$cs1{$cid1}=1;
										$nw=$cid1.".".sprintf("%02d",$cs1{$cid1});
										$ncid{$tmp[1]}=$nw;
										$f=1;
										$ncl{$nw}=$tmp1[1];
									}
									elsif ($cs1{$cid1} && $f==0){
										$cs1{$cid1}++;
										$f++;
										$nw=$cid1.".".sprintf("%02d",$cs1{$cid1});
										$ncid{$tmp[1]}=$nw;
										$f=1;
										$ncl{$nw}=$tmp1[1];
									}
									elsif ($cs1{$cid1} && $f==1){
										$nw=$cid1.".".sprintf("%02d",$cs1{$cid1});
										$ncid{$tmp[1]}=$nw;
										$f=1;
										$ncl{$nw}.=",".$tmp1[1];
									}
								}
							}				
							for ($y=0; $y<scalar(@p2); $y++){						
								@tmp2=split(/\|/,$p2[$y]);
								if ($ncid{$tmp2[1]}){ $nw=$ncid{$tmp2[1]}; }
								#print "4)$p2[$y]*$tmp2[1]\n";
								#print "2)$tmp1[1]##\t**$tmp2[1]**\t=$cgid{$tmp1[1]}=\t+$cgid{$tmp2[1]}+###\t$nw\t$cid1\n";
								if ($cgid{$tmp2[1]}){ 
									if ($cgid{$tmp2[1]} eq $cid1){
										#print "2)$tmp1[1]##\t**$tmp2[1]**\t=$cgid{$tmp1[1]}=\t+$cgid{$tmp2[1]}+###\t$nw\t$cid1\n";
										#$ncid{$tmp2[1]}=$nw;
										if (!$ncid{$tmp1[1]}){ $ncid{$tmp1[1]}=$nw; $ncl{$nw}.=",".$tmp1[1]; }
										if (!$ncid{$tmp2[1]}){ $ncid{$tmp2[1]}=$nw; $ncl{$nw}.=",".$tmp2[1]; }
										#print "3) $nw\t$ncl{$nw}\n";
									} 	
								}
							} ##
						}				
					}
					
				}

			}	
		}
	}
}##
my @keys=keys(%ncl);
@keys=sort(@keys);
foreach $el (@keys){
	@mb=split(/\,/,$ncl{$el});
	@mb=sort(@mb);
	for ($j=0; $j<scalar(@mb); $j++){
		$spid=$sp{$mb[$j]};
		if (!$spcl{$el}{$spid}){
			$spcl{$el}{$spid}=1;
			$nsps{$el}++;
		}
	}
}

open (OUT, ">$out");
my @keys=keys(%ncl);
@keys=sort(@keys);
foreach $el (@keys){
	if ($nsps{$el}>1){
		@mb=split(/\,/,$ncl{$el});
		@mb=sort(@mb);
		for ($j=0; $j<scalar(@mb); $j++){
			$spid=$sp{$mb[$j]};
			if (!$ch{$el}{$spid}{$mb[$j]}){
				$gid=$mb[$j];
				print OUT "$el\t$sp{$mb[$j]}\t$gid\t$nid{$gid}\n";		
				$rec{$spid}++;
				$ch{$el}{$spid}{$mb[$j]}=1;
			}
		}
		print OUT "\n";
	}
}

print STATS "Reclustered genes:\n";
my @k2=keys(%rec);
foreach $s (@k2){
	print STATS "$s\t$rec{$s}\n";
}
