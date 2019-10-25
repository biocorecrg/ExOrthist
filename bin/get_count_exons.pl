#!/usr/bin/env perl
#use warnings;
#use strict;

my $i1=$ARGV[0]; #Best_score_exon_hits_filtered.txt
my $out=$ARGV[1]; #Lift over exons
my $i2=$ARGV[2];
my $evs=$ARGV[3];##folder with event-gene ids from vastdb `ls /users/mirimia/ymarquez/ANALYSIS/ORTOLOGUES/EXON_ORTHOLOGY_V2/EXONS_DB/VAST_TOOLS_INFO/*Event-Gene.IDs.txt`;
my $refs=$ARGV[4];##folder with REFERENCE tables of vastdb`ls /users/mirimia/ymarquez/ANALYSIS/ORTOLOGUES/EXON_ORTHOLOGY_V2/EXONS_DB/VAST_TOOLS_INFO/REFERENCE-ALL_ANNOT*`;
 
my (@l);
my %hits;
my %sp;
my (%sc,%pair,%ex, %pr);
#BL00000	Sc0000095:745259-745367:-	SMAR011687	JH432088:105973-106084:-	Bla	Sma
open (OUT, ">$out");
open (INONE,"$i1");
while (<INONE>){
	chomp($_); 
	@l=split(/\t/,$_);
	$e1=$l[0]."\t".$l[1]."\t".$l[4];
	$e2=$l[2]."\t".$l[3]."\t".$l[5];	
	$ex{$e1}++;
	$ex{$e2}++;
	$gene{$l[0]}=1;

}
if ($i2){
open (INONE,"$i2");
while (<INONE>){
	chomp($_); 
	@l=split(/\t/,$_);
	$e1=$l[0]."\t".$l[1]."\t".$l[4];
	$e2=$l[2]."\t".$l[3]."\t".$l[5];	
	$ex{$e1}++;
	$ex{$e2}++;
	$gene{$l[0]}=1;
}
}
if ($evs && $refs){
	#HsaEX0042157	ENSG00000229653
	my @fs=split(/\n/,$evs);
	for ($n=0; $n<scalar(@fs); $n++){
		open (INONE,"$fs[$n]");
		while (<INONE>){ 
			chomp($_); 
			@l=split(/\t/,$_);
			if ($l[0]=~/EX/){
				$gid{$l[0]}=$l[1];	
			}
		}
	} 

	#KLHL17	HsaEX6027622	chr1:961826-962047	222	chr1:960921,961826-962047,962355	A_C3	chr1:960921,961826-962047,962355:+	222	chr1:960821-960921	chr1:961826-962047	chr1:962355-962471
	my @fs=split(/\n/,$refs);
	for ($n=0; $n<scalar(@fs); $n++){
		open (INONE,"$fs[$n]");
		@nm=split(/\-/,$fs[$n]);
		$sp=substr($nm[2],0,3); 
		while (<INONE>){ 
			chomp($_); 
			@l=split(/\t/,$_);
			if ($l[1]=~/EX/){
				$g=$gid{$l[1]};
				if ($gene{$g}){
					@k=split(/\:/,$l[6]);
					$coords=$l[9].":".$k[2];
					$e1=$g."\t".$coords."\t".$sp;
					$ex{$e1}+=1000; 
				}	
			}
		}
	}
} 

my @keys=keys(%ex);
@keys=sort (@keys);
foreach $el (@keys){
	print OUT "$el\t$ex{$el}\n";
}









