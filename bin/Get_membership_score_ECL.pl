#!/usr/bin/perl
use warnings;
use strict;

my $exons_db_folder=$ARGV[0]; #EXONS_DB folder
my $infile1=$ARGV[1]; ##Best exon hits filtered
my $infile2=$ARGV[2]; ##Exon clusters
my $outfile=$ARGV[3]; ##Exon clusters with membership score
my ($id1, $as, $ts, $l, $id2, $n);

$ef=`ls $exons_db_folder/*/*_prot_exons_overlap.txt`;
@f=split(/\n/,$ef);
for ($k=0; $k<scalar(@f); $k++){
    open (GCL, $f[$k]);
    #OV_EX_Hs2_23	ENSG00000000419	50948629-50948662
    while (<GCL>){
	chomp($_);
	@l=split(/\t/,$_);
	$oid{$l[1]."\t".$l[2]}=$l[0];
    }	
}

#ENSG00000203814	chr1:149783236-149783260:-	ENSMUSG00000069308	chr13:21789024-21789060:+	Hs2	Mm2	1.7204
open (BH, $infile1) || die "It cannot open the best exon hits filtered ($infile1)\n";
while (<BH>){
    chomp($_);
    @l=split(/\t/,$_);
    @t1=split(/\:/,$l[1]);
    @t2=split(/\:/,$l[3]);
    $r1=$oid{$l[0]."\t".$t1[1]};
    $r2=$oid{$l[2]."\t".$t2[1]};
    $id1=$r1."\t".$r2;
    $id2=$r2."\t".$r1;
    $pair{$id1}=1;
    $pair{$id2}=1;
}
close BH;

#GF0000007.001	ENSMUSG00000110439	chr4:60821320-60821418:-	Mm2
#GF0000007.001	FBgn0264775	chr3R:15509753-15509899:-	Dme
open (ECL, $infile2) || die "It cannot open the exon clusters file ($infile2)\n";
while (<ECL>){
    chomp($_);
    @l=split(/\t/,$_);
    @t1=split(/\:/,$l[2]);
    $idov=$oid{$l[1]."\t".$t1[1]};
    $eclsp{$l[0]}{$l[3]}++;	
    $sp{$idov}=$l[3];
    if (!$excl{$l[0]}){
	$excl{$l[0]}=$idov;
    }
    else {
	$excl{$l[0]}.=",".$idov;
    }  
}
close ECL;

## prepares score
my @keys=keys(%excl);
foreach $cl (@keys){
    @ne=split(/\,/,$excl{$cl});
    $nex=scalar(@ne);
    $exs{$cl}=$nex;
    for ($n=0; $n<scalar(@ne)-1; $n++){
	for ($m=0; $m<scalar(@ne); $m++){
	    if (!$check{$ne[$m]}{$ne[$n]} && !$check{$ne[$n]}{$ne[$m]}){
		$check{$ne[$m]}{$ne[$n]}=1;
		$check{$ne[$n]}{$ne[$m]}=1;
		if ($sp{$ne[$n]} ne $sp{$ne[$m]}){
		    if ($pair{$ne[$n]."\t".$ne[$m]}){
			$nhit{$ne[$n]}++;
			$nhit{$ne[$m]}++;
		    }	
		}
	    }
	}
    }
}

##GF0000007.001	FBgn0264775	chr3R:15509753-15509899:-	Dme
open (OUT, ">$outfile") || die "It cannot open the output file ($outfile)\n";
open (ECL, $infile2) || die "It cannot open the exon file (again)\n";
while (<ECL>){
    chomp($_);
    @l=split(/\t/,$_);
    @t1=split(/\:/,$l[2]);
    $idov=$oid{$l[1]."\t".$t1[1]};
    $maxhit=$exs{$l[0]}-($eclsp{$l[0]}{$l[3]});
    $memb=($nhit{$idov}/$maxhit);
    $memb=sprintf("%.2f",$memb);
    print OUT "$_\t$memb\n";
}
close ECL;
close OUT;








