#!/usr/bin/perl
#use warnings;
use strict;

my $ext1=$ARGV[0]; ##Folder vastdb info
my $in1=$ARGV[1]; ##exon clusters
my $out=$ARGV[2]; ##exon clusters with vastids
my ($id1, $sp, $sp1, $l, $id2, $n);
my $n1=$ext1."REFERENCE*ALL*ANNOT*.tab";;
my $n2=$ext1."*ID.names.txt";
my $list1=`ls $n1`;
my $list2=`ls $n2`;
my @files1=split(/\n/,$list1);
my @files2=split(/\n/,$list2);
my (%vid, %gid,%name, %sps);
my (@l, @l1, @line, @sp, @k, @t);
my ($posmic,$n, $sp);
my $tmpcid;
my %species;

#REFERENCE EXON TABLE
#Cdc45	MmuEX0010430	chr16:18811398-1881145760
my $cid=0;
for ($l=0; $l<scalar(@files1); $l++){
    open (INFILE,"$files1[$l]");
    while (<INFILE>){
	chomp($_); 
	@line=split(/\t/,$_);
	if ($line[1]=~/EX/){
	    $vid{$line[2]}=$line[1];
	    $vid{$line[9]}=$line[1];	
	}
    }
}
close (INFILE);

##GENE_IDs GENE_NAMES table
for ($l=0; $l<scalar(@files2); $l++){
    open (INFILE,"$files2[$l]");
    while (<INFILE>){
	chomp($_); 
	@line=split(/\t/,$_);
	$gid{$line[0]}=$line[1];
    }
}
close (INFILE);

my ($eid,$gn,$s, $tmp, $size, $cr, $t1);
my (@ts,@c);
open (OUT, ">$out");
print OUT "ExCID\tSps\tGeneID\tGene_name\tExon_coords\tExon_length\tVastID\n";x
open (INONE,"$in1");
# Format: 7	ENSDARG00000057688	chr6:12656155-12656227:-	Dre ???
$t1="";
while (<INONE>){
    chomp($_); 
    @line=split(/\t/,$_);
    @c=split(/\:/,$line[2]);
    $cr=$c[0].":".$c[1];
    if ($gid{$line[1]}){ $gn=$gid{$line[1]}; } else { $gn="no_gene_name"; }
    if ($vid{$cr}){ $eid=$vid{$cr}; } else { $eid="no_vast_id";  }
    $sp1=$line[3];
    #print "##$_\n";
    @ts=split(/\-/,$c[1]); $size=$ts[1]-$ts[0]+1;	
    if ($t1 && ($t1 ne $line[0])){ print OUT "\n"; }
    print OUT "$line[0]\t$sp1\t$line[1]\t$gn\t$cr\t$size\t$eid\n";
    $t1=$line[0];
}
close OUT;
close INONE;
