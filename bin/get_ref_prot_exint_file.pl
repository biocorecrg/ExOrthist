#!/usr/bin/perl
use warnings;
use strict;
my $i1=$ARGV[0]; 
my $i2=$ARGV[1]; 
my $out=$ARGV[2]; 
my ($i,$j,$key,$sp1, $sp2, $id, $c);
my %ref;
my %ids;
my $b=0;
my (@l, @l1, @line, @t1, @t2);
#BL61971	BL61971_cuf3	BL61971_cuf3
#ENSG00000000003	ENST00000373020	ENSP00000362111
open (INFILE,"$i1");
while (<INFILE>){
	chomp($_); 
	@line=split(/\s+/,$_);
	$ref{$line[1]}=1;
}
close (INFILE);
open (OUT, ">$out");
open (INTWO,"$i2");
#>ENSP00000379387|ENSG00000185219  144.0 200.2 232.0 274.2 311.2
while (<INTWO>){
	chomp($_);
	if ($_=~/\>/){ 
		$b=0;
		@line=split(/\|/,$_);
		$line[0]=~s/\>//;
		if ($ref{$line[0]}) { $b=1; $id=$_; }
	}
	elsif ($b==1) {
		print OUT "$id\n$_\n";
	}	
}
close (INTWO);
close (OUT);









