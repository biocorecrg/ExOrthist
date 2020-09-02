#!/usr/bin/env perl
use warnings;
use strict;

#my $p=1000;
my $dir=$ARGV[0];
my $p = $ARGV[1];

my $file=$dir."/tmp.txt";
`cat $dir/EXs_to_split_part_* > $file`;
#`rm $dir/exons_to_realign_part_*`; # temporary files? Could be removed?

my $h=0;
my $c=0;
my $t=1;
my $nrow = 0;
my $of;
my %hd;
open (FILE, $file) || die "It cannot open the merged file of exons to re-align\n";
while (<FILE>){
    $nrow++;
    chomp($_);
    if ($_=~/CID/) { $h=$_; }
    else { 
	$c++; 
	if ($c==1){
	    $of=$dir."/EXs_to_realign_part_".$t.".txt";
	    open (OUT, ">>$of");
	    if (!$hd{$t}){
		print OUT "$h\n";
		$hd{$t}=1;
	    }
	    print OUT "$_\n";
	}
	elsif ($c>$p){
	    $t++; 
	    $of=$dir."/EXs_to_realign_part_".$t.".txt";
	    open (OUT, ">>$of");
	    if (!$hd{$t}){
		print OUT "$h\n";
		$hd{$t}=1;
	    }
	    print OUT "$_\n";
	    $c=0;
	}
	else {
	    $of=$dir."/EXs_to_realign_part_".$t.".txt";
	    open (OUT, ">>$of");
	    print OUT "$_\n";	
	}	
    }
}
if ($nrow == 1) {
    $of=$dir."/EXs_to_realign_part_".$t.".txt";
    open (OUT, ">>$of");
    print OUT "$h\n";
}
