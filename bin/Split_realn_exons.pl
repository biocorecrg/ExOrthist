#!/usr/bin/env perl
use warnings;
#use strict;

my$p=1000;
my $dir=$ARGV[0];
$p=$ARGV[1];

$file=$dir."/tmp.txt";
`cat $dir/exons_to_split_part_* > $file`;
#`rm $dir/exons_to_realign_part_*`;

$h=0;
$c=0;
$t=1;
$nrow = 0;
open (FILE, "$file");
while (<FILE>){
    $nrow++;
	chomp($_);
	if ($_=~/CID/) { $h=$_; }
	else { $c++; 
	if ($c==1){
		$of=$dir."/exons_to_realign_part_".$t.".txt";
		open (OUT, ">>$of");
		if (!$hd{$t}){
			print OUT "$h\n";
			$hd{$t}=1;
		}
		print OUT "$_\n";
	}
	elsif ($c>$p){
		$t++; 
		$of=$dir."/exons_to_realign_part_".$t.".txt";
		open (OUT, ">>$of");
		if (!$hd{$t}){
			print OUT "$h\n";
			$hd{$t}=1;
		}
		print OUT "$_\n";
		$c=0;
	}
	else {
		$of=$dir."/exons_to_realign_part_".$t.".txt";
		open (OUT, ">>$of");
		print OUT "$_\n";	
	}	
	}

}
if ($nrow == 1) {
	$of=$dir."/exons_to_realign_part_".$t.".txt";
	open (OUT, ">>$of");
	print OUT "$h\n";
}


