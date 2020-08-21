#!/usr/bin/env perl
use warnings;
use strict;

##Program for joining all the score files for further processing##
my $folder=$ARGV[0];
my $log=$folder."/log_join_all_scores.txt";
open (OUT, ">$log");
print OUT "Joining all Score files\n\n";
print OUT "1. Protein scores\n\n";
my $i1=$folder."/score_proteins_";
my $o1=$folder."/h.txt";
my $tmp=$folder."/tmp.txt";
my $s1=$folder."/Aligned_proteins.txt";
system("cat $i1*_part* | head -1 > $o1"); ##getting header
system("cat $i1*_part* | grep -v 'CID' | sort -k1n | uniq  > $tmp"); ##getting the rest of the file
system("cat $o1 $tmp > $s1"); ##getting file of aligned proteins

print OUT "2. Exon scores\n\n";
my $i2=$folder."/score_exons_";
my $i3=$folder."/realigned_exons_";
my $o2=$folder."/h.txt";
my $tex=$folder."/tex.txt";
my $s2=$folder."/Score_all_exons.txt";
`head -1 $i2*part_1.txt > $o2.txt`; ##getting header
`cat $i2*_part_* $i3* | grep -v 'CID' | sort -k1 | uniq  > $tex`;
open (FILE, "$tex");
open (TMP, ">$tmp");
while (<FILE>){
    chomp($_);
    my @l=split(/\t/,$_);
    if ($l[7]==1){
	print TMP "$_\n";
    }
}
close (FILE);
close (TMP);

`cat $o2 $tmp > $s2`; ##getting file of exon scores
my $i4=$folder."/score_introns_";
my $s3=$folder."/Score_all_introns.txt";
print OUT "3. Intron scores\n\n";
`cat $i4*_part_1.txt | head -1 > $o2`; ##getting header
`cat $i4*_part_* | grep -v 'CID' | sort -k1 | uniq  > $tmp`;
`cat $o2  $tmp > $s3`;  ##getting file of intron scores 
`rm $o1 $o2 $tmp $tex`; ##removing temporal files
print OUT  "Done!!!\n\n";
