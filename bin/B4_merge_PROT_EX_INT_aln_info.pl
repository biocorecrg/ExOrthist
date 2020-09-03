#!/usr/bin/env perl
use warnings;
use strict;

##Program for joining all the score files for further processing##
my $folder=$ARGV[0];

print "Joining all Score files\n";
print "1. Protein scores\n";
my $i1=$folder."/PROT_aln_features_";
my $o1=$folder."/h1.txt";
my $tmp=$folder."/tmp.txt";
my $s1=$folder."/all_PROT_aln_features.txt";
system("cat $i1*_part* | head -1 > $o1"); ##getting header
system("cat $i1*_part* | grep -v 'CID' | sort -k1n | uniq  > $tmp"); ##getting the rest of the file
system("cat $o1 $tmp > $s1"); ##getting file of aligned proteins

print "2. Exon scores\n";
my $i2=$folder."/EX_aln_features_";
my $i3=$folder."/realigned_EXs_";
my $o2=$folder."/h2.txt";
my $tex=$folder."/tex.txt";
my $s2=$folder."/all_EX_aln_features.txt";
`head -1 $i2*part_1.txt > $o2`; ##getting header
`cat $i2*_part_* $i3* | grep -v 'CID' | sort -k1 | uniq  > $tex`;

open (FILE, "$tex") || die "It cannot open the file with aligned exons ($tex)\n";
open (TMP, ">$tmp") || die "It cannot open the temporary file ($tmp)\n"; # tmp.txt
while (<FILE>){
    chomp($_);
    my @l=split(/\t/,$_);
    if ($l[7]==1){ # it selects only the entries with a 1 in column 8...
	print TMP "$_\n"; 
    }
}
close (FILE);
close (TMP);
`cat $o2 $tmp > $s2`; ##getting file of exon scores

print "3. Intron scores\n";
my $i4=$folder."/INT_aln_features_";
my $s3=$folder."/all_INT_aln_features.txt";
my $o3=$folder."/h3.txt";
`cat $i4*_part_1.txt | head -1 > $o3`; ##getting header (new tmp file)
`cat $i4*_part_* | grep -v 'CID' | sort -k1 | uniq  > $tmp`; # re-creating tmp.txt
`cat $o3 $tmp > $s3`;  ## getting file of intron scores 

print "4. Merging EXINT ALNs\n";
system "cat $folder/*part_*.ALL.aln | gzip > $folder/EXINT_aln.gz";

print "5. Cleaning up intermediate files\n";
`rm $o1 $o2 $o3 $tmp $tex`; ## removing temporal files
`rm *part_* realigned_EXs_*`; # more temp files
