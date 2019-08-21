#!/usr/bin/perl
#use warnings;
#use strict;
use Cwd qw(chdir);
##
my $f=$ARGV[0]; ##Directory with all the parts of exon scores
my $np=$ARGV[1]; ##Number of parts
my $bin=$ARGV[2]; ##bin directory

if ($np<=9){
	print "1. Preparing files\n\n";
	for ($l=1; $l<=$np; $l++){
		$dir=$f."/PART_0$l/";
		#print "$dir\n";
		chdir ("$dir");
		#print  "\nCurrent Directory is $ENV{PWD} \n";
		`perl $bin/Get_Pre_cluster_exons_V2.pl`;
	}
	print "2. Submiting jobs for clustering\n\n";
	for ($l=1; $l<=$np; $l++){
		$dir=$f."/PART_0$l/";
		chdir ("$dir");
		`qsub array_job.sh`;
	}
}
else {	
	print "1. Preparing files\n\n";
	for ($l=1; $l=9; $l++){
		$dir=$f."/PART_0$l/";
		chdir ("$dir");
		`perl $bin/Get_Pre_cluster_exons_V2.pl`;
	}
	for ($l=10; $l<=$np; $l++){
		$dir=$f."/PART_$l/";
		chdir ("$dir");
		`perl $bin/Get_Pre_cluster_exons_V2.pl`;
	}
	print "2. Submiting jobs for clustering\n\n";		
	for ($l=1; $l=9; $l++){
		$dir=$f."/PART_0$l/";
		chdir ("$dir");
		`qsub array_job.sh`;
	}
	for ($l=0; $l<=$np; $l++){
		$dir=$f."/PART_0$l/";
		chdir ("$dir");
		`qsub array_job.sh`;
	}
}

