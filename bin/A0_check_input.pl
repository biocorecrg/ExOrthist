#! /usr/bin/env perl

=head1 NAME

=head1 SYNOPSIS

  perl check_input.pl [-e evodists.txt] [-g gtf files] [-f fasta files] [-h help]

=head1 DESCRIPTION

  This script check the correctness of the input for the ExOrthist pipeline
  It returns ok or exit

  
Typical usage is as follows:

  % perl check_input.pl -e evodists.txt -g test/GTF/*.gtf.gz -f test/GENOMES/*.fasta.gz
  
=head2 Options

The following options are accepted:

 --e=<evodists file>   	File path.

 --g=<GTF files>	Path to GTF files (glob)

 --f=<Fasta files>	Path to Fasta files (glob)

 --help                 This documentation.


=head1 AUTHOR

Luca Cozzuto <luca.cozzuto@crg.es> 

=cut
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Pod::Usage;
use Getopt::Long;

my $USAGE = "perl check_input.pl [-e evodists file] [-g gtf files] [-f fasta files] [-h help]";

my ($fastas,$gtf_files,$evodists,$show_help);

&GetOptions(    	
			'fasta|f=s'	=> \$fastas,
			'gtf|g=s'   	=> \$gtf_files,
			'evodists|e=s'	=> \$evodists,
			'help|h'        => \$show_help
			)
  or pod2usage(-verbose=>2);
pod2usage(-verbose=>2) if $show_help;

if (!$fastas) { die ("Please specify fasta files")}
if (!$gtf_files) { die ("Please specify gtf files")}
if (!$evodists) { die ("Please specify evodists file")}

my $fastafiles = getIds($fastas);
my $gtffiles = getIds($gtf_files);

open (my $handle, $evodists) or die ("Can't read $evodists");
my %evo_ids;
while(<$handle>) {
 	my $line=$_;
 	chomp($line);
 	my @fields = split("\t", $line);
	$evo_ids{$fields[0]} = $fields[0];
	$evo_ids{$fields[1]} = $fields[1];
}

my %all_ids;

foreach my $fid (keys(%$fastafiles)) {
	$all_ids{$fid} = $fid;
}
foreach my $gid (keys(%$gtffiles)) {
	$all_ids{$gid} = $gid;
}
foreach my $eid (keys(%evo_ids)) {
	$all_ids{$eid} = $eid;
}
my $status = 0;

foreach my $aid (keys(%all_ids)) {
	if (!$fastafiles->{$aid}) {
		print "*** WARNING fasta file corresponding to $aid is missing!!! *** \n";
		$status = 1;
	}

	if (!$gtffiles->{$aid}) {
		print "*** WARNING gtf file corresponding to $aid is missing!!! *** \n";
		$status = 1;
	}
	
	if (!$evo_ids{$aid}) {
		print "*** WARNING $aid is missing within $evodists file! *** \n";
		$status = 1;
	}
}

exit ($status);

print Dumper \%evo_ids;
print Dumper $fastafiles;
print Dumper $gtffiles;

sub getIds {
	my ($glob_pattern) = @_;
	my @files = glob($glob_pattern);
	my %ids;
	my @extension = split(/\*/, $glob_pattern); 
	foreach (@files) {
        my $id = basename($_,$extension[1]);
		$ids{$id} = $_;
	}
	return(\%ids);
}



