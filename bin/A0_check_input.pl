#!/usr/bin/perl
use Getopt::Long;
use warnings;
use strict;

my $evo_file;
my $gtfs;
my $gtfs_suf;
my $fastas;
my $fastas_suf;
my $gene_clusters;
my $helpFlag;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(               "e=s" => \$evo_file,
                          "g=s" => \$gtfs,
			  "gs=s" => \$gtfs_suf,
			  "f=s" => \$fastas,
			  "fs=s" => \$fastas_suf,
			  "c=s" => \$gene_clusters,
			  "h" => \$helpFlag
    );


if (!defined $evo_file || !defined $gtfs || !defined $fastas || !defined $gene_clusters || $helpFlag){
    die "
Usage: A0_check_input.pl -e evo_dist_file -g gtf_folder -gs gtf_suffix -f fasta_folder -fs fasta_suffix -c gene_clusters_file [-h]

OPTIONS
    -e evo_dist_file         File containing all pairwise exon distances
    -g gtf_folder/           Path to folder containing all gtf files
    -gs gtf_suffix	     Suffix at the end of the GTF filename (following the speciesID)
    -f fasta_folder/         Path to folder containing all fasta files
    -fs fasta_suffix	     Suffix at the end of the fasta files (following the speciesID)
    -c gene_cluster_file     File with gene clusters
    -h                       Prints this help message


Bugs: Manuel Irimia

";
}

#Generate whole path for gtf and fasta files
$gtfs = $gtfs."/*".$gtfs_suf;
$fastas = $fastas."/*".$fastas_suf;

my %species;
open (EVO, $evo_file) || die "It cannot open the file with pairwise evo distances\n";
while (<EVO>){
    chomp;
    my @t=split(/\t/);
    # stores species non-redundantly
    $species{$t[0]}++;
    $species{$t[1]}++;
}
close EVO;

### Checks if all species appear in the same number of pairwise comparisons
my $total_file_warnings=0;
my $first_value;
foreach my $sp (sort keys %species){
    $first_value = $species{$sp} if !defined $first_value;
    if ($first_value ne $species{$sp}){
	print "*** WARNING: the number of pairwise comparisons in the evodist file is not the same for all species\n";
	$total_file_warnings=1; # only counted once
    }
}

my @gtf_files=glob($gtfs);
my @fasta_files=glob($fastas);

my ($gtf_pre,$gtf_post)=split(/\*/,$gtfs);
my ($fasta_pre,$fasta_post)=split(/\*/,$fastas);
### Does the basic GTF/FASTA file presence test
foreach my $sp (sort keys %species){
    my ($test_gtf,$test_fasta);
    foreach my $temp_gtf (@gtf_files){
	$test_gtf=1 if $temp_gtf eq $gtf_pre.$sp.$gtf_post;
	$test_gtf=0 unless (-e $temp_gtf);
    }
    foreach my $temp_fasta (@fasta_files){
	$test_fasta=1 if $temp_fasta eq $fasta_pre.$sp.$fasta_post;
	$test_fasta=0 unless (-e $temp_fasta);
    }
    print "*** WARNING: GTF file is missing for $sp ($sp$gtf_post)\n" if !$test_gtf;
    print "*** WARNING: FASTA file is missing for $sp ($sp$fasta_post)\n" if !$test_fasta;
    $total_file_warnings++ if !$test_gtf;
    $total_file_warnings++ if !$test_fasta;
}

### Loads all the geneIDs with CDS for each species
my %valid_gene; 
foreach my $sp (sort keys %species){
    my $temp_gtf = $gtf_pre.$sp.$gtf_post;
    if ($temp_gtf=~/\.gz$/){
    	open (GTF, "gunzip -c $temp_gtf |");
    } else {
    	open (GTF, $temp_gtf);
    }
    while (<GTF>){
	chomp;
	my @t=split(/\t/);
	if ($t[2] eq "CDS"){
	    my ($gene_id) = /gene_id \"(.+?)\"/;
	    $valid_gene{$sp}{$gene_id}=1;
	}
    }
    close GTF;
}
#print "\n";

### Parses the cluster file
my $total_id_warnings=0;
if ($gene_clusters=~/\.gz$/){
    open (CLUSTER, "gunzip -c $gene_clusters |") || die "It cannot open file with gene clusters ($gene_clusters)\n";
} else {
    open (CLUSTER, $gene_clusters) || die "It cannot open file with gene clusters ($gene_clusters)\n";
}
while (<CLUSTER>){
    chomp;
    my @t=split(/\t/);
    my $sp=$t[1];
    my $gene_id=$t[2];

    if (!$valid_gene{$sp}{$gene_id}){
	print "*** WARNING: gene_id ($gene_id) is not a protein coding gene and/or is not present in GTF ($sp$gtf_post)\n";
	$total_id_warnings++;
    }
}
close CLUSTER;

print "Total file warnings: $total_file_warnings\n";
print "Total gene_id warnings: $total_id_warnings\n";

open (VERSION, "VERSION");
my $version = <VERSION>;
chomp($version);
my $current_time = &time;
print "\nVersion: $version $current_time\n";

die "\nExOrthist is aborting: file warnings are not allowed\n" if $total_file_warnings>0;



sub time {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
    $year += 1900;
    $mon += 1;
    my $datetime = sprintf "%04d-%02d-%02d (%02d:%02d)", $year, $mday, $mon, $hour, $min;
    return "[Date: ".$datetime."]";
}
