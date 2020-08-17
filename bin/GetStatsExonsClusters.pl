#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $f_gene_cluster;
my $f_exon_cluster;
my $exonsDB_folder;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions( "gene_clusters=s" => \$f_gene_cluster,
	    "exon_clusters=s" => \$f_exon_cluster,
	    "exons_db=s" => \$exonsDB_folder
    );

if (!defined $f_gene_cluster || !defined $f_exon_cluster || !defined $exonsDB_folder){
    die "
Usage: Get_stats_exon_cls.pl -gene_clusters FILE -exon_clusters FILE -exons_db FOLDER

Script to get basic stats from ExOrthist clusters.
    
";
}

my @exon_files=glob("$exonsDB_folder/*/*_prot_exons_overlap.txt");
die "[Aborted] It cannot find any files with exon info\n" if $#exon_files < 0;

#GF00001 Bta ENSBTAG00000045550 SPAN6
my %gene_to_cluster=(); 
my %gene_to_species=();

open (LIST, $f_gene_cluster) || die "It cannot open the gene clusters ($f_gene_cluster file)\n";
while (<LIST>){ ##checking list of consistent events
    chomp($_);
    my @l=split(/\t/,$_); 
#    if ($l[4] ne "Intronless"){ # what file is this??
	$gene_to_cluster{$l[2]} = $l[0];
	$gene_to_species{$l[2]} = $l[1];
#    }
}
close LIST;

#OV_EX_Mm2_186054 ENSMUSG00000067377 133892701-133892766
my %exon_catalog=();
my %tally_exons=(); # n of exons per species in genes in clusters
my %exon_conversion=(); # from gene=coordinate to OV_EX_Mm2_186054.

foreach my $exon_file (@exon_files){
    open (EXONS, $exon_file) || die "It cannot open the exon file $exon_file\n";
    while (<EXONS>){ ##checking list of consistent events
	chomp($_);
	my @l=split(/\t/,$_);
	my $gene=$l[1];

	if ($gene_to_cluster{$gene}){
	    my $species = $gene_to_species{$gene};
	    if (!defined $exon_catalog{$gene}{$l[0]}){
		$exon_catalog{$gene}{$l[0]} = 1;
		$tally_exons{$species}++;
	    }
	    my $exon_ID="$gene=$l[2]";
	    $exon_conversion{$exon_ID}=$l[0];
	} 
    }
    close EXONS;
}

#GF00002.003 ENSMUSG00000031250 chrX:133859723-133859863:+ Mm2
my %tally_exons_in_clusters=(); # n of exons per species in clusters
my %done_exon=(); my %done_cluster=();
my %tally_by_exon_cluster=(); # to get 1:1:1 strings
my $total_exon_clusters=0;

open (EX_CLUSTERS, $f_exon_cluster) || die "It cannot open the file with exon clusters ($f_exon_cluster)\n";
while (<EX_CLUSTERS>){ 
    chomp($_);
    my @l=split(/\t/,$_);
    my $exon_cluster = $l[0];
    my $gene = $l[1];
    my $coord = $l[2];
    my $species = $l[3];

    if ($gene_to_cluster{$gene}){ # if the gene has orthologs
	my @temp = split(/\:/,$coord);
	my $id2 = "$gene=$temp[1]";
	my $exon_name = $exon_conversion{$id2};
	
	if (!defined $done_exon{$exon_name}){
	    $tally_exons_in_clusters{$species}++;
	    $done_exon{$exon_name}=1;
	    # print OUT "$_\n";

	    ### count by cluster
	    $tally_by_exon_cluster{$exon_cluster}{$species}++;
	    $total_exon_clusters++ if !defined $done_cluster{$exon_cluster};
	    $done_cluster{$exon_cluster}=1;
	}
    }
}

### Gets the default cluster string:
my $default_string1="";
my $default_string2="";
foreach (sort keys %tally_exons){
    $default_string1.="1";
    $default_string2.="2";
}

my %tally_strings=();
my %tally_missing=();
my $total_strings=0;
foreach my $exon_cluster (sort keys %tally_by_exon_cluster){
    my $string="";
    foreach my $species (sort keys %tally_exons){
	$tally_by_exon_cluster{$exon_cluster}{$species}=0 if !defined $tally_by_exon_cluster{$exon_cluster}{$species};
	$tally_by_exon_cluster{$exon_cluster}{$species}=9 if $tally_by_exon_cluster{$exon_cluster}{$species}>9;
	$string.="$tally_by_exon_cluster{$exon_cluster}{$species}";
	$tally_missing{$species}++ if $tally_by_exon_cluster{$exon_cluster}{$species} == 0;
    }
    $tally_strings{$string}++;
    $tally_strings{MANY_3}++ if $string=~/3/;
    $tally_strings{MANY_4}++ if $string=~/4/;
    $tally_strings{INCOMPLETE}++ if $string=~/0/;
    $total_strings++;
}

print "Species\tAnnotated CDS exons\tExons in clusters\t%covered\n";
foreach my $species (sort keys %tally_exons){
    my $perc_covered = sprintf ("%.2f",100*($tally_exons_in_clusters{$species}/$tally_exons{$species}));
    print "$species\t$tally_exons{$species}\t$tally_exons_in_clusters{$species}\t$perc_covered\%\n";
}

my $perc_1=sprintf("%.2f",100*$tally_strings{$default_string1}/$total_strings);
my $perc_2=sprintf("%.2f",100*$tally_strings{$default_string2}/$total_strings);
my $perc_3=sprintf("%.2f",100*$tally_strings{MANY_3}/$total_strings);
my $perc_4=sprintf("%.2f",100*$tally_strings{MANY_4}/$total_strings);
my $perc_inc=sprintf("%.2f",100*$tally_strings{INCOMPLETE}/$total_strings);

print "\nTotal exon clusters\t$total_exon_clusters\n";
print "Clusters $default_string1\t$tally_strings{$default_string1}\t$perc_1\%\n";
print "Clusters $default_string2\t$tally_strings{$default_string2}\t$perc_2\%\n";
print "Clusters >=3 exons/species\t$tally_strings{MANY_3}\t$perc_3\%\n";
print "Clusters >=4 exons/species\t$tally_strings{MANY_4}\t$perc_4\%\n";
print "Incomplete clusters\t$tally_strings{INCOMPLETE}\t$perc_inc\%\n";
foreach my $species (sort keys %tally_exons){
    my $perc = sprintf("%.2f",100*$tally_missing{$species}/$tally_strings{INCOMPLETE});
    print "  - $species\t$tally_missing{$species}\t$perc\%\n";
}
