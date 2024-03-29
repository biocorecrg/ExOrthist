#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# needs to check the counts are right!

my $f_gene_cluster;
my $f_exon_cluster;
my $exonsDB_folder;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions( "gene_clusters=s" => \$f_gene_cluster,
	    "exon_clusters=s" => \$f_exon_cluster,
	    "main_output=s" => \$exonsDB_folder
    );

if (!defined $exonsDB_folder){
    die "
Usage: get_cluster_stats.pl -main_output FOLDER (-gene_clusters FILE -exon_clusters FILE) 

Script to get basic stats from ExOrthist clusters.

OPTIONS
   -main_output FOLDER/    Path to the output folder of main.nf [mandatory]

   -gene_clusters FILE     Alternative gene clusters file [optional]
   -exon_clusters FILE     Alternative exon clusters file [optional]

    
";
}

my @exon_files=glob("$exonsDB_folder/*/*_overlap_CDS_exons.txt"); # contains the ANNOTATED exons only (not lifted)
die "[Aborted] It cannot find any files with exon info\n" if $#exon_files < 0;

my $f_overlapping_EXs_by_species="$exonsDB_folder/overlapping_EXs_by_species.tab";

### defines the files:
if (!defined $f_gene_cluster){
    $f_gene_cluster = "$exonsDB_folder/gene_cluster_file.gz";
} 
if (!defined $f_exon_cluster){
    $f_exon_cluster = "$exonsDB_folder/EX_clusters.tab";
}
die "The gene cluster file ($f_gene_cluster) could not be found\n" unless (-e $f_gene_cluster);
die "The exon cluster file ($f_exon_cluster) could not be found\n" unless (-e $f_exon_cluster);


#GF00001 Bta ENSBTAG00000045550 SPAN6
my %gene_to_cluster=(); 
my %gene_to_species=();

if ($f_gene_cluster=~/\.gz$/){
    open (LIST, "gunzip -c $f_gene_cluster |") || die "It cannot open the gene clusters ($f_gene_cluster file)\n";
} else {
    open (LIST, $f_gene_cluster) || die "It cannot open the gene clusters ($f_gene_cluster file)\n";
}
while (<LIST>){ ##checking list of consistent events
    chomp($_);
    my @l=split(/\t/,$_); 
    $gene_to_cluster{$l[2]} = $l[0];
    $gene_to_species{$l[2]} = $l[1];
}
close LIST;

#OV_EX_Mm2_186054 ENSMUSG00000067377 133892701-133892766
my %exon_catalog=(); # catalog of exons by overlap
my %tally_exons_all=(); # n of all exons per species
my %tally_exons=(); # n of exons per species in genes in clusters
my %exon_conversion=(); # from gene=coordinate to OV_EX_Mm2_186054.
my %exon_conversion_partial=(); # from gene=coordinate1 and gene=coordinate2 to OV_EX_Mm2_186054.

foreach my $exon_file (@exon_files){
    open (EXONS, $exon_file) || die "It cannot open the exon file $exon_file\n";
    while (<EXONS>){ ##checking list of consistent events
	chomp($_);
	my @l=split(/\t/,$_);
	my $gene=$l[1];
	my ($species)=$l[0]=~/OV\_EX\_(.+)\_/;

	if (!defined $exon_catalog{$l[0]}){
	    $tally_exons_all{$species}++;
	    if ($gene_to_cluster{$gene}){
		$tally_exons{$species}++;
	    }
	    $exon_catalog{$l[0]} = 1;
	} 
	my $exon_ID="$gene=$l[2]";
	$exon_conversion{$exon_ID}=$l[0];

	# also keeps the match to the partial exon coordinates
	my ($co_start,$co_end)=$l[2]=~/(.+?)\-(.+)/;
	my $exon_ID_A="$gene=$co_start";
	my $exon_ID_B="$gene=$co_end";
	$exon_conversion_partial{$exon_ID_A}=$l[0];
	$exon_conversion_partial{$exon_ID_B}=$l[0];
    }
    close EXONS;
}
# EX overlap of hits
my %exon_catalog_hits=(); # catalog of exons by overlap
my %tally_exons_all_hits=(); # n of all exons per species
my %tally_exons_hits=(); # n of exons per species in genes in clusters
my %exon_conversion_hits=(); # from gene=coordinate to OV_EX_Mm2_186054.
unless (-e $f_overlapping_EXs_by_species){
    print "It cannot open the exon file $f_overlapping_EXs_by_species. Issues expected if --bonafide is used\n";
} else {
    open (EXON_HITS, $f_overlapping_EXs_by_species);
}
while (<EXON_HITS>){ ##checking list of hits
    chomp($_);
    my @l=split(/\t/,$_);
    my $gene=$l[1];
    next if $gene =~ /GeneID/; # lost header
    my ($species)=$l[0]=~/OV\_EX\_(.+?)\_/;
    
    if (!defined $exon_catalog_hits{$l[0]}){
	$tally_exons_all_hits{$species}++;
	if ($gene_to_cluster{$gene}){
	    $tally_exons_hits{$species}++;
	}
	$exon_catalog_hits{$l[0]} = 1;
    } 
    my ($co)=$l[2]=~/\:(.+)\:/;
    my $exon_ID="$gene=$co";
    $exon_conversion_hits{$exon_ID}=$l[0];
}
close EXON_HITS;

#GF00002.003 ENSMUSG00000031250 chrX:133859723-133859863:+ Mm2
my %tally_exons_in_clusters=(); # n of exons per species in clusters
my %done_exon=(); my %done_cluster=();
my %tally_by_exon_cluster=(); # to get 1:1:1 strings
my $total_exon_clusters=0;
my %tally_bonafide=();

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
	my $exon_name = $exon_conversion{$id2}; # this is the unique ID (OV_EX_dm6_1)
	
	if ($exon_name){ # excludes Lifted exons
	    $tally_exons_in_clusters{$species}++;
	}
	else {
	    my ($co_start,$co_end)=$coord=~/\:(.+?)\-(.+)\:/;
	    my $exon_ID_A="$gene=$co_start";
	    my $exon_ID_B="$gene=$co_end";

	    if (defined $exon_conversion_partial{$exon_ID_A}){
		$exon_name = $exon_conversion_partial{$exon_ID_A};
		$tally_exons_in_clusters{$species}++;
	    }
	    elsif (defined $exon_conversion_partial{$exon_ID_B}){
		$exon_name = $exon_conversion_partial{$exon_ID_B};
		$tally_exons_in_clusters{$species}++;
	    }
	    else {
		$exon_name = "$exon_conversion_hits{$id2}-HITS"; # for liftover hits
		$tally_bonafide{$species}++;
	    }
	}
	if (!defined $done_exon{$exon_name}){ # does not count redundant exons multiple times
	    $done_exon{$exon_name}=1;
	    ### count by cluster
	    $tally_by_exon_cluster{$exon_cluster}{$species}++;
	    $total_exon_clusters++ if !defined $done_cluster{$exon_cluster};
	    $done_cluster{$exon_cluster}=1;
	}
    }
}


### Gets the default cluster string:
my $default_string1;
my $default_string2;
foreach (sort keys %tally_exons){
    $default_string1.="1";
    $default_string2.="2";
}

my %tally_strings=();
my %tally_missing=();
my $total_strings=0;
foreach my $exon_cluster (sort keys %tally_by_exon_cluster){
    my $string;
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

print "Summary statistics of exon orthogroups (OGs)\n"; 
print "Species\tTotal annotated CDS exons\tExons in gene OGs\tExons in OGs\t% recovered\n";
foreach my $species (sort keys %tally_exons){
    my $perc_covered = sprintf ("%.2f",100*($tally_exons_in_clusters{$species}/$tally_exons{$species}));
    print "$species\t$tally_exons_all{$species}\t$tally_exons{$species}\t$tally_exons_in_clusters{$species}\t$perc_covered\%\n";
}


print "\nNon-annotated bonafide exons in OGs:\n";
foreach my $species (sort keys %tally_bonafide){
    print "$species\t$tally_bonafide{$species}\n";
}

$tally_strings{$default_string1}=0 if !defined $tally_strings{$default_string1};
$tally_strings{$default_string2}=0 if !defined $tally_strings{$default_string2};
$tally_strings{MANY_3}=0 if !defined $tally_strings{MANY_3};
$tally_strings{MANY_4}=0 if !defined $tally_strings{MANY_4};
$tally_strings{INCOMPLETE}=0 if !defined $tally_strings{INCOMPLETE};

my $perc_1=sprintf("%.2f",100*$tally_strings{$default_string1}/$total_strings);
my $perc_2=sprintf("%.2f",100*$tally_strings{$default_string2}/$total_strings);
my $perc_3=sprintf("%.2f",100*$tally_strings{MANY_3}/$total_strings);
my $perc_4=sprintf("%.2f",100*$tally_strings{MANY_4}/$total_strings);
my $perc_inc=sprintf("%.2f",100*$tally_strings{INCOMPLETE}/$total_strings);

my $tally_strings12; # for 1:2 or 2:1
my $perc_12;
if (length($default_string1)==2){
    $tally_strings{12}=0 if !defined $tally_strings{12};
    $tally_strings{21}=0 if !defined $tally_strings{21};
    $tally_strings12=$tally_strings{12}+$tally_strings{21};
    $perc_12=sprintf("%.2f",100*$tally_strings12/$total_strings);    
}

### Number of other clusters
my $other_strings = $total_strings-$tally_strings{$default_string1}-$tally_strings{$default_string2}-$tally_strings{MANY_3}-$tally_strings{MANY_4};
$other_strings = $other_strings - $tally_strings12 if length($default_string1)==2;
my $perc_others = sprintf("%.2f",100*$other_strings/$total_strings);

print "\nExon OG type\tNumber\t\% from total OGs\n";
print "Total exon OGs\t$total_exon_clusters\t100\%\n";
print "OGs 1:1\t$tally_strings{$default_string1}\t$perc_1\%\n";
print "OGs 1:2/2:1\t$tally_strings12\t$perc_12\%\n" if length($default_string1)==2;
print "OGs 2:2\t$tally_strings{$default_string2}\t$perc_2\%\n";
print "OGs >=3 exons/species\t$tally_strings{MANY_3}\t$perc_3\%\n";
print "OGs >=4 exons/species\t$tally_strings{MANY_4}\t$perc_4\%\n";
print "Other OGs\t$other_strings\t$perc_others\%\n" if length($default_string1)==2;

if (length($default_string1)>2){
    print "OGs missing species\t$tally_strings{INCOMPLETE}\t$perc_inc\%\n";
    print "\nMissing species\tNumber\t\% from missing\n";
    foreach my $species (sort keys %tally_exons){
	my $perc = "NA";
	$tally_missing{$species}=0 if !defined $tally_missing{$species};
	if ($tally_strings{INCOMPLETE}>0){
	    $perc = sprintf("%.2f",100*$tally_missing{$species}/$tally_strings{INCOMPLETE});
	    print "$species\t$tally_missing{$species}\t$perc\%\n";
	}
	else {
	    print "$species\t$tally_missing{$species}\tNA\n";
	}
    }
}
