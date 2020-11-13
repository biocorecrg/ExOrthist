#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $chain_file; # to do liftover
my $gene_clusters; 
my $annot_sp1; # GTF named Sp1.gtf
my $annot_sp2; # GTF named Sp2.gtf
my $canonical_ss; # requires at last 1 canonical splice site
my $line_to_parse="CDS"; # or exon
my $help;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(  "annot_sp1=s" => \$annot_sp1,
	     "annot_sp2=s" => \$annot_sp2,
	     "gene_clusters=s" => \$gene_clusters,
	     "chain_file=s" => \$chain_file,
	     "type=s" => \$line_to_parse,
	     "canonical_ss" => \$canonical_ss,
	     "help" => \$help
    );

if (!defined $annot_sp1 || defined $help){
    die "
Usage: GetLiftOverFile.pl OPTIONS

Options:
   -annot_sp1 GTF/LIST     GTF (Sp1.gtf) or list of exons (Sp1.exons) for species 1 (query).
                              If a list is provided (tsv, Sp1.exons): GeneID  chr:start-end:strand
   -annot_sp2 GTF          GTF for species 2 (target).
   -chain_file FILE        LiftOver chain alignment from Sp1ToSp2.
   -gene_clusters FILE     Gene orthology clusters including Sp1 and Sp2. 
#                              Format (tsv): ClusterID   Species    GeneID
   -type CDS/exon          If a GTF is provided for annot_sp1, parses either CDS or exon lines (def = CDS).
   -canonical_ss           If flag is active, it will require at last one cannonical splice site (def = off).
                              If active, it will need the filter Sp2.fasta in the same directory.


*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

# Checks all files exist
die "It cannot find $annot_sp1\n" unless (-e $annot_sp1);
die "It cannot find $annot_sp2\n" unless (-e $annot_sp2);
die "It cannot find $chain_file\n" unless (-e $chain_file);
die "It cannot find $gene_clusters\n" unless (-e $gene_clusters);

# Extracts roots
my $sp1 = $annot_sp1;
my $sp2 = $annot_sp2;
$sp1=~s/\..+//;
$sp1=~s/.+\///;
$sp2=~s/\.gtf//;
$sp2=~s/.+\///;

my $genome_sp2 = "$sp2.fasta";
if (defined $canonical_ss){
    die "It cannot find genome for species 2, and canonical_ss is active\n" unless (-e $genome_sp2);
}

# Processes cluster file
my %gene_to_cluster=();
my %sp_check=();
open (CLUSTERS, $gene_clusters) || die "It cannot open gene clusters\n";
while (<CLUSTERS>){
    chomp($_);
    my @t=split(/\t/,$_);
    my $species=$t[1];
    if ($species eq $sp1 || $species eq $sp2){
	$gene_to_cluster{$t[2]}=$t[0];
        $sp_check{$species}++;
    }
}
close CLUSTERS;

# Checks species
die "It did not identify any gene ortolog for $sp1 (correct file names? $annot_sp1)\n" if (!defined $sp_check{$sp1});
die "It did not identify any gene ortolog for $sp2 (correct file names? $annot_sp1)\n" if (!defined $sp_check{$sp2});

# Extracts exons from Sp1
print "Extracting exons from $annot_sp1\n";
my %done_exons=();
my $temp_exons_sp1 = "temp_exons1_$sp1-$sp2.bed";
open (TEMP1, ">$temp_exons_sp1") || die "It cannot open $temp_exons_sp1\n";
open (GTF_SP1, $annot_sp1) || die "It cannot open $annot_sp1\n";
while (<GTF_SP1>){
    chomp($_);
    my @t = split(/\t/,$_);
    
    if ((defined $t[3] && $t[3]=~/\d/) || $annot_sp1=~/\.gtf/){
	if (($t[2] eq "exon" && $line_to_parse eq "exon") || ($t[2] eq "CDS" && $line_to_parse eq "CDS")){
	    my ($gene)=$_=~/gene_id \"(.+?)\"/;
	    my $exon="$gene=$t[0]:$t[3]-$t[4]:$t[6]"; # chr:start-end:strand
	    
	    unless (defined $done_exons{$exon}){
		$t[3]--;
		print TEMP1 "$t[0]\t$t[3]\t$t[4]\t$exon\t1\t$t[6]\n";
		$done_exons{$exon}=1;
	    }	
	}
    }
    else {
        my $gene = $t[0];
        my $exon = "$gene=$t[1]";
	
        unless (defined $done_exons{$exon}){
            my ($chr,$i,$f,$strand) = $exon =~ /\=(.+?)\:(\d+)\-(\d+)\:([\+\-])/;
            $i--;
            print TEMP1 "$chr\t$i\t$f\t$exon\t1\t$strand\n";
            $done_exons{$exon}=1;
        }
    }
    
}
close GTF_SP1;

# Does the liftover to Sp2
print "Doing liftOver from $sp1 to $sp2\n";
my $temp_lifted_exons = "temp_exons1_$sp1-$sp2-in_$sp2.bed";
system "liftOver -minMatch=0.10 -multiple -minChainT=200 -minChainQ=200 $temp_exons_sp1 $chain_file $temp_lifted_exons /dev/null";

# If cannonical_ss, filters here
my $temp_lifted_exons_dinuc = "temp_exons1_$sp1-$sp2-in_$sp2-F_dinuc.bed"; # remaining lifted exons AFTER filter
if (defined $canonical_ss){
    # Loads genome sequence
    my %SEQ; my $current_chr;
    open (SEQ, $genome_sp2) || die "It cannot open the genome file ($genome_sp2)\n";
    while (<SEQ>){
	chomp;
	if (/\>/){
	    ($current_chr)=/\>(.+)/;
	} else {
	    $SEQ{$current_chr}.=$_ if $current_chr;
	}
    }
    close SEQ;
    
    # Gets the dinucleotide sequences and does the filtering
    my %done_chr;
    open (TEMP_LIFT, $temp_lifted_exons) || die "It cannot open the lifted exons ($temp_lifted_exons)\n";
    open (O_FILTERED, ">$temp_lifted_exons_dinuc") || die "It cannot open the output filtered file ($temp_lifted_exons_dinuc)\n";
    while (<TEMP_LIFT>){
	chomp;
	my @t = split(/\t/);
	my $chr = $t[0];
	my $start = $t[1];
	my $end = $t[2];
	my $strand = $t[5];
	my $ID = $t[3];
	
	my ($left_seq,$right_seq);
	if ($SEQ{$chr}){ # weird chromosomes
	    $left_seq = substr($SEQ{$chr},$start-2,2);
	    $right_seq = substr($SEQ{$chr},$end,2);
	    
	    if ($strand eq "-"){
		$left_seq = join ("", reverse split (//, $left_seq));
		$left_seq =~ tr/ACGTacgt/TGCATGCA/;
		$right_seq = join ("", reverse split (//, $right_seq));
		$right_seq =~ tr/ACGTacgt/TGCATGCA/;
	    }
	    
	    ### Fiters the file (liftover inverts strand sometimes, so this cases will be discarded
	    if (($strand eq "+" && ($left_seq eq "AG" || $right_seq=~/G[TC]/)) ||
		($strand eq "-" && ($left_seq=~/G[TC]/ || $right_seq eq "AG"))){
		print O_FILTERED "$_\n";
	    }
	}
	else {
	    print "No sequence data for $chr\n" unless $done_chr{$chr};
	    $done_chr{$chr}=1;
	}
    }
    close TEMP_LIFT;
    close O_FILTERED;
}
else { # no filter
    system "cp $temp_lifted_exons $temp_lifted_exons_dinuc"; # simply makes a copy
}

# Intersects with coordinates from Sp2 orthologs
print "Extracting exons from $annot_sp2\n";
my %gene2_start=();
my %gene2_end=();
my %gene2_chr=();
my %gene2_strand=();

my $temp_genes_sp2 = "temp_genes2_$sp1-$sp2.bed";
open (TEMP2, ">$temp_genes_sp2") || die "It cannot open $temp_genes_sp2\n";
open (GTF_SP2, $annot_sp2) || die "It cannot open $annot_sp2\n";
while (<GTF_SP2>){
    chomp($_);
    my @t = split(/\t/,$_);
    
    if ($t[2] eq "exon"){
	my ($gene)=$_=~/gene_id \"(.+?)\"/;
	if (defined $gene_to_cluster{$gene}){ # the gene is in a cluster
	    $gene2_start{$gene}=$t[3] if (!defined $gene2_start{$gene} || $t[3] <= $gene2_start{$gene});
	    $gene2_end{$gene}=$t[4] if (!defined $gene2_end{$gene} || $t[4] >= $gene2_end{$gene});
	    $gene2_chr{$gene}=$t[0];
	    $gene2_strand{$gene}=$t[6];
	}
    }
}
close GTF_SP2;

foreach my $gene (sort keys %gene2_start){
    print TEMP2 "$gene2_chr{$gene}\t$gene2_start{$gene}\t$gene2_end{$gene}\t$gene\t1\t$gene2_strand{$gene}\n";
}

print "Doing coordinate intersect for $sp2\n";
my $temp_intersect  = "temp_intersect_$sp1-$sp2.tab";
system "bedtools intersect -a $temp_lifted_exons_dinuc -b $temp_genes_sp2 -wao > $temp_intersect";

# Evaluates intersects
open (O, ">Liftover_exons-$sp1-$sp2.tab");
open (INTERSECTS, $temp_intersect) || die "It cannot open $temp_intersect\n";
while (<INTERSECTS>){
    chomp($_);
    my @t=split(/\t/,$_);
    
    #checks both gene IDs for same cluster
    if ($t[7] ne "-1"){
	my ($gene1,$exon1) = $t[3] =~ /(.+?)\=(.+)/;
	my $gene2 = $t[9];
        if (defined ($gene_to_cluster{$gene1}) && defined ($gene_to_cluster{$gene2})){
     	    if ($gene_to_cluster{$gene1} eq $gene_to_cluster{$gene2}){
	        $t[1]++;
	        my $exon2 = "$t[0]:$t[1]-$t[2]:$t[5]";
 	        print O "$gene1\t$exon1\t$gene2\t$exon2\t$sp1\t$sp2\n";
            }
	}
    }
}
close O;
close INTERSECTS;

system "rm temp*$sp1-$sp2*";
