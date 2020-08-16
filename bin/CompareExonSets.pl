#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;

my $f_gene_cluster;
my $f_exon_cluster;
my $f_exon_list_sp1;
my $f_exon_list_sp2;
my $sp1;
my $sp2;
my $dPSI_info = "auto";
my $min_dPSI = 15;
my $outFile;
my $helpFlag;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions( "gene_clusters=s" => \$f_gene_cluster,
	    "exon_clusters=s" => \$f_exon_cluster,
	    "exon_list_sp1=s" => \$f_exon_list_sp1,
	    "exon_list_sp2=s" => \$f_exon_list_sp2,
	    "sp1=s" => \$sp1,
	    "sp2=s" => \$sp2,
	    "dPSI_info=s" => \$dPSI_info,
	    "min_dPSI=i" => $min_dPSI,
	    "outfile=s" => \$outFile,
	    "help" => \$helpFlag
    );

### Help
if (!defined ($f_gene_cluster) || !defined($f_exon_cluster) || !defined($f_exon_list_sp1) || !defined($sp1) || !defined($sp2) ||  defined ($helpFlag)){
    die "
Usage: CompareExonsSets.pl -sp1 query_species -sp2 target_species -gene_clusters FILE -exon_clusters FILE -exon_list_sp1 FILE [OPTIONS]

Script to compare subsets of exons of interest between two species.

Compulsory Options:
     -sp1 query_species        Identifier of query species (it must match that used for the main module).
     -sp2 target_species       Identifier of target species (it must match that used for the main module).
     -gene_clusters FILE       File with clusters of gene orthology relationships (tsv).
     -exon_clusters FILE       File with clusters of exon orthology relationships (multi-species or pairwise [recommended])(tsv).
     -exon_list_sp1 FILE       File with exons from query species. It must contain qualitative or quantitative information on deltaPSIs.
                                  * Format (tsv): GeneID Exon_coord Info

Discretionary Options:
     -exon_list_sp2 FILE       File with exons from target species. It must contain qualitative or quantitative information on deltaPSIs.
                                  * If provided, the comparisons is done bidirectional (i.e. both species will be used as target and query).
     -dPSI_info auto           Type of information provided for each exon. Any of the following: dPSI, qual_call, none, auto.
                                  * none: no information is provided. All exons in the lists are treated as regulated.
                                  * dPSI: a value between -100 and 100. If a numeric value is provided (non NA), the exon is assumed to have sufficient coverage.
                                  * qual_call: a qualitative information. Valid values: UP, DOWN, NO_CHANGE, NO_COVERAGE (=NA or missing).
                                  * auto: decision among dPSI, qual_call and non will be done automatically [default].
     -min_dPSI int             Minimum absolute delta PSI used to make a qualitatitive UP or DOWN call if dPSI is provided as dPSI_info [def = 15].


*** Questions \& Bugs: mirimia\@gmail.com


";
}

###### Sanity checks
if ($dPSI_info ne "auto" && $dPSI_info ne "qual_call" && $dPSI_info ne "dPSI" && $dPSI_info ne "none"){
    die "[Aborted] dPSI_info must be dPSI, qual_call, none or auto\n";
}

## checks values
my %tally_info1=(); my $all_list1=0;
open (LIST1a, $f_exon_list_sp1) || die "It cannot open exon list for $sp1 ($f_exon_list_sp1)\n";
<LIST1a>; # not required to have a header, but just in case
while (<LIST1a>){
    chomp($_);
    my @t = split(/\t/,$_);
    my $temp_info="";

    if (defined $t[2]){
	$temp_info = $t[2];
	$temp_info="NUMERIC" if $temp_info=~/\d/;
    }
    else {
	$temp_info="EMPTY";
    }
    $tally_info1{$temp_info}++;
    $all_list1++;
}
close LIST1a;

$tally_info1{NUMERIC}=0 if !defined $tally_info1{NUMERIC}; $tally_info1{NA}=0 if !defined $tally_info1{NA};
$tally_info1{UP}=0 if !defined $tally_info1{UP}; $tally_info1{DOWN}=0 if !defined $tally_info1{DOWN};
$tally_info1{NO_CHANGE}=0 if !defined $tally_info1{NO_CHANGE}; $tally_info1{NO_COVERAGE}=0 if !defined $tally_info1{NO_COVERAGE};
$tally_info1{EMPTY}=0 if !defined $tally_info1{EMPTY};

if ($tally_info1{NUMERIC}+$tally_info1{NA} == $all_list1){
    $dPSI_info = "dPSI";
}
elsif ($tally_info1{UP}+$tally_info1{DOWN}+$tally_info1{NO_CHANGE}+$tally_info1{NO_COVERAGE}+$tally_info1{NA} == $all_list1){
    $dPSI_info = "qual_call";
}
elsif ($tally_info1{EMPTY} == $all_list1){
    $dPSI_info = "none";
}
else {
    die "[Aborted] The information provided does not meet the requirements\n".
	"[Aborted] dPSI: only numeric and NA values\n".
	"[Aborted] qual_call: only UP, DOWN, NO_CHANGE, NO_COVERAGE or NA values\n".
	"[Aborted] none: all values must be empty (no info column)\n";
}

### Check set 2 if available
if (defined $f_exon_list_sp2){
    my %tally_info2=(); my $all_list2=0;
    open (LIST2a, $f_exon_list_sp2) || die "It cannot open exon list for $sp2 ($f_exon_list_sp2)\n";
    <LIST2a>; # not required to have a header, but just in case
    while (<LIST2a>){
	chomp($_);
	my @t = split(/\t/,$_);
	my $temp_info="";
	
	if (defined $t[2]){
	    $temp_info = $t[2];
	    $temp_info="NUMERIC" if $temp_info=~/\d/;
	}
	else {
	    $temp_info="EMPTY";
	}
	$tally_info2{$temp_info}++;
	$all_list2++;
    }
    close LIST2a;
    
    $tally_info2{NUMERIC}=0 if !defined $tally_info2{NUMERIC}; $tally_info2{NA}=0 if !defined $tally_info2{NA};
    $tally_info2{UP}=0 if !defined $tally_info2{UP}; $tally_info2{DOWN}=0 if !defined $tally_info2{DOWN}; 
    $tally_info2{NO_CHANGE}=0 if !defined $tally_info2{NO_CHANGE}; $tally_info2{NO_COVERAGE}=0 if !defined $tally_info2{NO_COVERAGE};
    $tally_info2{EMPTY}=0 if !defined $tally_info2{EMPTY};
    
    if ($dPSI_info eq "dPSI"){
	die "[Aborted] Exon subset for $sp2 must have dPSI information\n" if $tally_info2{NUMERIC}+$tally_info2{NA} != $all_list2;
    }
    elsif ($dPSI_info eq "qual_call"){
	die "[Aborted] Exon subset for $sp2 must have qual_call information\n" if $tally_info2{UP}+$tally_info2{DOWN}+$tally_info2{NO_CHANGE}+$tally_info2{NO_COVERAGE}+$tally_info2{NA} != $all_list2;
    }
    elsif ($dPSI_info eq "none"){
	die "[Aborted] No info column should be provided\n" if $tally_info2{EMPTY} != $all_list2;
    }
}


# Format: GF0000002  Mmu  ENSMUSG00000069307
my %tally_genes=(); # hash to check for potentially incorrect species namings
my %gene_to_cluster=(); # conversion from geneID to clusterID
my %gene_cluster_has_sp=(); # exists if there is a gene for cluster X for Sp Y 

open (GCL, $f_gene_cluster) || die "It cannot open the gene orthology file ($f_gene_cluster)\n";
while (<GCL>){
    chomp($_);
    my @l=split(/\t/,$_);
    my $temp_sp = $l[1];
    if ($temp_sp eq $sp1 || $temp_sp eq $sp2) {
	$gene_to_cluster{$temp_sp}{$l[2]}=$l[0];
	$gene_cluster_has_sp{$l[0]}{$temp_sp}=1;
	$tally_genes{$l[1]}++;
    }
}
close GCL;

# Format: GF0000007.001  ENSMUSG00000096688  chr4:61595665-61595760:-  Mmu
my %tally_exons=(); # hash to check for potentially incorrect species namings
my %exon_to_cluster=(); # conversion from exonID to clusterID
my %exon_cluster_has_sp=(); # exists if there is an exon for cluster X for Sp Y 

open (ECL, $f_exon_cluster) || die "It cannot open the exon orthology file ($f_exon_cluster)\n";
while (<ECL>){
    chomp($_);
    my @l=split(/\t/,$_);
    my $temp_sp = $l[3];
    if ($temp_sp eq $sp1 || $temp_sp eq $sp2) {
	my ($chr,$i,$f)=$l[2]=~/(.+?)\:(\d+)\-(\d+)/;
	my $exon_id1 = "$l[1]=$i"; # geneID=exon_start
	my $exon_id2 = "$l[1]=$f"; # geneID=exon_end
	
	$exon_to_cluster{$temp_sp}{$exon_id1} = $l[0]; # must check for repetitions in the lists
	$exon_to_cluster{$temp_sp}{$exon_id2} = $l[0];
	$exon_cluster_has_sp{$l[0]}{$temp_sp}=1;
	$tally_exons{$temp_sp}++;
    }
}
close ECL;

### Checks species
die "[Abort] It could not find gene orthologs for $sp1; correct species 1 ID?\n" if !defined $tally_genes{$sp1};
die "[Abort] It could not find gene orthologs for $sp2; correct species 2 ID?\n" if !defined $tally_genes{$sp2};
die "[Abort] It could not find exon orthologs for $sp1; correct species 1 ID?\n" if !defined $tally_exons{$sp1};
die "[Abort] It could not find exon orthologs for $sp2; correct species 2 ID?\n" if !defined $tally_exons{$sp2};

### Options for last column: (i) dPSI, (ii) Qualitative calls, (iii) none
# ENSMUSG00000030279 chr6:143072543-143072694 DOWN/dPSI
my %exon_cluster_is_regulated=(); # ex_cluster sp
my %gene_cluster_is_regulated=(); # g_cluster sp
my $total_sp1_exons = 0;
my $total_sp1_genes = 0;
my $tally_sp1_exons_in_Gcons_genes = 0;
my $tally_sp1_genes_in_Gcons_genes = 0;
my $tally_sp1_exons_Gcons = 0;
my $tally_sp1_exons_in_Rcons_genes = 0;
my $tally_sp1_genes_in_Rcons_genes = 0;
my $tally_sp1_exons_Rcons = 0; # due to paralogy, it may be non-symmetrict to sp1

my %done_sp1_EX=(); # to store already seen PARTIAL exons
my %done_sp1_G=(); # count genes

#### Only useful if LIST2 is provided
my $total_sp2_exons = 0;
my $total_sp2_genes = 0;
my $tally_sp2_exons_in_Gcons_genes = 0;
my $tally_sp2_genes_in_Gcons_genes = 0;
my $tally_sp2_exons_Gcons = 0;
my $tally_sp2_exons_in_Rcons_genes = 0;
my $tally_sp2_genes_in_Rcons_genes = 0;
my $tally_sp2_exons_Rcons = 0; # due to paralogy, it may be non-symmetrict to sp1

my %done_sp2_EX=(); # to store already seen PARTIAL exons
my %done_sp2_G=(); # count genes


open (LIST1, $f_exon_list_sp1) || die "It cannot open exon list for $sp1 ($f_exon_list_sp1)\n";
while (<LIST1>){
    chomp($_);
    my @l=split(/\t/,$_);
    my $gene = $l[0];
    
    my ($i,$f)=$l[1]=~/\:(\d+)\-(\d+)/;
    my $exon_id1 = "$gene=$i"; # geneID=exon_start
    my $exon_id2 = "$gene=$f"; # geneID=exon_end
    
    ### to avoid multiple counting of exon variants (e.g. as in rMATS or SUPPA)
    next if ((defined $done_sp1_EX{$exon_id1}) || (defined $done_sp1_EX{$exon_id2}));
    $done_sp1_EX{$exon_id1}=1; $done_sp1_EX{$exon_id2}=1;
    $total_sp1_exons++;

    ##### Gene level conservation
    $total_sp1_genes++ unless (defined $done_sp1_G{$gene});
    ### Check first if the exon falls in an orthologous gene 
    if (defined $gene_to_cluster{$sp1}{$gene}){ # i.e. gene is in a cluster
	my $gene_cluster = $gene_to_cluster{$sp1}{$gene};
	$tally_sp1_exons_in_Gcons_genes++ if (defined $gene_cluster_has_sp{$gene_cluster}{$sp2});
	if (!defined $done_sp1_G{$gene}){
	    $tally_sp1_genes_in_Gcons_genes++ if (defined $gene_cluster_has_sp{$gene_cluster}{$sp2});
	}
	$gene_cluster_is_regulated{$gene_cluster}{$sp1}=1; 
    } # missing if Sp2 has a list => reg exons in the same genes
    $done_sp1_G{$gene}=1;

    if ((defined $exon_to_cluster{$sp1}{$exon_id1}) || (defined $exon_to_cluster{$sp1}{$exon_id2})){ 
	# it first checks the smaller coordinate. They could be inconsistent in a small number of cases.
	my $exon_cluster;
	$exon_cluster = $exon_to_cluster{$sp1}{$exon_id1} if (defined $exon_to_cluster{$sp1}{$exon_id1}); 
	$exon_cluster = $exon_to_cluster{$sp1}{$exon_id2} if (!defined $exon_to_cluster{$sp1}{$exon_id1}) && (defined $exon_to_cluster{$sp1}{$exon_id2}); 

	$tally_sp1_exons_Gcons++ if (defined $exon_cluster_has_sp{$exon_cluster}{$sp2});
	$exon_cluster_is_regulated{$exon_cluster}{$sp1}=1; 
    }
}
close LIST1;

if (defined $f_exon_list_sp2){
### Options for last column: (i) dPSI, (ii) Qualitative calls, (iii) none
# ENSMUSG00000030279 chr6:143072543-143072694 DOWN/dPSI
    open (LIST2, $f_exon_list_sp2) || die "It cannot open exon list for $sp2 ($f_exon_list_sp2)\n";
    while (<LIST2>){
	chomp($_);
	my @l=split(/\t/,$_);
	my $gene = $l[0];
	
	my ($i,$f)=$l[1]=~/\:(\d+)\-(\d+)/;
	my $exon_id1 = "$gene=$i"; # geneID=exon_start
	my $exon_id2 = "$gene=$f"; # geneID=exon_end
	
	### to avoid multiple counting of exon variants (e.g. as in rMATS or SUPPA)
	next if ((defined $done_sp2_EX{$exon_id1}) || (defined $done_sp2_EX{$exon_id2}));
	$done_sp2_EX{$exon_id1}=1; $done_sp2_EX{$exon_id2}=1;
	$total_sp2_exons++;

	##### Gene level conservation
	$total_sp2_genes++ unless (defined $done_sp2_G{$gene});
	### Check first if the exon falls in an orthologous gene 
	if (defined $gene_to_cluster{$sp2}{$gene}){ # i.e. gene is in a cluster
	    my $gene_cluster = $gene_to_cluster{$sp2}{$gene};
	    $tally_sp2_exons_in_Gcons_genes++ if (defined $gene_cluster_has_sp{$gene_cluster}{$sp1});
	    $tally_sp2_exons_in_Rcons_genes++ if (defined $gene_cluster_is_regulated{$gene_cluster}{$sp1});
	    
	    if (!defined $done_sp2_G{$gene}){
		$tally_sp2_genes_in_Gcons_genes++ if (defined $gene_cluster_has_sp{$gene_cluster}{$sp1});
		$tally_sp2_genes_in_Rcons_genes++ if (defined $gene_cluster_is_regulated{$gene_cluster}{$sp1});
	    }
	    $gene_cluster_is_regulated{$gene_cluster}{$sp2}=1; 
	} 	
	$done_sp2_G{$gene}=1;
	
	#### Exon levels conservation
	if ((defined $exon_to_cluster{$sp2}{$exon_id1}) || (defined $exon_to_cluster{$sp2}{$exon_id2})){ 
	    # it first checks the smaller coordinate. They could be inconsistent in a small number of cases.
	    my $exon_cluster;
	    $exon_cluster = $exon_to_cluster{$sp2}{$exon_id1} if (defined $exon_to_cluster{$sp2}{$exon_id1}); 
	    $exon_cluster = $exon_to_cluster{$sp2}{$exon_id2} if (!defined $exon_to_cluster{$sp2}{$exon_id1}) && (defined $exon_to_cluster{$sp2}{$exon_id2}); 
	    
	    $tally_sp2_exons_Gcons++ if (defined $exon_cluster_has_sp{$exon_cluster}{$sp1});
	    $tally_sp2_exons_Rcons++ if (defined $exon_cluster_is_regulated{$exon_cluster}{$sp1});
	    $exon_cluster_is_regulated{$exon_cluster}{$sp2}=1; 
	}
    }
    close LIST2;

    #### Re-open LIST1
    my %done_sp1_EX_b=(); my %done_sp1_G_b=();
    open (LIST1b, $f_exon_list_sp1) || die "It cannot open exon list for $sp1 ($f_exon_list_sp1)\n";
    while (<LIST1b>){
	chomp($_);
	my @l=split(/\t/,$_);
	my $gene = $l[0];
	
	my ($i,$f)=$l[1]=~/\:(\d+)\-(\d+)/;
	my $exon_id1 = "$gene=$i"; # geneID=exon_start
	my $exon_id2 = "$gene=$f"; # geneID=exon_end
	
	### to avoid multiple counting of exon variants (e.g. as in rMATS or SUPPA)
	next if ((defined $done_sp1_EX_b{$exon_id1}) || (defined $done_sp1_EX_b{$exon_id2}));
	$done_sp1_EX_b{$exon_id1}=1; $done_sp1_EX_b{$exon_id2}=1;
	
	### Gene level conservation
	if (defined $gene_to_cluster{$sp1}{$gene}){ # i.e. gene is in a cluster
	    my $gene_cluster = $gene_to_cluster{$sp1}{$gene};
	    $tally_sp1_exons_in_Rcons_genes++ if (defined $gene_cluster_is_regulated{$gene_cluster}{$sp2});
	    if (!defined $done_sp1_G_b{$gene}){
		$tally_sp1_genes_in_Rcons_genes++ if (defined $gene_cluster_is_regulated{$gene_cluster}{$sp2});
	    }
	}
	$done_sp1_G_b{$gene}=1;

	### Exon level
	if ((defined $exon_to_cluster{$sp1}{$exon_id1}) || (defined $exon_to_cluster{$sp1}{$exon_id2})){ 
	    # it first checks the smaller coordinate. They could be inconsistent in a small number of cases.
	    my $exon_cluster;
	    $exon_cluster = $exon_to_cluster{$sp1}{$exon_id1} if (defined $exon_to_cluster{$sp1}{$exon_id1}); 
	    $exon_cluster = $exon_to_cluster{$sp1}{$exon_id2} if (!defined $exon_to_cluster{$sp1}{$exon_id1}) && (defined $exon_to_cluster{$sp1}{$exon_id2}); 

	    $tally_sp1_exons_Rcons++ if (defined $exon_cluster_is_regulated{$exon_cluster}{$sp2});
	}
    }
    close LIST1b;
}

### Prints basic stats
if (defined $f_exon_list_sp2){
    ### gene-level
    # perc of genes with an ortholog
    my $perc_sp1_genes_Gcons_genes = sprintf ("%.2f", 100*$tally_sp1_genes_in_Gcons_genes/$total_sp1_genes);
    my $perc_sp2_genes_Gcons_genes = sprintf ("%.2f", 100*$tally_sp2_genes_in_Gcons_genes/$total_sp2_genes);
    # perc of genes with an ortholog with a R-cons exon
    my $perc_sp1_genes_Rcons_genes = sprintf ("%.2f", 100*$tally_sp1_genes_in_Rcons_genes/$total_sp1_genes);
    my $perc_sp2_genes_Rcons_genes = sprintf ("%.2f", 100*$tally_sp2_genes_in_Rcons_genes/$total_sp2_genes);

    ### exon-level
    # perc of exons with a gene ortholog
    my $perc_sp1_exons_Gcons_genes = sprintf ("%.2f", 100*$tally_sp1_exons_in_Gcons_genes/$total_sp1_exons);
    my $perc_sp2_exons_Gcons_genes = sprintf ("%.2f", 100*$tally_sp2_exons_in_Gcons_genes/$total_sp2_exons);
    # perc of exons whose gene orthologs have regulated exons
    my $perc_sp1_exons_Rcons_genes = sprintf ("%.2f", 100*$tally_sp1_exons_in_Rcons_genes/$total_sp1_exons);
    my $perc_sp2_exons_Rcons_genes = sprintf ("%.2f", 100*$tally_sp2_exons_in_Rcons_genes/$total_sp2_exons);
    # perc of exons with exon ortholog (Gcons)
    my $perc_sp1_exons_Gcons_exons = sprintf ("%.2f", 100*$tally_sp1_exons_Gcons/$total_sp1_exons);    
    my $perc_sp2_exons_Gcons_exons = sprintf ("%.2f", 100*$tally_sp2_exons_Gcons/$total_sp2_exons);    
    # perc of Gcons exons out of the total N of exons with gene orthologs
    my $perc_sp1_exons_Gcons_exons_OrthGenes = sprintf ("%.2f", 100*$tally_sp1_exons_Gcons/$tally_sp1_exons_in_Gcons_genes);
    my $perc_sp2_exons_Gcons_exons_OrthGenes = sprintf ("%.2f", 100*$tally_sp2_exons_Gcons/$tally_sp2_exons_in_Gcons_genes);
    # perc of Rcons exons
    my $perc_sp1_exons_Rcons_exons = sprintf ("%.2f", 100*$tally_sp1_exons_Rcons/$total_sp1_exons);
    my $perc_sp2_exons_Rcons_exons = sprintf ("%.2f", 100*$tally_sp2_exons_Rcons/$total_sp2_exons);
    # perc of Rcons exons out of the total N of exons with gene orthologs      
    my $perc_sp1_exons_Rcons_exons_OrthGenes = sprintf ("%.2f", 100*$tally_sp1_exons_Rcons/$tally_sp1_exons_in_Gcons_genes);
    my $perc_sp2_exons_Rcons_exons_OrthGenes = sprintf ("%.2f", 100*$tally_sp2_exons_Rcons/$tally_sp2_exons_in_Gcons_genes);
    # perc of Rcons exons out of the total Gcons exons
    my $perc_sp1_exons_Rcons_exons_Gcons = sprintf ("%.2f", 100*$tally_sp1_exons_Rcons/$tally_sp1_exons_Gcons) if $tally_sp1_exons_Gcons>0;
    my $perc_sp2_exons_Rcons_exons_Gcons = sprintf ("%.2f", 100*$tally_sp2_exons_Rcons/$tally_sp2_exons_Gcons) if $tally_sp2_exons_Gcons>0;

    print "
- Gene-level stats:
   - $sp1 => $sp2
Genes with exons from $sp1 in the exon lists\t$total_sp1_genes
Genes with exons from $sp1 with gene orthologs in $sp2\t$tally_sp1_genes_in_Gcons_genes\t$perc_sp1_genes_Gcons_genes\%
Genes with exons from $sp1 with gene orthologs with regulated exons in $sp2\t$tally_sp1_genes_in_Rcons_genes\t$perc_sp1_genes_Rcons_genes\%

   - $sp2 => $sp1
Genes with exons from $sp2 in the exon lists\t$total_sp2_genes
Genes with exons from $sp2 with gene orthologs in $sp1\t$tally_sp2_genes_in_Gcons_genes\t$perc_sp2_genes_Gcons_genes\%
Genes with exons from $sp2 with gene orthologs with regulated exons in $sp1\t$tally_sp2_genes_in_Rcons_genes\t$perc_sp2_genes_Rcons_genes\%


- Exon-level stats:
   - $sp1 => $sp2
Exons from $sp1 in exon list\t$total_sp1_exons
Exons from $sp1 with gene orthologs in $sp2\t$tally_sp1_exons_in_Gcons_genes\t$perc_sp1_exons_Gcons_genes\%
Exons from $sp1 with gene orthologs with regulated exons in $sp2\t$tally_sp1_exons_in_Rcons_genes\t$perc_sp1_exons_Rcons_genes\%
Exons from $sp1 with exon orthologs in $sp2 (G-conserved)\t$tally_sp1_exons_Gcons\t$perc_sp1_exons_Gcons_exons\%\t$perc_sp1_exons_Gcons_exons_OrthGenes\% (in Orth genes)
Exons from $sp1 with regulated exon orthologs in $sp2 (R-conserved)\t$tally_sp1_exons_Rcons\t$perc_sp1_exons_Rcons_exons\%\t$perc_sp1_exons_Rcons_exons_OrthGenes\% (in Orth genes)
    Percent of R-conserved exons from $sp1 out of G-conserved exons\t$perc_sp1_exons_Rcons_exons_Gcons\%

   - $sp2 => $sp1
Exons from $sp2 in exon list\t$total_sp2_exons
Exons from $sp2 with gene orthologs in $sp1\t$tally_sp2_exons_in_Gcons_genes\t$perc_sp2_exons_Gcons_genes\%
Exons from $sp2 with gene orthologs with regulated exons in $sp1\t$tally_sp2_exons_in_Rcons_genes\t$perc_sp2_exons_Rcons_genes\%
Exons from $sp2 with exon orthologs in $sp1 (G-conserved)\t$tally_sp2_exons_Gcons\t$perc_sp2_exons_Gcons_exons\%\t$perc_sp2_exons_Gcons_exons_OrthGenes\% (in Orth genes)
Exons from $sp2 with regulated exon orthologs in $sp1 (R-conserved)\t$tally_sp2_exons_Rcons\t$perc_sp2_exons_Rcons_exons\%\t$perc_sp2_exons_Rcons_exons_OrthGenes\% (in Orth genes)
    Percent of R-conserved exons from $sp2 out of G-conserved exons\t$perc_sp2_exons_Rcons_exons_Gcons\%


";
# it could do the "unique cluster" counts for exon and genes
}
else {
    ### gene-level
    # perc of genes with an ortholog
    my $perc_sp1_genes_Gcons_genes = sprintf ("%.2f", 100*$tally_sp1_genes_in_Gcons_genes/$total_sp1_genes);
    # perc of genes with an ortholog with a R-cons exon
    my $perc_sp1_genes_Rcons_genes = "NA";

    ### exon-level
    # perc of exons with a gene ortholog
    my $perc_sp1_exons_Gcons_genes = sprintf ("%.2f", 100*$tally_sp1_exons_in_Gcons_genes/$total_sp1_exons);
    # perc of exons whose gene orthologs have regulated exons
    my $perc_sp1_exons_Rcons_genes = "NA";
    # perc of exons with exon ortholog (Gcons)
    my $perc_sp1_exons_Gcons_exons = sprintf ("%.2f", 100*$tally_sp1_exons_Gcons/$total_sp1_exons);    
    # perc of Gcons exons out of the total N of exons with gene orthologs
    my $perc_sp1_exons_Gcons_exons_OrthGenes = sprintf ("%.2f", 100*$tally_sp1_exons_Gcons/$tally_sp1_exons_in_Gcons_genes);
    # perc of Rcons exons
    my $perc_sp1_exons_Rcons_exons = "NA";
    # perc of Rcons exons out of the total N of exons with gene orthologs      
    my $perc_sp1_exons_Rcons_exons_OrthGenes = "NA";
    # perc of Rcons exons out of the total Gcons exons
    my $perc_sp1_exons_Rcons_exons_Gcons = "NA";

    # NA's some variables
    $tally_sp1_genes_in_Rcons_genes="NA";
    $tally_sp1_exons_in_Rcons_genes="NA";
    $tally_sp1_exons_Rcons="NA";    

    print "
- Gene-level stats ($sp1 => $sp2):
Genes with exons from $sp1 in the exon lists\t$total_sp1_genes
Genes with exons from $sp1 with gene orthologs in $sp2\t$tally_sp1_genes_in_Gcons_genes\t$perc_sp1_genes_Gcons_genes\%
Genes with exons from $sp1 with gene orthologs with regulated exons in $sp2\t$tally_sp1_genes_in_Rcons_genes\t$perc_sp1_genes_Rcons_genes\%

- Exon-level stats ($sp1 => $sp2):
Exons from $sp1 in exon list\t$total_sp1_exons
Exons from $sp1 with gene orthologs in $sp2\t$tally_sp1_exons_in_Gcons_genes\t$perc_sp1_exons_Gcons_genes\%
Exons from $sp1 with gene orthologs with regulated exons in $sp2\t$tally_sp1_exons_in_Rcons_genes\t$perc_sp1_exons_Rcons_genes\%
Exons from $sp1 with exon orthologs in $sp2 (G-conserved)\t$tally_sp1_exons_Gcons\t$perc_sp1_exons_Gcons_exons\%\t$perc_sp1_exons_Gcons_exons_OrthGenes\% (in Orth genes)
Exons from $sp1 with regulated exon orthologs in $sp2 (R-conserved)\t$tally_sp1_exons_Rcons\t$perc_sp1_exons_Rcons_exons\%\t$perc_sp1_exons_Rcons_exons_OrthGenes\% (in Orth genes)


";
# it could do the "unique cluster" counts for exon and genes
}
