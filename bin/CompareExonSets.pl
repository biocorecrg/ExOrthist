#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

# thoughts: 
# - upload exons_overlap to have ALL exons.
# - upload the file with the best hits? 
#     - check 1) where the best is and how far it falls from the targt
#     - If it’s the best or close: overlapping/orth seq.

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
my $max_dif_le_ratio = 1.5; # max (le1-le2)/le1 or (le1-le2)/le2 allowed

Getopt::Long::Configure("no_auto_abbrev");
GetOptions( "gene_clusters=s" => \$f_gene_cluster,
	    "exon_clusters=s" => \$f_exon_cluster,
	    "exon_list_sp1=s" => \$f_exon_list_sp1,
	    "exon_list_sp2=s" => \$f_exon_list_sp2,
	    "sp1=s" => \$sp1,
	    "sp2=s" => \$sp2,
	    "dPSI_info=s" => \$dPSI_info,
	    "min_dPSI=i" => $min_dPSI,
	    "outFile" => \$outFile,
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
                                  * qual_call: a qualitative information. Valid values: UP, DOWN, REGULATED, NO_CHANGE, NO_COVERAGE (=NA or missing).
                                  * auto: decision among dPSI, qual_call and non will be done automatically [default].
     -min_dPSI int             Minimum absolute delta PSI used to make a qualitatitive UP or DOWN call if dPSI is provided as dPSI_info [def = 15].
     -outFile                  It creates an output file with the exons in conserved clusters (otherwise, it does NOT create it).


*** Questions \& Bugs: mirimia\@gmail.com


";
}

# Changes to implement:
# - load the lists from the begining. Loop through hash instead of list
# - load the clusters first, and then the lists (in case redundancy is found).
# - convert the info already: dPSI => X. none => REGULATED. Add into %info_by_exon.
# - if ALL exons are provided => add all exons to the gene string: more accurate
# - only test exons that are REG UP DOWN.


########### Sanity checks and detection of dPSI_info format
if ($dPSI_info ne "auto" && $dPSI_info ne "qual_call" && $dPSI_info ne "dPSI" && $dPSI_info ne "none"){
    die "[Aborted] dPSI_info must be dPSI, qual_call, none or auto\n";
}
## checks values automatically
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
############################ Finish detecting format of dPSI



###### Starts loading data
### Gene clusters
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

### Exons clusters
# Format: GF0000007.001  ENSMUSG00000096688  chr4:61595665-61595760:-  Mmu
my %tally_exons=(); # hash to check for potentially incorrect species namings
my %exon_to_cluster=(); # conversion from exonID to clusterID
my %exon_cluster_has_sp=(); # exists if there is an exon for cluster X for Sp Y 
my %exons_by_cluster=(); # all the exons in a cluster organized by species
my %partial_coord_conversion=(); # to match full coordinates of clusters and lists
my %array_of_exons_per_gene=(); # an array of coordinates per gene, to derive a first non-cons score
my %gene_strand=(); # keeps the strand of each gene (currently, only if it has exon clusters

open (ECL, $f_exon_cluster) || die "It cannot open the exon orthology file ($f_exon_cluster)\n";
while (<ECL>){
    chomp($_);
    my @l=split(/\t/,$_);
    my $exon_cluster = $l[0];
    my $gene = $l[1];
    my $temp_sp = $l[3];
    if ($temp_sp eq $sp1 || $temp_sp eq $sp2) {
	my ($chr,$i,$f,$strand)=$l[2]=~/(.+?)\:(\d+)\-(\d+)\:(.)/;
	my $exon_id1 = "$l[1]=$i"; # geneID=exon_start
	my $exon_id2 = "$l[1]=$f"; # geneID=exon_end
	my $exon_id = "$l[1]=$chr:$i-$f";
	
	$gene_strand{$temp_sp}{$gene}=$strand;
	$tally_exons{$temp_sp}++; # sanity check for species provided
	$exon_to_cluster{$temp_sp}{$exon_id1} = $exon_cluster; # to check for variants in the list
	$exon_to_cluster{$temp_sp}{$exon_id2} = $exon_cluster; # to check for variants in the list
	$exon_cluster_has_sp{$exon_cluster}{$temp_sp}=1;
	push(@{$exons_by_cluster{$exon_cluster}{$temp_sp}},$exon_id); # array of all exons per cluster for each species
	$partial_coord_conversion{$temp_sp}{$exon_id1}=$exon_id; # to avoid inconsistencies between full coord in list and clusters
	$partial_coord_conversion{$temp_sp}{$exon_id2}=$exon_id; # to avoid inconsistencies between full coord in list and clusters
	push(@{$array_of_exons_per_gene{$temp_sp}{$gene}},"$i-$f=$exon_cluster"); # all coordinates of the exons per gene
    }
}
close ECL;

### Checks species
die "[Abort] It could not find gene orthologs for $sp1; correct species 1 ID?\n" if !defined $tally_genes{$sp1};
die "[Abort] It could not find gene orthologs for $sp2; correct species 2 ID?\n" if !defined $tally_genes{$sp2};
die "[Abort] It could not find exon orthologs for $sp1; correct species 1 ID?\n" if !defined $tally_exons{$sp1};
die "[Abort] It could not find exon orthologs for $sp2; correct species 2 ID?\n" if !defined $tally_exons{$sp2};

### Genes with regulated exons
if (defined $outFile && defined $f_exon_list_sp2){
    open (OUT_GENES, ">OrthoGenes_with_reg_exons-$sp1-$sp2.tab");
    print OUT_GENES "OrthoCluster\tGeneID_1\tExon_coord_1\tExon_length_1\tEx_cluster_1\t".
	"GeneID_2\tExon_coord_2\tExon_length_2\tEx_cluster_2\tCONSERV_CALL\tSpecies_1\tSpecies_2\n";
}

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
my %tally_sp1_exons_in_Rcons_genes_by_type=(); # when the exon is in a conserved gene with a reg exon => by type of conservation

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

my %info_by_exon=(); # all necessary info
my %exact_exon_with_a_match=(); # this stores for each cluster all the exons matched in lists [the coord may not be exact]

open (LIST1, $f_exon_list_sp1) || die "It cannot open exon list for $sp1 ($f_exon_list_sp1)\n";
while (<LIST1>){
    chomp($_);
    my @l=split(/\t/,$_);
    my $gene = $l[0];
    
    my ($i,$f)=$l[1]=~/\:(\d+)\-(\d+)/;
    my $exon_id1 = "$gene=$i"; # geneID=exon_start
    my $exon_id2 = "$gene=$f"; # geneID=exon_end
    my $exon_id = "$gene=$l[1]"; # fullID: geneID=coord

    ### to avoid multiple counting of exon variants (e.g. as in rMATS or SUPPA)
    next if ((defined $done_sp1_EX{$exon_id1}) || (defined $done_sp1_EX{$exon_id2}));
    $done_sp1_EX{$exon_id1}=1; $done_sp1_EX{$exon_id2}=1;
    $total_sp1_exons++;

    ### tries to assign to ex cluster
    my $t_exon_cl1; # redundant with some code below
    if (defined $exon_to_cluster{$sp1}{$exon_id1}){$t_exon_cl1=$exon_to_cluster{$sp1}{$exon_id1};}
    elsif (defined $exon_to_cluster{$sp1}{$exon_id2}){$t_exon_cl1=$exon_to_cluster{$sp1}{$exon_id2};}
    else {$t_exon_cl1 = "NO_CLUSTER";}

    ##### Gene level conservation
    $total_sp1_genes++ unless (defined $done_sp1_G{$gene});
    ### Check first if the exon falls in an orthologous gene 
    if (defined $gene_to_cluster{$sp1}{$gene}){ # i.e. gene is in a cluster
	my $gene_cluster = $gene_to_cluster{$sp1}{$gene};
	$tally_sp1_exons_in_Gcons_genes++ if (defined $gene_cluster_has_sp{$gene_cluster}{$sp2});
	if (!defined $done_sp1_G{$gene}){
	    $tally_sp1_genes_in_Gcons_genes++ if (defined $gene_cluster_has_sp{$gene_cluster}{$sp2});
	}

	if (!defined $gene_cluster_is_regulated{$gene_cluster}{$sp1}){
	    $gene_cluster_is_regulated{$gene_cluster}{$sp1}.="$exon_id=$t_exon_cl1,";
	}
	elsif ($gene_cluster_is_regulated{$gene_cluster}{$sp1}!~/$exon_id\=$t_exon_cl1\,/){
	    $gene_cluster_is_regulated{$gene_cluster}{$sp1}.="$exon_id=$t_exon_cl1,";
	} 	
    } # missing if Sp2 has a list => reg exons in the same genes
    $done_sp1_G{$gene}=1;

    if ((defined $exon_to_cluster{$sp1}{$exon_id1}) || (defined $exon_to_cluster{$sp1}{$exon_id2})){ 
	# it first checks the smaller coordinate. They could be inconsistent in a small number of cases.
	my $exon_cluster;
	if (defined $exon_to_cluster{$sp1}{$exon_id1}){
	    $exon_cluster = $exon_to_cluster{$sp1}{$exon_id1};
	    $exon_id = $partial_coord_conversion{$sp1}{$exon_id1}; # redefines exon_id!
	}
	elsif (defined $exon_to_cluster{$sp1}{$exon_id2}){
	    $exon_cluster = $exon_to_cluster{$sp1}{$exon_id2};
	    $exon_id = $partial_coord_conversion{$sp1}{$exon_id2}; # redefines exon_id!
	}

	$tally_sp1_exons_Gcons++ if (defined $exon_cluster_has_sp{$exon_cluster}{$sp2});
	$exon_cluster_is_regulated{$exon_cluster}{$sp1}=1; 
	$exact_exon_with_a_match{$exon_cluster}{$exon_id}=1;
    }
    else { # added to the gene's exon array as non Gcons
	push(@{$array_of_exons_per_gene{$sp1}{$gene}},"$i-$f=NO_CLUSTER"); # all coordinates of the exons per gene
    }
    
    #### Defines the info here, in case exon_id is redefined by cluster match
    if ($dPSI_info eq "none"){
	$info_by_exon{$exon_id}="REGULATED";
    }
    elsif ($dPSI_info eq "qual_call"){
	$info_by_exon{$exon_id}="$l[2]";
    }
    # to be developed further
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
	my $exon_id = "$gene=$l[1]"; # fullID: geneID=coord
	
	### to avoid multiple counting of exon variants (e.g. as in rMATS or SUPPA)
	next if ((defined $done_sp2_EX{$exon_id1}) || (defined $done_sp2_EX{$exon_id2}));
	$done_sp2_EX{$exon_id1}=1; $done_sp2_EX{$exon_id2}=1;
	$total_sp2_exons++;

	### tries to assign to ex cluster
	my $t_exon_cl2; # redundant with some code below
	if (defined $exon_to_cluster{$sp2}{$exon_id1}){$t_exon_cl2=$exon_to_cluster{$sp2}{$exon_id1};}
	elsif (defined $exon_to_cluster{$sp2}{$exon_id2}){$t_exon_cl2=$exon_to_cluster{$sp2}{$exon_id2};}
	else {$t_exon_cl2 = "NO_CLUSTER";}

	##### Gene level conservation
	$total_sp2_genes++ unless (defined $done_sp2_G{$gene});
	### Check first if the exon falls in an orthologous gene 
	if (defined $gene_to_cluster{$sp2}{$gene}){ # i.e. gene is in a cluster
	    my $gene_cluster = $gene_to_cluster{$sp2}{$gene};
	    $tally_sp2_exons_in_Gcons_genes++ if (defined $gene_cluster_has_sp{$gene_cluster}{$sp1});

	    if (defined $gene_cluster_is_regulated{$gene_cluster}{$sp1}){ # gene orth is Sp1 have regulated exons (cons or convergent?)
		$tally_sp2_exons_in_Rcons_genes++;
		
		### This bit is to print the genes with Reg events (it's redundant with the reciprocal Sp1 => Sp2)
#                $gene_cluster_is_regulated{$gene_cluster}{$sp1}=~s/\,$//;
#                my ($t_g2,$t_co2)=$exon_id=~/(.+?)\=(.+)/;
#                my @t_ex_sp1=split(/\,/,$gene_cluster_is_regulated{$gene_cluster}{$sp1});
#                foreach my $t_exon_id1 (@t_ex_sp1){
#                    my ($t_g1,$t_co1,$t_exon_cl1)=$t_exon_id1=~/(.+?)\=(.+)\=(.+)/;
#                    print OUT_GENES "$gene_cluster\t$t_g2\t$t_co2\t$t_exon_cl2\t$t_g1\t$t_co1\t$t_exon_cl1\t$sp2\t$sp1\n";
#                }
	    }

	    if (!defined $done_sp2_G{$gene}){
		$tally_sp2_genes_in_Gcons_genes++ if (defined $gene_cluster_has_sp{$gene_cluster}{$sp1});
		$tally_sp2_genes_in_Rcons_genes++ if (defined $gene_cluster_is_regulated{$gene_cluster}{$sp1});
	    }
	    
	    if (!defined $gene_cluster_is_regulated{$gene_cluster}{$sp2}){
		$gene_cluster_is_regulated{$gene_cluster}{$sp2}.="$exon_id=$t_exon_cl2,";
	    }
	    elsif ($gene_cluster_is_regulated{$gene_cluster}{$sp2}!~/$exon_id\=$t_exon_cl2\,/){
		$gene_cluster_is_regulated{$gene_cluster}{$sp2}.="$exon_id=$t_exon_cl2,";
	    } 	
	}
	$done_sp2_G{$gene}=1;
	
	#### Exon levels conservation
	if ((defined $exon_to_cluster{$sp2}{$exon_id1}) || (defined $exon_to_cluster{$sp2}{$exon_id2})){ 
	    # it first checks the smaller coordinate. They could be inconsistent in a small number of cases.
	    my $exon_cluster;

	    if (defined $exon_to_cluster{$sp2}{$exon_id1}){
		$exon_cluster = $exon_to_cluster{$sp2}{$exon_id1};
		$exon_id = $partial_coord_conversion{$sp2}{$exon_id1}; # redefines exon_id!
	    }
	    elsif (defined $exon_to_cluster{$sp2}{$exon_id2}){
		$exon_cluster = $exon_to_cluster{$sp2}{$exon_id2};
		$exon_id = $partial_coord_conversion{$sp2}{$exon_id2}; # redefines exon_id!
	    }
	    
	    $tally_sp2_exons_Gcons++ if (defined $exon_cluster_has_sp{$exon_cluster}{$sp1});
	    $tally_sp2_exons_Rcons++ if (defined $exon_cluster_is_regulated{$exon_cluster}{$sp1});
	    $exon_cluster_is_regulated{$exon_cluster}{$sp2}=1; 
	    $exact_exon_with_a_match{$exon_cluster}{$exon_id}=1;
	}
	else { # added to the gene's exon array as non Gcons
	    push(@{$array_of_exons_per_gene{$sp2}{$gene}},"$i-$f=NO_CLUSTER"); # all coordinates of the exons per gene
	}

	#### Defines the info here, in case exon_id is redefined by cluster match
	if ($dPSI_info eq "none"){
	    $info_by_exon{$exon_id}="REGULATED";
	}
	elsif ($dPSI_info eq "qual_call"){
	    $info_by_exon{$exon_id}="$l[2]";
	}
	# to be developed further
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
	my $exon_id = "$gene=$l[1]"; # fullID: geneID=coord
	
	### to avoid multiple counting of exon variants (e.g. as in rMATS or SUPPA)
	next if ((defined $done_sp1_EX_b{$exon_id1}) || (defined $done_sp1_EX_b{$exon_id2}));
	$done_sp1_EX_b{$exon_id1}=1; $done_sp1_EX_b{$exon_id2}=1;

	### tries to assign to ex cluster
	my $t_exon_cl1; # redundant with some code below
	if (defined $exon_to_cluster{$sp1}{$exon_id1}){$t_exon_cl1=$exon_to_cluster{$sp1}{$exon_id1};}
	elsif (defined $exon_to_cluster{$sp1}{$exon_id2}){$t_exon_cl1=$exon_to_cluster{$sp1}{$exon_id2};}
	else {$t_exon_cl1 = "NO_CLUSTER";}
	
	### Gene level conservation
	if (defined $gene_to_cluster{$sp1}{$gene}){ # i.e. gene is in a cluster
	    my $gene_cluster = $gene_to_cluster{$sp1}{$gene};
	    if (defined $gene_cluster_is_regulated{$gene_cluster}{$sp2}){ # gene orth have regulated exon Sp2 (cons or convergent?)
		$tally_sp1_exons_in_Rcons_genes++;

		### This bit is to print the genes with Reg events
		$gene_cluster_is_regulated{$gene_cluster}{$sp2}=~s/\,$//; # gets rid of final comma
		my ($t_g1,$t_co1)=$exon_id=~/(.+?)\=(.+)/;
		my @t_ex_sp2=split(/\,/,$gene_cluster_is_regulated{$gene_cluster}{$sp2});
		foreach my $t_exon_id2 (@t_ex_sp2){ # all exons regulated in gene orthologs from Sp2
		    my ($t_g2,$t_co2,$t_exon_cl2)=$t_exon_id2=~/(.+?)\=(.+?)\=(.+)/;
		    my ($i2,$f2) = $t_co2 =~ /\:(\d+)\-(\d+)/;
		    ### Gets exons lengths
		    my $le_ex1 = ($f-$i)+1;
		    my $le_ex2 = ($f2-$i2)+1;

		    # Tries to check conservation vs convergence		    
		    my $cons_conv_call = "UNCLEAR_A";
		    if ($t_exon_cl1 eq $t_exon_cl2 && $t_exon_cl1 ne "NO_CLUSTER" && $t_exon_cl2 ne "NO_CLUSTER"){
			$cons_conv_call = "CONSERVED";
		    }
		    elsif ((abs($le_ex1-$le_ex2)/$le_ex1) > $max_dif_le_ratio || (abs($le_ex1-$le_ex2)/$le_ex2) > $max_dif_le_ratio){
			$cons_conv_call = "NON_CONSERVED_L";
		    }
		    elsif ($t_exon_cl1 ne $t_exon_cl2 && $t_exon_cl1 ne "NO_CLUSTER" && $t_exon_cl2 ne "NO_CLUSTER"){
			$cons_conv_call= "NON_CONSERVED_A";
		    }
		    elsif ($t_exon_cl1 ne "NO_CLUSTER" && $exon_cluster_has_sp{$t_exon_cl1}{$sp2}){ # there is another Sp2 exon in the Sp1 cluster
			$cons_conv_call= "NON_CONSERVED_B";
		    }
		    elsif ($t_exon_cl2 ne "NO_CLUSTER" && $exon_cluster_has_sp{$t_exon_cl2}{$sp1}){ # there is another Sp1 exon in the Sp2 cluster
			$cons_conv_call= "NON_CONSERVED_C";
		    }
		    elsif (!defined $gene_strand{$sp1}{$t_g1} || !defined $gene_strand{$sp2}{$t_g2}){ # either of them has no clusters
			$cons_conv_call = "UNCLEAR_B";
		    }
		    else { # both genes have exon clusters
			my @temp_sp1; my @temp_sp2;
			@temp_sp1 = sort {($a=~/(\d+)\-/)[0]<=>($b=~/(\d+)\-/)[0]}(@{$array_of_exons_per_gene{$sp1}{$t_g1}}) if $gene_strand{$sp1}{$t_g1} eq "+"; 
			@temp_sp1 = sort {($b=~/(\d+)\-/)[0]<=>($a=~/(\d+)\-/)[0]}(@{$array_of_exons_per_gene{$sp1}{$t_g1}}) if $gene_strand{$sp1}{$t_g1} eq "-"; 
			@temp_sp2 = sort {($a=~/(\d+)\-/)[0]<=>($b=~/(\d+)\-/)[0]}(@{$array_of_exons_per_gene{$sp2}{$t_g2}}) if $gene_strand{$sp2}{$t_g2} eq "+"; 
			@temp_sp2 = sort {($b=~/(\d+)\-/)[0]<=>($a=~/(\d+)\-/)[0]}(@{$array_of_exons_per_gene{$sp2}{$t_g2}}) if $gene_strand{$sp2}{$t_g2} eq "-"; 

			my $match1; my $match2; # the position of the actual exons
			my %cl_to_index1; my %cl_to_index2; # from the cluster ID to index in array
			my $anchor_ups; my $anchor_downs;
			my %anchor_conversion;
			foreach my $t_ex (0..$#temp_sp1){
			    $match1 = $t_ex if ($temp_sp1[$t_ex]=~/$i/ || $temp_sp1[$t_ex]=~/$f/);
			    my ($value) = $temp_sp1[$t_ex] =~ /\=(.+)/;
			    $cl_to_index1{$value}=$t_ex;
			}
			foreach my $t_ex (0..$#temp_sp2){
			    $match2 = $t_ex if ($temp_sp2[$t_ex]=~/$i2/ || $temp_sp2[$t_ex]=~/$f2/);
			    my ($value) = $temp_sp2[$t_ex] =~ /\=(.+)/;
			    $cl_to_index2{$value}=$t_ex;
			    if (defined ($cl_to_index1{$value})){ # it can take 0
				if ($value ne "NO_CLUSTER"){ #i.e. it's a cluster = anchor
				    $anchor_conversion{$cl_to_index1{$value}}=$t_ex;  # Sp1 => Sp2: e.g 23 => 2, 24 => 4, etc
				    if (!defined $anchor_ups){
					$anchor_ups = $cl_to_index1{$value} if $cl_to_index1{$value} < $match1; 
				    } else {
					$anchor_ups = $cl_to_index1{$value} if $cl_to_index1{$value} < $match1 && $cl_to_index1{$value} > $anchor_ups; 
				    }
				    if (!defined $anchor_downs){
					$anchor_downs = $cl_to_index1{$value} if $cl_to_index1{$value} > $match1; 
				    } else {
					$anchor_downs = $cl_to_index1{$value} if $cl_to_index1{$value} > $match1 && $cl_to_index1{$value} < $anchor_downs; 
				    }
				}
			    }
			}
			### does the anchors assessment
			if (!defined $anchor_ups && !defined $anchor_downs){ # no common clusters
			    $cons_conv_call= "UNCLEAR_C";
			}
			elsif (defined $anchor_ups && !defined $anchor_downs){ # anchor upstream
			    $cons_conv_call= "UNCLEAR_D" if $match2 > $anchor_conversion{$anchor_ups}; # it's downstream the anchor
			    $cons_conv_call= "NON_CONSERVED_D" if $match2 <= $anchor_conversion{$anchor_ups}; # it's upstream
			}
			elsif (!defined $anchor_ups && defined $anchor_downs){ # anchor downstream
			    $cons_conv_call= "UNCLEAR_E" if $match2 < $anchor_conversion{$anchor_downs}; # it's upstream the anchor
			    $cons_conv_call= "NON_CONSERVED_E" if $match2 >= $anchor_conversion{$anchor_downs}; # it's downstream
			}
			else { # both anchors
			    if ($match2 > $anchor_conversion{$anchor_ups} && $match2 < $anchor_conversion{$anchor_downs}){
				$cons_conv_call= "UNCLEAR_F";
			    }
			    else {
				$cons_conv_call= "NON_CONSERVED_F";
			    }
			}
##### Testing and debugging
#			my $array_sp1 = join(" ", @temp_sp1); my $array_sp2 = join(" ", @temp_sp2);
#			$anchor_ups="NA" if !defined $anchor_ups;
#			$anchor_downs="NA" if !defined $anchor_downs;
#			$anchor_conversion{$anchor_ups}="NA" if !defined $anchor_conversion{$anchor_ups};
#			$anchor_conversion{$anchor_downs}="NA" if !defined $anchor_conversion{$anchor_downs};
#			print "SP1\t$t_co1\t$gene_strand{$sp1}{$t_g1}\t$match1\t$anchor_ups\t$anchor_downs\t$cons_conv_call\t$array_sp1\n";
#			print "SP2\t$t_co2\t$gene_strand{$sp2}{$t_g2}\t$match2\t$anchor_conversion{$anchor_ups}\t$anchor_conversion{$anchor_downs}\t$cons_conv_call\t$array_sp2\n";
		    }
		    print OUT_GENES "$gene_cluster\t$t_g1\t$t_co1\t$le_ex1\t$t_exon_cl1\t$t_g2\t$t_co2\t$le_ex2\t$t_exon_cl2\t$cons_conv_call\t$sp1\t$sp2\n";
		    $cons_conv_call=~s/\_[A-Z]$//;
		    $tally_sp1_exons_in_Rcons_genes_by_type{$cons_conv_call}++;
		}
	    }
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

### Prints the conserved exons if option is activated
if (defined $outFile){
    open (OUT, ">Conserved_exons-$sp1-$sp2.tab");
    print OUT "Exon_cluster\tSpecies\tGeneID\tEx_Coord\tLength\tINFO\n";
    
    foreach my $exon_cluster (sort keys %exon_cluster_is_regulated){ # do also by gene (to check convergence)
	if ((defined $exon_cluster_has_sp{$exon_cluster}{$sp1}) && (defined $exon_cluster_has_sp{$exon_cluster}{$sp2})){
	    # here, there may be exons in clusters from sp1 and sp2 non-reg
	    foreach my $exon_sp1 (@{$exons_by_cluster{$exon_cluster}{$sp1}}){
		$info_by_exon{$exon_sp1}="NO_INFO" if (!defined $info_by_exon{$exon_sp1});
		my ($gene_sp1,$coord_sp1) = $exon_sp1 =~ /(.+?)\=(.+)/;
		my ($start,$end) = $coord_sp1 =~ /\:(\d+)\-(\d+)/;
		my $length = $end-$start+1;
		print OUT "$exon_cluster\t$sp1\t$gene_sp1\t$coord_sp1\t$length\t$info_by_exon{$exon_sp1}\n";
	    }
	    foreach my $exon_sp2 (@{$exons_by_cluster{$exon_cluster}{$sp2}}){
		$info_by_exon{$exon_sp2}="NO_INFO" if (!defined $info_by_exon{$exon_sp2});
		my ($gene_sp2,$coord_sp2) = $exon_sp2 =~ /(.+?)\=(.+)/;
		my ($start,$end) = $coord_sp2 =~ /\:(\d+)\-(\d+)/;
		my $length = $end-$start+1;
		print OUT "$exon_cluster\t$sp2\t$gene_sp2\t$coord_sp2\t$length\t$info_by_exon{$exon_sp2}\n";
	    }
	}
    }
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
    # perc of exons whose gene orthologs have regulated exons by type
    $tally_sp1_exons_in_Rcons_genes_by_type{CONSERVED}=0 if !defined $tally_sp1_exons_in_Rcons_genes_by_type{CONSERVED};
    $tally_sp1_exons_in_Rcons_genes_by_type{NON_CONSERVED}=0 if !defined $tally_sp1_exons_in_Rcons_genes_by_type{NON_CONSERVED};
    $tally_sp1_exons_in_Rcons_genes_by_type{UNCLEAR}=0 if !defined $tally_sp1_exons_in_Rcons_genes_by_type{UNCLEAR};
    my $total_exons_in_Rcons_genes = $tally_sp1_exons_in_Rcons_genes_by_type{CONSERVED}+$tally_sp1_exons_in_Rcons_genes_by_type{NON_CONSERVED}+$tally_sp1_exons_in_Rcons_genes_by_type{UNCLEAR};
    my $perc_sp1_exons_Rcons_genes_cons = sprintf ("%.2f", 100*$tally_sp1_exons_in_Rcons_genes_by_type{CONSERVED}/$total_exons_in_Rcons_genes);
    my $perc_sp1_exons_Rcons_genes_not = sprintf ("%.2f", 100*$tally_sp1_exons_in_Rcons_genes_by_type{NON_CONSERVED}/$total_exons_in_Rcons_genes);
    my $perc_sp1_exons_Rcons_genes_unclear = sprintf ("%.2f", 100*$tally_sp1_exons_in_Rcons_genes_by_type{UNCLEAR}/$total_exons_in_Rcons_genes);
    
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
    Total pairwise exon comparisons between regulated exons in $sp1 and $sp2 orthologs\t$total_exons_in_Rcons_genes
         Conserved exons in $sp2\t$tally_sp1_exons_in_Rcons_genes_by_type{CONSERVED}\t$perc_sp1_exons_Rcons_genes_cons\%
         Not conserved exons in $sp2\t$tally_sp1_exons_in_Rcons_genes_by_type{NON_CONSERVED}\t$perc_sp1_exons_Rcons_genes_not\%
         Unclear in $sp2\t$tally_sp1_exons_in_Rcons_genes_by_type{UNCLEAR}\t$perc_sp1_exons_Rcons_genes_unclear\%
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
