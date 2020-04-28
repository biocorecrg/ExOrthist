#!/usr/bin/perl
#use warnings;
#use strict;
use Getopt::Long;
use Data::Dumper;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions( "gene_cluster=s" => \$f_gene_cluster,
	         "exon_cluster=s" => \$f_exon_cluster,
	         "exon_list=s" => \$f_exon_list,
	         "s1=s" => \$sp1,
	         "s2=s" => \$sp2,
	         "outfile=s" => \$out,
	         "help" => \$help
    );

### Help
if (!defined ($f_gene_cluster) || !defined($f_exon_cluster) || !defined($f_exon_list) || !defined($out) || !defined($sp1) || !defined($sp2) ||  defined ($help)){
    die "
Usage: Get_cons_exons2species.pl -s1 sp1 [species_query] -s2 sp2 [species_subject] --gene_cluster FILE --exon_cluster FILE --exon_list FILE [exon_list_species_query] --outfile FILE
";
}

#GF0000002Mm2ENSMUSG00000069307Hist1h2bq
open (GCL, $f_gene_cluster);
while (<GCL>){
    chomp($_);
    @l=split(/\t/,$_);
    if ($l[1] eq $sp1 || $l[1] eq $sp2) {
	$gcl{$l[1]}{$l[2]}=$l[0];
	$gcid{$l[0]}{$l[1]}=1;
	$cid{$l[2]}=$l[0];
    }
}
#GF0000007.001ENSMUSG00000096688chr4:61595665-61595760:-Mm2
open (ECL, $f_exon_cluster);
while (<ECL>){
    chomp($_);
    @l=split(/\t/,$_);
    if($l[3] eq $sp1 || $l[3] eq $sp2) {
	$l[2]=~s/\-/\:/;
	@tmp=split(/\:/,$l[2]);
	$eid=$l[1].":".$tmp[1];
	$ecl{$eid}=$l[0];
	$eid=$l[1].":".$tmp[2];
	$ecl{$eid}=$l[0];
	$esp{$l[3]}{$l[0]}=1;
    }
}
#ENSMUSG00000030279chr6:143072543-143072694DOWN
open (LIST, $f_exon_list);
while (<LIST>){
    chomp($_);
    @l=split(/\t/,$_);
    $l[1]=~s/\-/\:/;
    $exsp1++;
    $rcl="";
    $gsp1{$l[0]}=1;
    @tmp=split(/\:/,$l[1]);
    $eid1=$l[0].":".$tmp[1];
    $eid2=$l[0].":".$tmp[2];
    ##1 Check first if the exon falls in an orthologous gene 
    if ($cid{$l[0]}){
	$gc=$cid{$l[0]};
	if ($gcid{$gc}{$sp1}){
	    $cons_gcl{$gc}=1; ##saving those gene clusters in common
	    $gcons{$l[0]}=$gc; ##conserved genes
	    $gex++;
	}
    }
    if ($ecl{$eid1}){
	$rcl=$ecl{$eid1};
	if ($esp{$sp2}{$rcl}){
	    $cons_ecl{$rcl}=1; ##saving exon clusters in common
	    $excons{$l[0]."\t".$l[1]}=$rcl;
	    $cex++;
	}
    }
    elsif ($ecl{$eid2}){
	if ($esp{$sp2}{$rcl}){
	    $cons_ecl{$rcl}=1; ##saving exon clusters in common
	    $excons{$l[0]."\t".$l[1]}=$rcl;
	    $cex++;
	}
    }
}

##printing 
@k0=keys(%gsp1); ##number of conserved genes
@k1=keys(%cons_gcl); ##number of gene clusters
@k2=keys(%gcons); ##number of conserved genes
@k3=keys(%cons_ecl); ##number of conserved genes
$tg=scalar(@k2);
$tgc=scalar(@k1);
$tgs1=scalar(@k0);
$tex=scalar(@k3);
open (OUT, ">$out"); ##outfile with stats
print OUT "\nTotal genes of species $sp1 in exon list\t$tgs1\n";
print OUT "Total exons of species $sp1 in exon list\t$exsp1\n\n";
print OUT "### STATS OF GENES ###\n\n";
print OUT "Total genes hosting exons in $sp1 having a gene orthologue in $sp2\t$tg\n";
print OUT "Total gene clusters with exons of $sp1 with $sp2 genes\t$tgc\n\n";
print OUT "### STATS OF EXONS ###\n\n";
print OUT "Total exons in $sp1 conserved with exons in $sp2\t$cex\n";
print OUT "Total exon clusters with exons of $sp1 with $sp2 exons\t$tex\n\n";

open (GENES, ">Conserved_genes_$sp1-$sp2.tab");
open (GCL, $f_gene_cluster);
while (<GCL>){
    chomp($_);
    @l=split(/\t/,$_);
    if ($l[1] eq $sp1 || $l[1] eq $sp2) {
	if ($cons_gcl{$l[0]}){
	    print GENES "$_\n";
	}
    }
}

open (EXONS, ">Conserved_exons_$sp1-$sp2.tab");
#GF0000007.001ENSMUSG00000096688chr4:61595665-61595760:-Mm2
open (ECL, $f_exon_cluster);
while (<ECL>){
    chomp($_);
    @l=split(/\t/,$_);
    if($l[3] eq $sp1 || $l[3] eq $sp2) {
	if ($cons_ecl{$l[0]}){
	    print EXONS "$_\n";
	}
    }
}






