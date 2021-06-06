#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $main_folder;
my ($gene1,$gene2,$exon1,$exon2);
my ($sp1,$sp2);
my ($prot1,$prot2);
my ($ex_number1,$ex_number2);
my $outFile;
my $helpFlag;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions( "main_folder=s" => \$main_folder,
	    "g1=s" => \$gene1,
	    "gene1=s" => \$gene1,
	    "g2=s" => \$gene2,
	    "gene2=s" => \$gene2,
	    "ex1=s" => \$exon1, 
	    "exon1=s" => \$exon1, 
	    "ex2=s" => \$exon2, 
	    "exon2=s" => \$exon2, 
	    "sp1=s" => \$sp1,
	    "sp2=s" => \$sp2,
	    "outFile" => \$outFile,
	    "helpFlag" => \$helpFlag
    );


### Help
if (!defined ($main_folder) || !defined ($gene1) || !defined ($gene2) || !defined ($exon1) || !defined ($exon2) || !defined($sp1) || !defined($sp2) ||  defined ($helpFlag)){
    die "
Usage: retrieve_IPA.pl -g1 gene_sp1 -g2 gene_sp2 -ex1 exon_sp1 -ex2 exon_sp2 -sp1 species1 -sp2 species2 -main_folder main_folder/

Script to find and retrive the most representative Intron Position Aware alignment for a pair of exons

Options:
     -sp1 species 1             Identifier of the first species.
     -sp2 species 2             Identifier of the second species.
     -gene1/g1 geneID           Gene ID for Sp1
     -gene2/g2 geneID           Gene ID for Sp2
     -exon1/ex1 exon_coord      Exon coordinate for Sp1 (chr:start-end)
     -exon2/ex2 exon_coord      Exon coordinate for Sp2 (chr:start-end)  
     -main_folder folder        Main folder with the ouput of ExOrthist 
     -outFile                   Creates an output file with IPA (IPA-g1:ex1-g2:ex2-sp1-sp2.aln)
     -helpFlag                  Prints this help message
  

";
}

# removes strand in exon if provided
$exon1=~s/\:[\+\-]//;
$exon2=~s/\:[\+\-]//;

# sets species order
my $t_spA = (sort($sp1,$sp2))[0];
my $t_spB = (sort($sp1,$sp2))[1];

# input files
my $best_hits="$main_folder/$t_spA-$t_spB/best_scored_EX_matches_by_targetgene.txt";
my $ALN="$main_folder/$t_spA-$t_spB/EXINT_aln.gz";
open (BEST, $best_hits) || die "It cannot open $best_hits\n";
open (ALN, "gunzip -c $ALN |") || die "It cannot find $ALN\n";

# output file
if (defined $outFile){
    $outFile="IPA-$gene1:$exon1-$gene2:$exon2-$sp1-$sp2.aln";
    open (O, ">$outFile") || die "It cannot open the output file ($outFile)\n";
}

<BEST>;
LBL1:while (<BEST>){
    chomp;
    my @t=split(/\t/);
    my ($temp_g1)=$t[3]=~/.+?\|(.+)/;
    my $t_ex1=$t[5];
    $t_ex1=~s/\:[\+\-]//;
    my ($temp_g2)=$t[6]=~/.+?\|(.+)/;
    my $t_ex2=$t[8];
    $t_ex2=~s/\:[\+\-]//;
    
    if ($temp_g1 eq $gene1 && $temp_g2 eq $gene2 && $t_ex1 eq $exon1 && $t_ex2 eq $exon2){
	$prot1=$t[3];
	$prot2=$t[6];
	$ex_number1=$t[4];
	$ex_number2=$t[7];
    }
    last LBL1 if $prot1 && $prot2;
}
close BEST;

die "It could not find the match for $gene1:$exon1 and $gene2:$exon2 (try reverting species query/target)\n" if !$prot1 || !$prot2;
#print "Protein pair identified: $prot1 and $prot2\n" if $prot1 && $prot2;

# inverts protein order if needed
my ($t_protA, $t_protB) = ($prot1,$prot2);
$t_protA=$prot2 if $t_spB eq $sp1;
$t_protB=$prot1 if $t_spA eq $sp2;

# >>> dm6 mm10 FBpp0072009|FBgn0011296 ENSMUSP00000005077|ENSMUSG00000004951
my $aln_data;
my $activate;
LBL2:while (<ALN>){
    if (/\>\>\>/ && !$aln_data){
	if (/\>\>\> $t_spA $t_spB $t_protA $t_protB/){
	    $aln_data=">>> $sp1: $prot1 ($ex_number1) <=> $sp2: $prot2 ($ex_number2)\n";
	}
    }
    elsif (/\>\>\>/ && $aln_data){
	last LBL2;
    }
    elsif ($aln_data){
	if ($activate){
	    $aln_data.=$_;
	}
	elsif (/> /){
	    $activate=1;
	}
    }
}
close ALN;

if ($aln_data){
    print "\n$aln_data\n";
    
    if (defined $outFile){
	print O "$aln_data\n";
	close O;
    }
} else {
    print "It could not identify an alignment for $t_protA and $t_protB\n";
}
