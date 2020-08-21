#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;

### AIMS OF SCRIPT
#* to split the gene clusters into smaller bits (e.g. 500 clusters per file) and creates a folder GENE_CLUSTERS in the species pair folder
#* prepare folders for a pair of species
#* submit jobs to score exons for each sub-file of clusters for the given species pair
#* submit jobs to score introns for each sub-file of clusters for the given species pair
#* score introns and exons [these are the two main scripts]: Score_exons_pair_sp.pl and Score_introns_pair_sp.pl
         

#### DECLARATION AND GENERATION OF VARIABLES
my $clean; # forces full deletion of pre-existing runs
my $f_gene_cluster; # file with the gene clusters for the give species pair => it may need additional and specific rules...
my $f_extra_cluster; # a different cluster file for some species (hard to generalize); eventually allow it to be comma separated
my $sp1; ##species1 
my $sp2; ##species2
my $EX_DB = "./EXONS_DB"; # default where the EXONS_DB/ folder is
my $project_dir  = "./"; ##output folder
my $verbose; # prints steps
my $N=10000;
my $w=0;
my $help;

### Get options:
Getopt::Long::Configure("no_auto_abbrev");
GetOptions( "gene_cluster=s" => \$f_gene_cluster, ##gene clusters for sp1 and sp2
	    "sp1=s" => \$sp1,
	    "sp2=s" => \$sp2,
	    "clean" => \$clean,
	    "N_split=i" => \$N,
	    "expath=s" => \$EX_DB,
	    "project_dir=s" => \$project_dir,
	    "verbose" => \$verbose,
	    "help" => \$help
    );

### Help
if (!defined ($f_gene_cluster) || !defined($project_dir) || !defined($sp1) || !defined($sp2)|| defined ($help)){
    die "
Usage: Prepare_and_Submit_Aln_sp_pair.pl -s1 sp1 -s2 sp2 --gene_cluster FILE --expath Exons_DB_folder  --project_dir Output_folder [options]

OPTIONS
    --clean            Force removal of previous data [def OFF]
    --N_split int      Number of alignments in subfile [def 10000]
    --sp1	       Species 1
    --sp2	       Species 2
    --expath PATH      EXONS_DB folder with the output obtained in Module I for each species (def = EXONS_DB/)
    --project_dir PATH Output folder: the species pair folder will be created here (def = ./)
    --verbose          Prints commands in screen [def OFF]
    
";
}

my @species=($sp1,$sp2);
@species=sort(@species); # sort alphanumerically
$sp1=$species[0];
$sp2=$species[1];
### Stores the relevant files from EXONS_DB/Sp
my %f_exint;
my %f_protIDs;
my %f_exon_pos;
my %f_intron_pos;

##Files Sp1
$f_exint{$sp1}="$EX_DB/$sp1/$sp1".".exint";
$f_protIDs{$sp1}="$EX_DB/$sp1/$sp1"."_annot_exons_prot_ids.txt"; # before _annot_fake_protein_ids_exons.txt"
$f_exon_pos{$sp1}="$EX_DB/$sp1/$sp1"."_protein_ids_exons_pos.txt";
$f_intron_pos{$sp1}="$EX_DB/$sp1/$sp1"."_protein_ids_intron_pos_CDS.txt";

die "Cannot find exint file for $sp1\n" unless (-e $f_exint{$sp1});
die "Cannot find prot IDs file for $sp1\n" unless (-e $f_protIDs{$sp1});
die "Cannot find exon pos file for $sp1\n" unless (-e $f_exon_pos{$sp1});
die "Cannot find intron pos file for $sp1\n" unless (-e $f_intron_pos{$sp1});

##Files Sp2

$f_exint{$sp2}="$EX_DB/$sp2/$sp2".".exint";
$f_protIDs{$sp2}="$EX_DB/$sp2/$sp2"."_annot_exons_prot_ids.txt"; # before _annot_fake_protein_ids_exons.txt"
$f_exon_pos{$sp2}="$EX_DB/$sp2/$sp2"."_protein_ids_exons_pos.txt";
$f_intron_pos{$sp2}="$EX_DB/$sp2/$sp2"."_protein_ids_intron_pos_CDS.txt";

die "Cannot find exint file for $sp2\n" unless (-e $f_exint{$sp2});
die "Cannot find prot IDs file for $sp2\n" unless (-e $f_protIDs{$sp2});
die "Cannot find exon pos file for $sp2\n" unless (-e $f_exon_pos{$sp2});
die "Cannot find intron pos file for $sp2\n" unless (-e $f_intron_pos{$sp2});

### Setting up and creating PROJECT directory

my $pair_folder=$sp1."_".$sp2;

if (defined ($verbose)){
    print "\nrm -r $project_dir/$pair_folder\n" if $clean && (-e "$project_dir/$pair_folder"); #Forces full deletion of data
    print "\nmkdir $project_dir/$pair_folder\n" unless (-e "$project_dir/$pair_folder");
}
system "rm -rf $project_dir/$pair_folder" if $clean && (-e "$project_dir/$pair_folder"); #Forces full deletion of data
system "mkdir $project_dir/$pair_folder" unless (-e "$project_dir/$pair_folder");

my %N_parts; # stores the number of splits for each file of clusters
my $gcl=0;

### Splitting clusters

&split_cluster($f_gene_cluster,$f_exint{$sp1},$f_exint{$sp2},$N,$sp1,$sp2,$project_dir,$pair_folder) if ($gcl==0); 
print "\n$f_gene_cluster parts: $N_parts{$f_gene_cluster}\n" if defined ($verbose);

#### SUBROUTINES

sub split_cluster {
    my @input = @_;
    my $full_cluster = $input[0];
    my $ex1=$input[1]; # exint sp1
    my $ex2=$input[2]; # exint sp2
    my $N=$input[3];
    my $sp1=$input[4];
    my $sp2=$input[5];
    my $project_dir=$input[6];
    my $pair_folder=$input[7];
    my $cluster_root = $full_cluster;
    my (%trs, %nt, %warn, %nalns);
    $cluster_root =~ s/.+\///; # removes the path
    $cluster_root = "$project_dir/$pair_folder/$cluster_root";
    my %size;
    my $part=1;
    my %cid;
    my $n=0;
    open (ONE, "$ex1") || die "It cannot open the exint file for $sp1 ($ex1)\n";
    while (<ONE>){
    	chomp($_);
    	if ($_=~/\>/){
	    my @l=split(/\s+/,$_);
	    my @t=split(/\|/,$l[0]);
	    $trs{$t[1]}++; # the exint file is the one that decides what to align/how to split
    	}
    }
    close ONE;
    open (TWO, "$ex2") || die "It cannot open the exint file for $sp2 ($ex2)\n";
    while (<TWO>){
	chomp($_);
    	if ($_=~/\>/){
	    my @l=split(/\s+/,$_);
	    my @t=split(/\|/,$l[0]);
	    $trs{$t[1]}++;
    	}
    }
    close TWO;

    #GF000001 Mmu ENSMUSG00000090353	Gm17555O	G0000001	O=0,R=0
    open (CL, "$full_cluster") || die "It cannot open full clusters\n";
    while (<CL>){
    	chomp($_);
    	my @l=split(/\t/,$_);
    	$nt{$l[0]}{$l[1]}+=$trs{$l[2]};
    	$cid{$l[0]}=1;
    }
    close CL;
    my @k=sort(keys(%cid));
    my $el;
    my $taln=0;
    my $max;
    open (NALN, ">$project_dir/$pair_folder/Num_alns_by_cl_$sp1-$sp2.log");
    foreach $el (@k){
    	my $naln=$nt{$el}{$sp1}*$nt{$el}{$sp2};
    	$taln+=$naln;
    	if ($naln>$N){  
	    print NALN "$el\t$naln\tThis cluster exceed $N alignments: Reduce number of total transcripts in exint file for this cluster or increase number of alignments!!!\n";
	    $warn{$el}=1;
	    $w++;
    	}
    	else { print NALN "$el\t$naln\n"; }
    	$nalns{$el}=$naln;
    	if (!$max){ $max=$naln; }
    	elsif ($naln>$max){ $max=$naln; }
    }
    print NALN "\nTotal alignments\t$taln\t$max\nNumber of clusters with warnings\t$w\n";
    if ($w){ 
	print STDERR "\n\nWARNING!!! $w clusters exceed the maximum number of alignments: $N  !!!!\n\n"; 
	#print STDERR "WARNING!!! Skipping $w clusters in the splitted files !!!\n\n"; 
    }	
    my $sum=0;
    my %prcid;
    $part=1;
    open (CL, "$full_cluster") || die "It cannot open clusters file\n";
    while (<CL>){
	chomp($_);
	my @l=split(/\t/,$_);
	if (!$warn{$l[0]}){
	    open (OUT, ">>$cluster_root-part_"."$part"); # first part already defined
	    if ($sum<$N){
		if (!exists($prcid{$l[0]})){
		    $sum+=$nalns{$l[0]};
		    $prcid{$l[0]}=1;
		}
		print OUT "$_\n";
	    }				
	    else {       
		if ($prcid{$l[0]}){
		    print OUT "$_\n";
		}
		elsif (!exists($prcid{$l[0]})){
		    $sum=$nalns{$l[0]};
		    $part++;
		    open (OUT, ">>$cluster_root-part_"."$part");
		    print OUT "$_\n";
		    $prcid{$l[0]}=1;
		}
	    }
	}
    }
    close CL;

    ##printing files of clusters with > N alignments
    my %print;
    open (CL, "$full_cluster") || die "It cannot open clusters file\n";
    while (<CL>){
	chomp($_);
	my @l=split(/\t/,$_);
	if ($warn{$l[0]}){
	    if (!$print{$l[0]}){
		$part++;
		$print{$l[0]}=1
	    }
	    open (OUT, ">>$cluster_root-part_"."$part"); # first part already defined
	    print OUT "$_\n";
	}
    }		
    close CL;
    $N_parts{$full_cluster}=$part;
    close OUT;
}
