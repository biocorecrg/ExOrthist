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
my $f_gene_cluster; # file with the gene clusters for the give species pair => it may need additional and specific rules...
my $f_extra_cluster; # a different cluster file for some species (hard to generalize); eventually allow it to be comma separated
my $sp1; ##species1 
my $sp2; ##species2
my $submit_aln; # whether to submit the jobs
my $expath;
my $project_dir; ##output folder
my $verbose; # prints steps
my $N_split = 500; # number of clusters per split
my $help;

### Get options:
Getopt::Long::Configure("no_auto_abbrev");
GetOptions( "gene_cluster=s" => \$f_gene_cluster, ##gene clusters for sp1 and sp2
	     "sp1=s" => \$sp1,
	     "sp2=s" => \$sp2,
	     "submit_aln" => \$submit_aln,
	     "expath=s" => \$expath,
	     "project_dir=s" => \$project_dir,
	     "verbose" => \$verbose,
	     "N_split=i" => \$N_split,
	     "help" => \$help
    );

### Help
if (!defined ($f_gene_cluster) || !defined($expath) || !defined($project_dir) || !defined($sp1) || !defined($sp2)|| defined ($help)){
    die "
Usage: Prepare_and_Submit_Aln_sp_pair.pl -sp1 sp1 -sp2 sp2 --gene_cluster FILE --expath Exons_DB_folder  --project_dir Output_folder 

OPTIONS
    --submit_aln              Submit the jobs of for scoring ex/introns [def OFF]
    --sp1		      Species 1
    --sp2		      Species 2
    --expath PATH             EXONS_DB folder with the output obtained in Module I for each species.
    --project_dir PATH        Output folder: the species pair folder will be created here
    --verbose                 Prints commands in screen [def OFF]
    --N_split int             Number of cluster per subfile in jobs [def 500]       
    
*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

my @species=($sp1,$sp2);
@species=sort(@species); # sort alphanumerically
##
$sp1=$species[0];
$sp2=$species[1];
### Stores the relevant files from EXONS_DB/Sp
my %f_exint;
my %f_protIDs;
my %f_exon_pos;
my %f_intron_pos;

##Files Sp1

$f_exint{$sp1}="$expath/$sp1/$sp1".".exint";
$f_protIDs{$sp1}="$expath/$sp1/$sp1"."_annot_exons_prot_ids.txt"; # before _annot_fake_protein_ids_exons.txt"
$f_exon_pos{$sp1}="$expath/$sp1/$sp1"."_protein_ids_exons_pos.txt";
$f_intron_pos{$sp1}="$expath/$sp1/$sp1"."_protein_ids_intron_pos_CDS.txt";

die "Cannot find exint file for $sp1\n" unless (-e $f_exint{$sp1});
die "Cannot find prot IDs file for $sp1\n" unless (-e $f_protIDs{$sp1});
die "Cannot find exon pos file for $sp1\n" unless (-e $f_exon_pos{$sp1});
die "Cannot find intron pos file for $sp1\n" unless (-e $f_intron_pos{$sp1});

##Files Sp2

$f_exint{$sp2}="$expath/$sp2/$sp2".".exint";
$f_protIDs{$sp2}="$expath/$sp2/$sp2"."_annot_exons_prot_ids.txt"; # before _annot_fake_protein_ids_exons.txt"
$f_exon_pos{$sp2}="$expath/$sp2/$sp2"."_protein_ids_exons_pos.txt";
$f_intron_pos{$sp2}="$expath/$sp2/$sp2"."_protein_ids_intron_pos_CDS.txt";

die "Cannot find exint file for $sp2\n" unless (-e $f_exint{$sp2});
die "Cannot find prot IDs file for $sp2\n" unless (-e $f_protIDs{$sp2});
die "Cannot find exon pos file for $sp2\n" unless (-e $f_exon_pos{$sp2});
die "Cannot find intron pos file for $sp2\n" unless (-e $f_intron_pos{$sp2});

### Setting up and creating PROJECT directory

my $pair_folder=$sp1."_".$sp2;


my %N_parts; # stores the number of splits for each file of clusters
my $gcl=0;
system "mkdir $project_dir/$pair_folder/"; 

### Splitting clusters

&split_cluster($f_gene_cluster) if ($gcl==0); 

print "\n$f_gene_cluster parts: $N_parts{$f_gene_cluster}\n" if defined ($verbose);

### Prepares submissions
### here it could read the file with ALL previous aln to weight the hours in job submission
	    
my $temp_cl_file = $f_gene_cluster; # used by default
for my $part (1..$N_parts{$temp_cl_file}){
		$temp_cl_file=~s/.+\///;
		$temp_cl_file="$project_dir/$pair_folder/$temp_cl_file";
		my $cl_part = $temp_cl_file."-part_".$part;
		my $tb="";
		$tb=~s/bin/files/;		
		my $bl62 = "$tb/blosum62.txt";
		my $of=$project_dir."/".$pair_folder;
		### Submits for exons and introns
		if (defined ($verbose)){
		    print "\nsubmitjob long -l h_rt=50:00:00,virtual_free=30G ".
			"perl /Score_exons_introns_pair_sp.pl $sp1 $sp2 $cl_part ".
			"$f_protIDs{$sp1} $f_protIDs{$sp2} $f_exon_pos{$sp1} $f_exon_pos{$sp2} $f_intron_pos{$sp1} $f_intron_pos{$sp2} $f_exint{$sp1} $f_exint{$sp2} $part  $bl62 $of\n";
		}
		
		system "submitjob long -l h_rt=50:00:00,virtual_free=30G ". # job conditions
		    "perl /Score_exons_introns_pair_sp.pl $sp1 $sp2 $cl_part ". # species and gene_cluster
		    "$f_protIDs{$sp1} $f_protIDs{$sp2} $f_exon_pos{$sp1} $f_exon_pos{$sp2} $f_intron_pos{$sp1} $f_intron_pos{$sp2} $f_exint{$sp1} $f_exint{$sp2} $part  $bl62 $of\n" if $submit_aln; # input files and outputs

}
    
#### SUBROUTINES

sub split_cluster {
    my @input = @_;
    my $full_cluster = $input[0];
    my $cluster_root = $full_cluster;
    $cluster_root =~ s/.+\///; # removes the path
    $cluster_root = "$project_dir/$pair_folder/$cluster_root";
    my %size;
    my $part=1;
    my %cid;
    my $n=0;
    open (CL, $full_cluster) || die "Cannot open $full_cluster file\n";
    while (<CL>){ 
		chomp($_);
		if ($_ ne "") { 
			my @line = split(/\t/,$_);    	
				open (OUT, ">>$cluster_root-part_"."$part"); # first part already defined
			if ($n<$N_split){
				if (!exists($cid{$line[0]})){
					$n++;
					$cid{$line[0]}=1;
				}
				print OUT "$_\n";
			}
			else {
				#$n++;
				$cid{$line[0]}=1;
				print OUT "$_\n";
				$n=0;
				$part++;
	
			}
		}
		close (OUT);	
	}
    close (CL);
    $N_parts{$full_cluster}=$part;
}
