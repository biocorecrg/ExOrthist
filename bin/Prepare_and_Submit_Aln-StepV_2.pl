#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;

### AIMS OF SCRIPT
#* to split the gene clusters into smaller bits (e.g. 1000 clusters per file)
#* prepare folders for each pair of species
#* makes several symbolic links that will be needed in downstream scripts => should be able to get them from EXONS_DB/
#* submit jobs to score exons for each sub-file of clusters and pair of species
#* submit jobs to score introns for each sub-file of clusters and pair of species
#* score introns and exons [these are the two main scripts]: Score_exons_pair_sp.pl and Score_introns_pair_sp.pl

### FOLDER STRUCTURE
#* Root:                           $dir/
#* Common to different projects:   $dir/EXONS_DB/ with all the general info per species
#* Project specific:               $dir/PROJECT_Sp1_Sp2_Sp3/ => made alphabetically [will avoid bugs]
#   * by pairs:                    $dir/PROJECT_Sp1_Sp2_Sp3/Sp1_Sp2/ etc
#   * for gene cluster parts:      $dir/PROJECT_Sp1_Sp2_Sp3/GENE_CLUSTERS/          


#### DECLARATION AND GENERATION OF VARIABLES
my $clean; # forces full deletion of pre-existing runs
my $f_gene_cluster; # file with the gene clusters => it may need additional and specific rules...
my $f_extra_cluster; # a different cluster file for some species (hard to generalize); eventually allow it to be comma separated
my $species_array; ##species separated by comma 
my $submit_aln; # whether to submit the jobs
my $bin; # default (folder where files and program of ORTH_PIPE are located)
my $dir;
my $verbose; # prints steps
my $N_split = 500; # number of clusters per split
my $help;

### Get options:
Getopt::Long::Configure("no_auto_abbrev");
GetOptions(  "gene_cluster=s" => \$f_gene_cluster,
	     "extra_cluster=s" => \$f_extra_cluster,
	     "sp=s" => \$species_array,
	     "clean" => \$clean,
	     "submit_aln" => \$submit_aln,
	     "dir=s" => \$dir,
	     "bin=s" => \$bin,
	     "verbose" => \$verbose,
	     "N_split=i" => \$N_split,
	     "help" => \$help
    );

### Help
if (!defined ($f_gene_cluster) || !defined($bin) || !defined($dir) || !defined($species_array) || defined ($help)){
    die "
Usage: Prepare_and_Submit_Aln-StepV_2.pl -sp Sp1,Sp2,Sp3 --gene_cluster FILE [options]

OPTIONS
    --clean                   Force removal of previous data [def OFF]
    --submit_aln              Submit the jobs of for scoring ex/introns [def OFF]
    --dir PATH                EXONS_DB folder with the output obtained in Module I for each species.
    --bin PATH                Bin directory where the exon orthology pipeline is located                 
    --verbose                 Prints commands in screen [def OFF]
    --N_split int             Number of cluster per subfile in jobs [def 1000]       
    --extra_cluster  FILE     Additional file for certain species [def OFF]
    

*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

my @species=split(/\,/,$species_array);
@species=sort(@species); # sort alphanumerically

### Stores the relevant files from EXONS_DB/Sp
my %f_exint;
my %f_protIDs;
my %f_exon_pos;
my %f_intron_pos;

foreach my $sp (@species){ # technically speaking, the downstream scripts should handle with just providing Sp1 and Sp2
    # I removed the hardcoded EXONS_DB
    $f_exint{$sp}="$dir/$sp/$sp".".exint";
    $f_protIDs{$sp}="$dir/$sp/$sp"."_annot_exons_prot_ids.txt"; # before _annot_fake_protein_ids_exons.txt"
    $f_exon_pos{$sp}="$dir/$sp/$sp"."_protein_ids_exons_pos.txt";
    $f_intron_pos{$sp}="$dir/$sp/$sp"."_protein_ids_intron_pos_CDS.txt";

    die "Cannot find exint file for $sp\n" unless (-e $f_exint{$sp});
    die "Cannot find prot IDs file for $sp\n" unless (-e $f_protIDs{$sp});
    die "Cannot find exon pos file for $sp\n" unless (-e $f_exon_pos{$sp});
    die "Cannot find intron pos file for $sp\n" unless (-e $f_intron_pos{$sp});
}

### Setting up and creating PROJECT directory
my $temp_sp = join("_",@species);
my $project_dir = "PROJECT_".$temp_sp;
if (defined $verbose){
    print "\nrm -r $dir/$project_dir\n" if (defined ($clean) && -e "$dir/$project_dir");
    print "mkdir $dir/$project_dir\n" unless (-e "$dir/$project_dir");
    print "mkdir $dir/$project_dir/GENE_CLUSTERS\n" unless (-e "$dir/$project_dir/GENE_CLUSTERS");
}
my %N_parts; # stores the number of splits for each file of clusters
my $gcl=0;
system "rm -r $dir/$project_dir" if (defined ($clean) && -e "$dir/$project_dir");
system "mkdir $dir/$project_dir" unless (-e "$dir/$project_dir");
if (-e "$dir/$project_dir/GENE_CLUSTERS"){ 
	$gcl=1; 
	my $parts=`ls $dir/$project_dir/GENE_CLUSTERS/*part* | wc -l`;
	$N_parts{$f_gene_cluster}=$parts;

} else { system "mkdir $dir/$project_dir/GENE_CLUSTERS"; }

### Splitting clusters

&split_cluster($f_gene_cluster) if ($gcl==0); 
&split_cluster($f_extra_cluster) if (defined $f_extra_cluster);

print "\n$f_gene_cluster parts: $N_parts{$f_gene_cluster}\n" if defined ($verbose);

### Prepares folders and submissions
for my $i (0..$#species){
    my $sp1 = $species[$i];

    for my $j (0..$#species){
	my $sp2 = $species[$j];

	if ($i < $j){ # only does Sp1_Sp2 (not Sp2_Sp1)
	    my $pair_folder=$sp1."_".$sp2;
	    
	    if (defined ($verbose)){
		print "\nrm -r $dir/$project_dir/$pair_folder\n" if $clean && (-e "$dir/$project_dir/$pair_folder"); #Forces full deletion of data
		print "\nmkdir $dir/$project_dir/$pair_folder\n" unless (-e "$dir/$project_dir/$pair_folder");
	    }
	    system "rm -rf $dir/$project_dir/$pair_folder" if $clean && (-e "$dir/$project_dir/$pair_folder"); #Forces full deletion of data
	    system "mkdir $dir/$project_dir/$pair_folder" unless (-e "$dir/$project_dir/$pair_folder");
	    
	    ### here it could read the file with ALL previous aln to weight the hours in job submission
	    
	    my $temp_cl_file = $f_gene_cluster; # used by default
	    
	    ### hard-coded behaviour for some species:
	    if (defined ($f_extra_cluster)){
		if (($sp1 eq "Bla") || ($sp2 eq "Bla")) {$temp_cl_file = $f_extra_cluster;}
		elsif (($sp1 eq "Spu") || ($sp2 eq "Spu")){$temp_cl_file = $f_extra_cluster;}
		elsif (($sp1 eq "Sma") || ($sp2 eq "Sma")){$temp_cl_file = $f_extra_cluster;}
		elsif (($sp1 eq "Obi") || ($sp2 eq "Obi")){$temp_cl_file = $f_extra_cluster;}
		elsif (($sp1 eq "Dme") || ($sp2 eq "Dme")){$temp_cl_file = $f_extra_cluster;}
	    }
	    
	    for my $part (1..$N_parts{$temp_cl_file}){
		$temp_cl_file=~s/.+\///;
		$temp_cl_file="$dir/$project_dir/GENE_CLUSTERS/$temp_cl_file";
		my $cl_part = $temp_cl_file."-part_".$part;
		my $tb=$bin;
		$tb=~s/bin/files/;		
		my $bl62 = "$tb/blosum62.txt";
		my $of=$dir."/".$project_dir."/".$pair_folder;
		### Submits for exons and introns
		if (defined ($verbose)){
		    print "\nsubmitjob long -l h_rt=50:00:00,virtual_free=30G ".
			"perl $bin/Score_exons_introns_pair_sp.pl $sp1 $sp2 $cl_part ".
			"$f_protIDs{$sp1} $f_protIDs{$sp2} $f_exon_pos{$sp1} $f_exon_pos{$sp2} $f_intron_pos{$sp1} $f_intron_pos{$sp2} $f_exint{$sp1} $f_exint{$sp2} $part $bin $bl62 $of\n";
		}
		
		system "submitjob long -l h_rt=50:00:00,virtual_free=30G ". # job conditions
		    "perl $bin/Score_exons_introns_pair_sp.pl $sp1 $sp2 $cl_part ". # species and gene_cluster
		    "$f_protIDs{$sp1} $f_protIDs{$sp2} $f_exon_pos{$sp1} $f_exon_pos{$sp2} $f_intron_pos{$sp1} $f_intron_pos{$sp2} $f_exint{$sp1} $f_exint{$sp2} $part $bin $bl62 $of\n" if $submit_aln; # input files and outputs
		#}
	    }
	}	
    }
}
    


#### SUBROUTINES

sub split_cluster {
    my @input = @_;
    my $full_cluster = $input[0];
    my $cluster_root = $full_cluster;
    $cluster_root =~ s/.+\///; # removes the path
    $cluster_root = "$dir/$project_dir/GENE_CLUSTERS/$cluster_root";
    my %size;
    my $part=1;
    my %cid;
    my $n=0;
    open (CL, $full_cluster) || die "Cannot open $full_cluster file\n";
    while (<CL>){ 
		chomp($_);
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

			$n++;
			$cid{$line[0]}=1;
			print OUT "$_\n";
			$n=0;
			$part++;

		}
		close (OUT);	

		}
    close (CL);
    $N_parts{$full_cluster}=$part;

}
