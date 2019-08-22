#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

### AIMS OF SCRIPT
#* to prepare and submit the jobs to obtain the final score of an exon between two species.

### FOLDER STRUCTURE
#* Root:                           $dir/
#* Common to different projects:   $dir/EXONS_DB/ with all the general info per species
#* Project specific:               $dir/PROJECT_Sp1_Sp2_Sp3/ => made alphabetically [will avoid bugs]
#   * by pairs:                    $dir/PROJECT_Sp1_Sp2_Sp3/Sp1_Sp2/ etc
#   * for gene cluster parts:      $dir/PROJECT_Sp1_Sp2_Sp3/GENE_CLUSTERS/ 

#### DECLARATION AND GENERATION OF VARIABLES
my $species_array; ##species separated by comma 
my $submit_jobs; # whether to submit the jobs
my $bin;
my $dir;
my $verbose; # prints steps
my $help;

### Get options:
Getopt::Long::Configure("no_auto_abbrev");
GetOptions(  "sp=s" => \$species_array,
	     "submit_jobs" => \$submit_jobs,
	     "dir=s" => \$dir,
	     "bin=s" => \$bin,
	     "verbose" => \$verbose,
	     "help" => \$help
    );


if (!defined ($dir) || !defined($species_array) || defined ($help)){
    die "
Usage: Run_score_by_exon.pl -sp Sp1,Sp2,Sp3 --dir PATH [options]

OPTIONS
    --clean                   Force removal of previous data [def OFF]
    --submit_jobs              Submit the jobs of for scoring ex/introns [def OFF]
    --dir PATH                Root directory where project will be created. 
    --bin PATH               Bin directory where scripts are located
    --verbose                 Prints commands in screen [def OFF]   
";
}

my @species=split(/\,/,$species_array);
@species=sort(@species); # sort alphanumerically
my ($in,$of,$l, $last,$p, $jid, $tx);
my @pt;
my $temp_sp = join("_",@species);
my $project_dir = "PROJECT_".$temp_sp;


### Prepares folders and submissions
for my $i (0..$#species){
    my $sp1 = $species[$i];

    for my $j (0..$#species){
	my $sp2 = $species[$j];

	if ($i < $j){ # only does Sp1_Sp2 (not Sp2_Sp1)
	    my $pair_folder=$dir."/".$project_dir."/".$sp1."_".$sp2;
	    my $inf=$pair_folder."/Final_aln_scores_".$sp1."_".$sp2.".txt";
	    if (defined ($verbose)){
			print "submitjob long -l h_rt=10:00:00,virtual_free=10G perl $bin/get_score_by_exon.pl $inf $pair_folder\n";
	    }
	    if (defined ($submit_jobs)){
			system "submitjob long -l h_rt=10:00:00,virtual_free=10G perl $bin/get_score_by_exon.pl $inf $pair_folder";
	    }
                                                                   
	}
   }

}


