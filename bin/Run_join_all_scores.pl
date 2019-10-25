#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

### AIMS OF SCRIPT
#* to prepare and submit the jobs to join all the scores resulting from the aligments of exons and introns in each pair of species folder

### FOLDER STRUCTURE
#* Root:                           $dir/
#* Common to different projects:   $dir/EXONS_DB/ with all the general info per species
#* Project specific:               $dir/PROJECT_Sp1_Sp2_Sp3/ => made alphabetically [will avoid bugs]
#   * by pairs:                    $dir/PROJECT_Sp1_Sp2_Sp3/Sp1_Sp2/ etc
#   * for gene cluster parts:      $dir/PROJECT_Sp1_Sp2_Sp3/GENE_CLUSTERS/ 

#### DECLARATION AND GENERATION OF VARIABLES
my $species_array; ##species separated by comma 
my $submit_jobs; # whether to submit the jobss
my $dir;
my $bin;
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


if (!defined ($dir) || !defined ($bin) || !defined($species_array) || defined ($help)){
    die "
Usage: Run_join_all_scores.pl -sp Sp1,Sp2,Sp3 --dir PATH --bin PATH [options]

OPTIONS
    --clean                   Force removal of previous data [def OFF]
    --submit_jobs             Submit the jobs of for scoring ex/introns [def OFF]
    --dir PATH                Directory where project was created.
    --bin PATH                Bin directory where scripts are located
    --verbose                 Prints commands in screen [def OFF]   
    --extra_cluster  FILE     Additional file for certain species [def OFF]
";
}

my @species=split(/\,/,$species_array);
@species=sort(@species); # sort alphanumerically
my %f_exint;
foreach my $sp (@species){ # technically speaking, the downstream scripts should handle with just providing Sp1 and Sp2
    $f_exint{$sp}="$dir/EXONS_DB/$sp/$sp"."_annot_fake.exint";
}

my $tb=$bin;
$tb=~s/bin/files/;
my $bl62 = "$tb/blosum62.txt";
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
	    
	    system `cd $pair_folder`;
	    if (defined ($verbose)){
			print "cd $pair_folder\n";
			print "submitjob perl $bin/join_score_files.pl $pair_folder\n";
	    }
			if (defined ($submit_jobs)){
				system "submitjob perl $bin/join_score_files.pl $pair_folder";
			}
                                                                   
	}
   }

}





