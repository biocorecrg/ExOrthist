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
my $join; # whether to join the files
my $bin;
my $dir;
my $verbose; # prints steps
my $help;

### Get options:
Getopt::Long::Configure("no_auto_abbrev");
GetOptions(  "sp=s" => \$species_array,
	          "join" => \$join,
	          "dir=s" => \$dir,
	          "bin=s" => \$bin,
	          "verbose" => \$verbose,
	          "help" => \$help
    );


if (!defined ($dir) || !defined ($bin) || !defined($species_array) || defined ($help)){
    die "
Usage: Join_best_scores_files.pl -sp Sp1,Sp2,Sp3 --dir PATH [options]

OPTIONS
    --clean                   Force removal of previous data [def OFF]
    --join                      Do joining of the files [def OFF]
    --dir PATH                Root directory where project will be created 
    --bin PATH                Root directory where scripts are located
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
	    my $inf=$pair_folder."/Best_score_hits_exons.txt ";
	    my $out_folder=$dir."/".$project_dir."/Best_score_hits_exons.txt";
	    if (defined ($verbose)){
		print "cat $inf >> $out_folder\n";
	    }
	    if (defined ($join)){
		system "cat $inf >> $out_folder\n";
	    }
                                                                   
	}
    }

}
