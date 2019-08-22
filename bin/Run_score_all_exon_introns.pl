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
my $submit_jobs; # whether to submit the jobs
my $bin; # default (folder where files and program of ORTH_PIPE are located)
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
Usage: Run_score_all_exon_introns.pl -sp Sp1,Sp2,Sp3 --dir PATH [options]

OPTIONS
    --clean                   Force removal of previous data [def OFF]
    --submit_jobs             Submit the jobs of for scoring ex/introns [def OFF]
    --dir PATH                Directory where project will be created. 
    --bin PATH                Bin directory where files and scripts of ORTH_PIPE ARE LOCATED
    --verbose                 Prints commands in screen [def OFF]   
    --extra_cluster  FILE     Additional file for certain species [def OFF]
";
}

my @species=split(/\,/,$species_array);
@species=sort(@species); # sort alphanumerically
my ($in,$of,$l, $last,$p, $jid, $tx);
my @pt;
my $temp_sp = join("_",@species);
my $project_dir = "PROJECT_".$temp_sp;

### Stores the relevant files from EXONS_DB/Sp
my %f_exint;
my %int_pos;
foreach my $sp (@species){ # technically speaking, the downstream scripts should handle with just providing Sp1 and Sp2
    $f_exint{$sp}="$dir/EXONS_DB/$sp/$sp"."_annot_fake.exint";
    $int_pos{$sp}="$dir/EXONS_DB/$sp/$sp"."_protein_ids_intron_pos_CDS.txt";
}


### Prepares folders and submissions
for my $i (0..$#species){
    my $sp1 = $species[$i];

    for my $j (0..$#species){
	my $sp2 = $species[$j];

	if ($i < $j){ # only does Sp1_Sp2 (not Sp2_Sp1)
	    my $pair_folder=$dir."/".$project_dir."/".$sp1."_".$sp2;
	    my $in1=$pair_folder."/Aligned_proteins.txt";
	    my $in2=$pair_folder."/Best_scores_pair_exons.txt";
	    my $in3=$pair_folder."/Score_all_introns.txt";
	    my $out=$pair_folder."/Final_aln_scores_".$sp1."_".$sp2.".txt";
	    my $e1=$f_exint{$sp1};
	    my $e2=$f_exint{$sp2};
	    my $if1=$int_pos{$sp1};
	    my $if2=$int_pos{$sp2};
	    if (defined ($verbose)){
			print "submitjob long -l h_rt=3:00:00,virtual_free=10G perl $bin/get_scores_exons_introns.pl $sp1 $sp2 $in1 $in2 $in3 $e1 $e2 $if1 $if2 $out\n";
	    }
	    if (defined ($submit_jobs)){
			system "submitjob long -l h_rt=3:00:00,virtual_free=10G perl $bin/get_scores_exons_introns.pl $sp1 $sp2 $in1 $in2 $in3 $e1 $e2 $if1 $if2 $out";
	    }
                                                                   
	}
   }

}





