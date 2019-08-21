#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

### AIMS OF SCRIPT
#* Prepare and submit the jobs to realign exons that result for the first round of alignments for each pair of species folder 

### FOLDER STRUCTURE
#* Root:                           $dir/
#* Common to different projects:   $dir/EXONS_DB/ with all the general info per species
#* Project specific:               $dir/PROJECT_Sp1_Sp2_Sp3/ => made alphabetically [will avoid bugs]
#   * by pairs:                    $dir/PROJECT_Sp1_Sp2_Sp3/Sp1_Sp2/ etc


#### DECLARATION AND GENERATION OF VARIABLES
my $species_array; ##species separated by comma 
my $submit_aln; # whether to submit the jobs
my $bin; # default (folder where files and program of ORTH_PIPE are located)
my $dir;
my $verbose; # prints steps
my $help;

### Get options:
Getopt::Long::Configure("no_auto_abbrev");
GetOptions(  "sp=s" => \$species_array,
	     "submit_aln" => \$submit_aln,
	     "dir=s" => \$dir,
	     "bin=s" => \$bin,
	     "verbose" => \$verbose,
	     "help" => \$help
    );


if (!defined ($dir) || !defined($species_array) || !defined($bin)|| defined ($help)){
    die "
Usage: Run_realn_exons.pl -sp Sp1,Sp2,Sp3 --dir PATH [options]

OPTIONS
    --clean                   Force removal of previous data [def OFF]
    --submit_aln              Submit the jobs of for scoring ex/introns [def OFF]
    --dir PATH                Root directory where project will be created. This directory has to have the EXONS_DB folder with the output obtained in Module I for each species.
    --bin PATH                Bin directory where the exon orthology pipeline is located                 
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
	    
	    my $parts=`ls $pair_folder/exons_to_realign_part_*`;
	    my @files=split(/\n/,$parts);

	    for ($l=0; $l<scalar(@files); $l++){
			$in=$files[$l];
			@pt=split(/\_/,$in);
			$last=scalar(@pt)-1;
			$of=$pair_folder."/realigned_exons_".$sp1."_".$sp2."_".$pt[$last];
			$p=$pt[$last];
			$tx=".txt";
			$p=~s/$tx//;
			$jid="R".$p."_".$sp1.",".$sp2;
			if (defined ($verbose)){
				print "system \"submitjob long -l h_rt=70:00:00,virtual_free=20G perl $bin/Realign_target_exons_by_part.pl $sp1 $sp2 $in $f_exint{$sp1} $f_exint{$sp2} $p $of $pair_folder $bl62 $bin\";\n";

			}
			if (defined ($submit_aln)){
				system "submitjob long -l h_rt=70:00:00,virtual_free=20G perl $bin/Realign_target_exons_by_part.pl $sp1 $sp2 $in $f_exint{$sp1} $f_exint{$sp2} $p $of $pair_folder $bl62 $bin\n";
			}
	    }		

	}
   }

}





