#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

### AIMS OF SCRIPT
#* submit jobs to map and scores exons for each sub-file of clusters between pair of species proteins

### FOLDER STRUCTURE
#* Root:                           $dir/
#* Common to different projects:   $dir/EXONS_DB/ with all the general info per species
#* Project specific:               $dir/PROJECT_Sp1_Sp2_Sp3/ => made alphabetically [will avoid bugs]
#   * by pairs:                    $dir/PROJECT_Sp1_Sp2_Sp3/Sp1_Sp2/ etc
#   * for gene cluster parts:      $dir/PROJECT_Sp1_Sp2_Sp3/GENE_CLUSTERS/          


#### DECLARATION AND GENERATION OF VARIABLES
my $f_gene_cluster; # file with the gene clusters => it may need additional and specific rules...
my $f_extra_cluster; # a different cluster file for some species (hard to generalize); eventually allow it to be comma separated
my $species_array; ##species separated by comma 
my $submit_aln; # whether to submit the jobs
#my $root = "/users/mirimia/ymarquez/ANALYSIS/ORTOLOGUES/ORTH_PIPE/EXON_ORTH"; # default (folder where files and program of ORTH_PIPE are located)
my $bin;
my $dir;
my $verbose; # prints steps
my $help;

### Get options:
Getopt::Long::Configure("no_auto_abbrev");
GetOptions(  "gene_cluster=s" => \$f_gene_cluster,
	     "extra_cluster=s" => \$f_extra_cluster,
	     "sp=s" => \$species_array,
	     "submit_aln" => \$submit_aln,
	     "dir=s" => \$dir,
	     "bin=s" => \$bin,
	     "verbose" => \$verbose,
	     "help" => \$help
    );

### Help
if (!defined ($dir) || !defined($bin) || !defined($species_array) || !defined($species_array) || defined ($help)){
    die "
Usage: Run_map_each_exon_prot.pl -sp Sp1,Sp2,Sp3 --dir PATH project  --gene_cluster File_gene_clusters --bin PATH script files [options]

OPTIONS
    --clean                   Force removal of previous data [def OFF]
    --submit_aln              Submit the jobs of for scoring ex/introns [def OFF]
    --dir PATH                Root directory where project will be created.
    --gene_cluster	      Gene cluster file
    --bin PATH               bin directory where the scripts are located
    --verbose                 Prints commands in screen [def OFF]     
    --extra_cluster  FILE     Additional file for certain species [def OFF]

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
    $f_exint{$sp}="$dir/EXONS_DB/$sp/$sp"."_annot_fake.exint";
    $f_protIDs{$sp}="$dir/EXONS_DB/$sp/$sp"."_annot_fake_protein_ids_exons.txt";
    $f_exon_pos{$sp}="$dir/EXONS_DB/$sp/$sp"."_annot_fake_protein_ids_exons_pos.txt";
    $f_intron_pos{$sp}="$dir/EXONS_DB/$sp/$sp"."_annot_fake_protein_ids_intron_pos_CDS.txt";

    die "Cannot find exint file for $sp\n" unless (-e $f_exint{$sp});
    die "Cannot find prot IDs file for $sp\n" unless (-e $f_protIDs{$sp});
    die "Cannot find exon pos file for $sp\n" unless (-e $f_exon_pos{$sp});
    die "Cannot find intron pos file for $sp\n" unless (-e $f_intron_pos{$sp});
}

### Setting up and creating PROJECT directory
my $temp_sp = join("_",@species);
my $project_dir = "PROJECT_".$temp_sp;
my $gcl=0;
my %N_parts;
if (-e "$dir/$project_dir/GENE_CLUSTERS"){ 
	$gcl=1; 
	my $parts=`ls $dir/$project_dir/GENE_CLUSTERS/*part* | wc -l`;
	$N_parts{$f_gene_cluster}=$parts;

} else { print "Wrong path to gene cluster parts, check your working directory\n\n"; exit; }

print "\n$f_gene_cluster parts: $N_parts{$f_gene_cluster}\n" if defined ($verbose);

### Prepares folders and submissions
for my $i (0..$#species){
    my $sp1 = $species[$i];

    for my $j (0..$#species){
	my $sp2 = $species[$j];

	if ($i < $j){ # only does Sp1_Sp2 (not Sp2_Sp1)
	    my $pair_folder=$sp1."_".$sp2;	    
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

		my $tb = $bin;
		$tb=~s/bin/files/;
		my $bl62 = "$tb/blosum62.txt";
		my $of=$dir."/".$project_dir."/".$pair_folder;
		### Submits for exons and introns
		if (defined ($verbose)){
		    print "\nsubmitjob long -l h_rt=300:00:00,virtual_free=30G ".
			"$bin/map_each_exon_prot.pl $sp1 $sp2 $cl_part ".
			"$f_protIDs{$sp1} $f_protIDs{$sp2} $f_exon_pos{$sp1} $f_exon_pos{$sp2} $f_intron_pos{$sp1} $f_intron_pos{$sp2} $f_exint{$sp1} $f_exint{$sp2} $part $bin $bl62 $of\n";
		}
		system "submitjob long -l h_rt=300:00:00,virtual_free=30G ". # job conditions
		    "perl $bin/map_each_exon_prot.pl $sp1 $sp2 $cl_part ". # species and gene_cluster
		    "$f_protIDs{$sp1} $f_protIDs{$sp2} $f_exon_pos{$sp1} $f_exon_pos{$sp2} $f_intron_pos{$sp1} $f_intron_pos{$sp2} $f_exint{$sp1} $f_exint{$sp2} $part $bin $bl62 $of\n" if $submit_aln; # input files and outputs
	    }
	}	
    }
}
    



