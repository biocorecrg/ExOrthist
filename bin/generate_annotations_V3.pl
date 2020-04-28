#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Cwd qw(abs_path cwd);

#Declaration of variables
my $gtf_file;
my $genome_file;
my $exons_db_folder="./";
#my $species_string;
my $species_id;
my $verboseFlag=1;
my $help;
my $add_exons="NA";

#Arguments

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(        "GTF=s" => \$gtf_file,
		              "G=s" => \$genome_file,
		              "EX_DB=s" => \$exons_db_folder,
		              "sp=s" => \$species_id,
		              "add_exons=s" => \$add_exons, #before $vastdb_refs
		              "verbose=s" => \$verboseFlag,
		              "h" => \$help,
		              "help" => \$help
    );


sub verbPrint {
    my $verbMsg = shift;
    unless ($verboseFlag == 0 || $verboseFlag eq "F" || $verboseFlag eq "FALSE") {
	chomp($verbMsg);
	print STDERR "[exorter annotation]: $verbMsg\n";
    }
}
if (!defined $genome_file || !defined $gtf_file || !defined $species_id || defined $help){
    die "\nUsage: generate_annotations.pl -GTF path_to_gtfs/ -G path_to_genomes/ -sp Sp1 [-EX_DB path_to_EXONS_DB/ -add_exons exon_file ]
Script that creates all annotation files needed for the second module of the pipeline
COMPULSORY
     -GTF              Path where GTFs are stored (they should be named Sp1_annot.gtf)
     -G                Path where gDNAs are stored (they should be named Sp1_gDNA.fasta)
     -sp Sp1           Species id.
OPTIONAL
     -add_exons        File with additional exons to get their orthology
                           If a species is missing the reference file, \"NA\" should be provided
     -EX_DB            Path to EXONS_DB/ folder (default ./; i.e. working directory)
                           If it does not exit, it will create a EXONS_DB/ folder in the working directory
     -verbose T/F      Verbose (default TRUE) 
     -h/--help         This help message.
";    
}

### setting EXONS_DB
system "mkdir $exons_db_folder" unless (-e $exons_db_folder);
my $full_path_exons_db = abs_path($exons_db_folder);
verbPrint("EXONS_DB path set to $full_path_exons_db\n");
#my %VASTDB_files;
my $species = $species_id;
### loops for each species (Loop 1: previous get_ref_proteins.pl)
#foreach my $species (@SPECIES){
 # In this first loop, it checks if the folders for each species exist or create them
system "mkdir $exons_db_folder/$species" unless (-e "$exons_db_folder/$species");
my (%seq);
my $chr;
my @g;
if ($genome_file =~ /.gz$/) {
    open(GENOME, "gunzip -c $genome_file |") || die "t open pipe to $genome_file";
}
else {
    open(GENOME, $genome_file) || die "t open $genome_file";
}    #open (GENOME, $genome_file) || die "Cannot open $genome_file for $species (loop 1)\n";
verbPrint("Parsing gDNA file for $species\n");
while (<GENOME>){
    chomp($_);
    if ($_=~/\>/){ $_=~s/\>//; @g=split(/\s+/,$_); $chr=$g[0]; }
    else { $seq{$chr}.=$_;  }    
}
close (GENOME);
verbPrint("Generating exint file for $species\n");
##Opening GTF file##
# Format GTF:
#Scaffold1348protein_codingexFon162899162997.+.gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "5";
#Scaffold1348protein_codingCDS162899162997.+.gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "5"; protein_id "WHL22.100348.0";

my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
my (%pid, %header, %ntseq);
#my $b;
my (@l, @l1, @line, @l2, @l3, @s);
my %rseq;
my ($res, $size, $tmpseq);
my %fex; ##saving the phase of the first exon
if ($gtf_file =~ /.gz$/) {
    open(GTF, "gunzip -c $gtf_file |") || die "t open pipe to $gtf_file";
}
else {
    open(GTF, $gtf_file) || die "t open $gtf_file";
}
while (<GTF>){
    chomp($_); 
    @line=split(/\t/,$_);
    if ($seq{$line[0]}){
	if ($line[2] eq "CDS"){
	    @l1=split(/\"/,$line[8]);
	    $gid=$l1[1]; ##change here
	        #print "GID: $gid\n";
	    if ($line[8]=~/protein_id/){
		$tmp="protein_id";
	    } 
	    else { $tmp="transcript_id"; }
	    @l2=split(/$tmp/,$line[8]);
	    @l3=split(/\"/,$l2[1]);
	    $prot=$l3[1];
	    if (!$pid{$prot}){
		$m=0;
		if ($line[7] ne "."){  $phase=$line[7]; $fex{$prot}=$phase; $m=$phase; } ##
		else { 
		    $phase=0; 
		    $fex{$prot}=$phase;
		}
		$nuc=0;
		$size=$line[4]-$line[3]+1;
		$tmpseq=substr($seq{$line[0]},($line[3]-1),($size));
		$size=$size-$m;
		if ($line[6] eq "-"){  $tmpseq=~tr/ATCG/TAGC/; $tmpseq=reverse($tmpseq); }
		$ntseq{$prot}=$tmpseq;
		$rseq{$prot}=$tmpseq;
		$res=$size%3;
		$nuc=$nuc+$size;
		$aa=int(1+$nuc/3); 
		$pid{$prot}=1;
		$header{$prot}="$prot|$gid ";
	    }
	    else {
		$header{$prot}.=" $aa.";
		if ($line[7] ne "."){  $phase=$line[7]; }
		else {  
		    if ($res==0){ $phase=0; } elsif ($res==1 ){ $phase=2; } elsif ($res==2) { $phase=1; } 
		}
		$size=($line[4]-$line[3])+1;
		$tmpseq=substr($seq{$line[0]},($line[3]-1),($size));
		if ($line[6] eq "-"){  $tmpseq=~tr/ATCG/TAGC/; $tmpseq=reverse($tmpseq); }
		$ntseq{$prot}.=$tmpseq;
		$rseq{$prot}.="|".$tmpseq;
		$nuc=$nuc+$size;
		$size=$size+$res;
		$res=$size%3;
		$aa=int(1+$nuc/3); ##adding the residues that were left
		$header{$prot}.="$phase";
	    }
	}
    }
}
close (GTF);
    
###Getting protein sequence####
my %codons;
$codons{"TTT"}="F";    $codons{"TTC"}="F";    $codons{"TTA"}="L";    $codons{"TTG"}="L";
$codons{"CTT"}="L";    $codons{"CTC"}="L";    $codons{"CTA"}="L";    $codons{"CTG"}="L";
$codons{"TCT"}="S";    $codons{"TCC"}="S";    $codons{"TCA"}="S";    $codons{"TCG"}="S";
$codons{"TAT"}="Y";    $codons{"TAC"}="Y";    $codons{"TAA"}="stop"; $codons{"TAG"}="stop";
$codons{"TGT"}="C";    $codons{"TGC"}="C";    $codons{"TGA"}="stop"; $codons{"TGG"}="W";    
$codons{"CCT"}="P";    $codons{"CCC"}="P";    $codons{"CCA"}="P";    $codons{"CCG"}="P";    
$codons{"CAT"}="H";    $codons{"CAC"}="H";    $codons{"CAA"}="Q";    $codons{"CAG"}="Q";    
$codons{"CGT"}="R";    $codons{"CGC"}="R";    $codons{"CGA"}="R";    $codons{"CGG"}="R";    
$codons{"AGT"}="S";    $codons{"AGC"}="S";    $codons{"AGA"}="R";    $codons{"AGG"}="R";    
$codons{"ATT"}="I";    $codons{"ATC"}="I";    $codons{"ATA"}="I";    $codons{"ATG"}="M";
$codons{"ACT"}="T";    $codons{"ACC"}="T";    $codons{"ACA"}="T";    $codons{"ACG"}="T";
$codons{"AAT"}="N";    $codons{"AAC"}="N";    $codons{"AAA"}="K";    $codons{"AAG"}="K";    
$codons{"GTT"}="V";    $codons{"GTC"}="V";    $codons{"GTA"}="V";    $codons{"GTG"}="V";
$codons{"GCT"}="A";    $codons{"GCC"}="A";    $codons{"GCA"}="A";    $codons{"GCG"}="A";    
$codons{"GAT"}="D";    $codons{"GAC"}="D";    $codons{"GAA"}="E";    $codons{"GAG"}="E";
$codons{"GGT"}="G";    $codons{"GGC"}="G";    $codons{"GGA"}="G";    $codons{"GGG"}="G";

my @keys=keys(%ntseq);
my $el;
my $p;
my ($l, $triplet, $protein, $nseq, $name);
my $c=0;
### output files (former -F 1 run)
my $exint_output = "$exons_db_folder/$species/$species.exint";
open (EXINT_OUT, ">$exint_output") || die "Cannot open $exint_output (loop 1)\n";
foreach $el(@keys){ 
    $nseq=$ntseq{$el};
    $name=$header{$el};
    $p=$fex{$el}; ##checking phase of the first exon
    $protein="";
    $c=0;
    for($l=$p; $l<length($nseq); $l=$l+3){
	$triplet=substr($nseq,$l,3);
	if ($c==0){
	    if ($codons{$triplet}){
		if ($l>=(length($nseq)-3) && $codons{$triplet} eq "stop"){
		    $c=1;
		}
		elsif ($codons{$triplet} eq "stop"){ $c=1; } #$protein.=$codons{$triplet}; } #$c=1; }
		else {  $protein.=$codons{$triplet};   }
	    }
	    else {
		$protein.="X";
	    }
	}
    }
    print EXINT_OUT ">$name\n$protein\n";
}
close (EXINT_OUT);
verbPrint("Generating size and references files for $species\n");
### output files (former -F 2 run) 
my $size_output = "$exons_db_folder/$species/$species"."_prot_sizes.txt";
open (SIZE_OUT, ">$size_output") || die "Cannot open $size_output (loop 1)\n";
open (EXINT_IN, $exint_output) || die "Cannot open $exint_output (loop 1)\n";
my ($exint_prot_name,$exint_gene_name, $header);
my %length_gene;
my %longest_prot;
my %exint_data;
while (<EXINT_IN>){
    chomp($_);
    if ($_=~/\>/){ 
	$header = $_;
	my @line=split(/\s+/,$_); # >protID|geneID
	$line[0]=~s/\>//; 
	($exint_prot_name,$exint_gene_name) = split(/\|/,$line[0]);
    }
    else { 
	my $s2=length($_);
	$exint_data{$exint_prot_name}="$header\n$_";
	    
	print SIZE_OUT "$exint_prot_name\t$exint_gene_name\t$s2\n";
	    
	if (!defined $length_gene{$exint_gene_name}){
	    $length_gene{$exint_gene_name}=$s2;
	    $longest_prot{$exint_gene_name} = $exint_prot_name;
	}
	elsif ($s2 > $length_gene{$exint_gene_name}){
	    $length_gene{$exint_gene_name}=$s2;
	    $longest_prot{$exint_gene_name} = $exint_prot_name;
	}
    }
}
close EXINT_IN;
close SIZE_OUT;
my $longest_prot_output = "$exons_db_folder/$species/$species"."_ref_proteins.txt";
my $longest_exint_output = "$exons_db_folder/$species/$species"."_ref_proteins.exint";
open (L_P_OUTPUT, ">$longest_prot_output") || die "Cannot open $longest_prot_output (loop 1)\n";
open (L_E_OUTPUT, ">$longest_exint_output") || die "Cannot open $longest_exint_output (loop 1)\n";
#    print L_P_OUTPUT "GeneID\tProteinID\n"; # maintain without header, as the original
foreach my $temp_exint (sort keys %longest_prot){
    print L_P_OUTPUT "$temp_exint\t$longest_prot{$temp_exint}\n";
    print L_E_OUTPUT "$exint_data{$longest_prot{$temp_exint}}\n";
}
close L_P_OUTPUT;
close L_E_OUTPUT;

###IntroduceEXON_to_GTF
##Introducing missing exons in GTF
if ($add_exons ne "NA"){##if added exons ne "NA" --> if additional exons are provided by the user 
    verbPrint("Introducing additional exons for $species\n");
#ENSG00000029534HsaEX0004110chr8:41758036-41758137chr8:41797512-41797622chr8:41733971-41734069chr8:41797512,41758036-41758137,41734069:-
#### Get the exons and their info
    my (@t);
    my ($g, $ev,$i,$f,$coA,$coB);
    my (%Ev_Gene,%C1_ref,%A_ref,%C2_ref,%coA_ev,%co_ev,%coB_ev,%C1_inc,%C2_inc);

    open (EXONS, $add_exons) || die "*** DIE: Cannot open $add_exons\n";
    <EXONS>;
    while (<EXONS>){
	chomp;
	@t=split(/\t/);
	$g=$t[0];
	$ev=$t[1];
	$Ev_Gene{$ev}=$g;
	$C1_ref{$ev}=$t[3];
	$A_ref{$ev}=$t[2];
	$C2_ref{$ev}=$t[4];
	$C1_inc{$ev}=$t[3]; ##same as C1 reference
	$C2_inc{$ev}=$t[4]; ##same as C2 reference 	
	$co_ev{$t[2]}=$ev;
	($chr,$i,$f)=$t[2]=~/(.+?)\:(.+?)\-(.+)/;
	$coA="$chr:$i";
	$coB="$chr:$f";
	$coA_ev{$coA}=$ev;
	$coB_ev{$coB}=$ev;
   #print "$ev\t$t[2]\t$t[3]\t$t[4]\t$coA\t$coB\n";
    }
    close EXONS;

    my ($tr, $co,$exon_number,$g_2,$tr_2,);
    my (%gene_name,%annotated,%partial,%array_tr_co,%array_tr_coA,%array_tr_coB,%index_coA,%index_coB,%TR,%index_co,%tr_co,%tr_coA,%tr_coB,%str,%done);
    my (%cds_lines,%cds_ini,%cds_end,%offset,%start_lines,%stop_lines,%tr_lines);
#### Starts checking for annotation
    open (GTF, $gtf_file) || die "*** DIE: Cannot open GTF\n";
    while (<GTF>){
	chomp;
	@t=split(/\t/);    
	if ($t[2] eq "exon"){
	    ($g)=$t[8]=~/gene_id \"(.+?)\"/;
	    ($tr)=$t[8]=~/transcript_id \"(.+?)\"/;
	    ($name)=$t[8]=~/gene_name \"(.+?)\"/;
	    $gene_name{$g}=$name;
	    
	    $co="$t[0]:$t[3]-$t[4]";
	    $coA="$t[0]:$t[3]";
	    $coB="$t[0]:$t[4]";
	    #print "$co\t$coA\t$coB\n";
	    $annotated{$co_ev{$co}}=1 if $co_ev{$co} && $Ev_Gene{$co_ev{$co}} eq $g; # keeps track of annotated events
	    ## any event not here is a non-annotated event (although subversions may be in GTFs
	    $partial{$coA_ev{$coA}}=1 if $coA_ev{$coA} && $Ev_Gene{$coA_ev{$coA}} eq $g;
	    $partial{$coB_ev{$coB}}=1 if $coB_ev{$coB} && $Ev_Gene{$coB_ev{$coB}} eq $g;
	    
	    ### now stores data for C1 and C2s
	    push(@{$array_tr_co{$tr}},$co); # for each tr, the whole array of exons
	    $index_co{$tr}{$co}=$#{$array_tr_co{$tr}}; # keeps the index of the exon; C1-C2 MUST be index+1 of the other
	    
	    push(@{$array_tr_coA{$tr}},$coA); # for each tr, the whole array of exons
	    $index_coA{$tr}{$coA}=$#{$array_tr_coA{$tr}}; # keeps the index of the exon; C1-C2 MUST be index+1 of the other
	    push(@{$array_tr_coB{$tr}},$coB); # for each tr, the whole array of exons
	    $index_coB{$tr}{$coB}=$#{$array_tr_coB{$tr}}; # keeps the index of the exon; C1-C2 MUST be index+1 of the other
	    #print "##$index_coA{$tr}{$coA}\t$index_coB{$tr}{$coB}\n";
#push(@{$TR{$g}},$tr); # for each g, gets all tr
	    $TR{$g}{$tr}=1;
	    $tr_co{$tr}{$co}=1; # co exists in that tr
	    $tr_coA{$tr}{$coA}=1; # coA exists in that tr
	    $tr_coB{$tr}{$coB}=1; # coB exists in that tr
	    
	    push(@{$tr_lines{$tr}},$_); # all lines of a transcript. Index 0 should be exon 1
	    $exon_number=$#{$tr_lines{$tr}};

	    $str{$g}=$t[6];
	    $done{$tr}=1;
	}
	elsif ($t[2] eq "CDS"){ # In principle, CDS comes always after exon, and then start or stop
	    ($g_2)=$t[8]=~/gene_id \"(.+?)\"/;
	    ($tr_2)=$t[8]=~/transcript_id \"(.+?)\"/;
	    $cds_lines{$tr_2}[$exon_number]=$_;
	    $cds_ini{$tr_2}[$exon_number]=$t[3];
	    $cds_end{$tr_2}[$exon_number]=$t[4];
	    $offset{$tr_2}[$exon_number]=$t[7];
	    die "$tr_2 ne $tr\n" if $tr_2 ne $tr; # ie. no match between exon and next line CDS
	}
	elsif ($t[2] eq "start_codon"){
	    ($g_2)=$t[8]=~/gene_id \"(.+?)\"/;
	    ($tr_2)=$t[8]=~/transcript_id \"(.+?)\"/;
	    $start_lines{$tr_2}[$exon_number]=$_;
	    die "$tr_2 ne $tr\n" if $tr_2 ne $tr; # ie. no match between exon and next line CDS
	}
	elsif ($t[2] eq "stop_codon"){
	    ($g_2)=$t[8]=~/gene_id \"(.+?)\"/;
	    ($tr_2)=$t[8]=~/transcript_id \"(.+?)\"/;
	    $stop_lines{$tr_2}[$exon_number]=$_;
	    die "$tr_2 ne $tr\n" if $tr_2 ne $tr; # ie. no match between exon and next line CDS
	}
    }
    close GTF;

#
### v3 => vB to avoid confusions with Ensembl versions (10/02/17)
    open (O, ">$exons_db_folder/$species/FakeTranscripts-$species-vB.gtf");
    open (LOG, ">$exons_db_folder/$species/LOG_FakeTranscripts-$species-vB.tab");
    print LOG "EVENT\tGENE_ID\tEXON\tFAKE_TR\tTYPE\tSOLUTION\n";
    open (LOG2, ">$exons_db_folder/$species/STATS_FakeTranscripts-$species-vB.tab");

#### Now loops through every exon, and focuses on the non-annotated
# Rules:
# 1) if C1inc and C2inc are present OK.
# ...
# In some cases, the exact exon is wrongly annotated, but there is no skipping form
    my ($tally_annotated,$selected_tr,$C1_inc,$C2_inc,$C1_ref,$C2_ref,$C1i,$C1r,$C1_incA,$C1_refA,$C1f,$C2f,$C2_incA,$C2_refA,$C2i,$C2r,$Ai,$Af,$coA_A,$coB_A,$temp_tr,$exons_tr,$selected_tr,$ref_tr,$OK,$t_C2,$exons_sel,$t_C1,$tally_non_annot_hit,$rescue_type,$final_tr,$exN,$out_frame,$started,$finished,$comment,$CDS_seq,$CDS_ini,$CDS_end);
    my($CDS_end,$CDS_ini,$CDS_seq,$STOP_detected,$bit,$co_Af_CDS,$co_Ai_CDS,$codon,$comment,$exN,$extra,$final_tr,$finished,$first_offset,$full_line,$full_line2,$full_line3,$full_line4);
    my ($le_CDS_ex,$length_A,$lineN,$loop,$new_line,$new_line2,$new_line3,$new_offset,$offset_A,$out_frame,$rescue_type,$seq_bit,$seq_bit_test,$started,$stop_co_f,$stop_co_i,$str,$tally_non_annot_hit,$tally_non_annot_no_hit,$tally_partial,$type);
    my (@t_line,@t_line2,@t_line3,@t_line4,@t_line_C2);
    my (%C1_accepted_index,%tally_solutions);
    foreach $ev (sort keys %A_ref){
	$g=$Ev_Gene{$ev};
	if ($annotated{$ev}){
	    $tally_annotated++;
	}
	else {
	    $selected_tr="";
	    $C1_inc=$C1_inc{$ev};
	    $C2_inc=$C2_inc{$ev};
	    $C1_ref=$C1_ref{$ev};
	    $C2_ref=$C2_ref{$ev};
	    ### only the C1do and C2ac
	    ($chr,$C1i,$C1f)=$C1_inc=~/(.+?)\:(.+?)\-(.+)/;
	    $C1_incA="$chr:$C1f" if $str{$g} eq "+";
	    $C1_incA="$chr:$C1i" if $str{$g} eq "-";
	    ($chr,$C1i,$C1f)=$C1_ref=~/(.+?)\:(.+?)\-(.+)/;
	    $C1_refA="$chr:$C1f" if $str{$g} eq "+";
	    $C1_refA="$chr:$C1i" if $str{$g} eq "-";
	    ($chr,$C2i,$C2f)=$C2_inc=~/(.+?)\:(.+?)\-(.+)/;
	    $C2_incA="$chr:$C2f" if $str{$g} eq "-";
	    $C2_incA="$chr:$C2i" if $str{$g} eq "+";
	    ($chr,$C2i,$C2f)=$C2_ref=~/(.+?)\:(.+?)\-(.+)/;
	    $C2_refA="$chr:$C2f" if $str{$g} eq "-";
	    $C2_refA="$chr:$C2i" if $str{$g} eq "+";
	    
	    ### exon coordinates
	    ($Ai,$Af)=$A_ref{$ev}=~/\:(.+?)\-(.+)/;
$coA_A="$chr:$Ai";
$coB_A="$chr:$Af";
#print "#$ev\t$C1_inc\t$C2_inc\t$C1_ref\t$C2_ref\t$C1_incA\t$C1_refA\t$C2_incA\t$C2_refA\t$coA_A\t$coB_A#\n";
	    foreach $tr (sort keys %{$TR{$g}}){
		next if $tr_coA{$tr}{$coA_A} || $tr_coB{$tr}{$coB_A}; # if the transcript contains any junction of the A exon

		if ($tr_co{$tr}{$C1_inc} && $tr_co{$tr}{$C2_inc} && $index_co{$tr}{$C1_inc}==$index_co{$tr}{$C2_inc}-1){
		    ($temp_tr)=$selected_tr=~/1\=(.+)/; # it may be empty
		    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
		    $exons_tr=$#{$cds_lines{$tr}}+1;
		    $selected_tr="1=$tr" unless ($selected_tr eq "1=$ref_tr" || $exons_sel>=$exons_tr);
		    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_inc};
		}
		elsif ($tr_co{$tr}{$C1_inc} && $tr_co{$tr}{$C2_ref} && $index_co{$tr}{$C1_inc}==$index_co{$tr}{$C2_ref}-1){
		    ($temp_tr)=$selected_tr=~/[12]\=(.+)/; # it may be empty
		    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
		    $exons_tr=$#{$cds_lines{$tr}}+1;
		    $selected_tr="2=$tr" unless ($selected_tr eq "2=$ref_tr" || $selected_tr=~/1\=/ || $exons_sel>=$exons_tr);
		    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_inc};
		}
		elsif ($tr_co{$tr}{$C1_ref} && $tr_co{$tr}{$C2_inc} && $index_co{$tr}{$C1_ref}==$index_co{$tr}{$C2_inc}-1){
		    ($temp_tr)=$selected_tr=~/[123]\=(.+)/; # it may be empty
		    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
		    $exons_tr=$#{$cds_lines{$tr}}+1;
		    $selected_tr="3=$tr" unless ($selected_tr eq "3=$ref_tr" || $selected_tr=~/[12]\=/ || $exons_sel>=$exons_tr);
		    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_ref};
		}
		elsif ($tr_co{$tr}{$C1_ref} && $tr_co{$tr}{$C2_ref} && $index_co{$tr}{$C1_ref}==$index_co{$tr}{$C2_ref}-1){
		    ($temp_tr)=$selected_tr=~/[1234]\=(.+)/; # it may be empty
		    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
		    $exons_tr=$#{$cds_lines{$tr}}+1;
		    $selected_tr="4=$tr" unless ($selected_tr eq "4=$ref_tr" || $selected_tr=~/[123]\=/ || $exons_sel>=$exons_tr);
		    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_ref};
		}
		elsif ($tr_co{$tr}{$C1_inc}){
		    ($temp_tr)=$selected_tr=~/[12345]\=(.+)/; # it may be empty
		    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
		    $exons_tr=$#{$cds_lines{$tr}}+1;
		    $t_C2=$array_tr_co{$tr}[$index_co{$tr}{$C1_inc}+1];
		    ($i,$f)=$t_C2=~/\:(.+?)\-(.+)/;
		    $OK="";
		    $OK=1 if $i>$Af && $str{$g} eq "+";
		    $OK=1 if $f<$Ai && $str{$g} eq "-";
		    $selected_tr="5=$tr" unless ($selected_tr eq "5=$ref_tr" || $selected_tr=~/[1234]\=/ || $exons_sel>=$exons_tr || !$OK);
		    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_inc};
		}
		elsif ($tr_co{$tr}{$C1_ref}){
		    ($temp_tr)=$selected_tr=~/[123456]\=(.+)/; # it may be empty
		    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
		    $exons_tr=$#{$cds_lines{$tr}}+1;
		    $t_C2=$array_tr_co{$tr}[$index_co{$tr}{$C1_ref}+1];
		    ($i,$f)=$t_C2=~/\:(.+?)\-(.+)/;
		    $OK="";
		    $OK=1 if $i>$Af && $str{$g} eq "+";
		    $OK=1 if $f<$Ai && $str{$g} eq "-";
		    $selected_tr="6=$tr" unless ($selected_tr eq "6=$ref_tr" || $selected_tr=~/[12345]\=/ || $exons_sel>=$exons_tr || !$OK);
		    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_ref};
		}
		elsif (($tr_coA{$tr}{$C1_incA} && $str{$g} eq "-") || ($tr_coB{$tr}{$C1_incA} && $str{$g} eq "+") ){ # C1do exists in transcript
		    ($temp_tr)=$selected_tr=~/[1234567]\=(.+)/; # it may be empty
		    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
		    $exons_tr=$#{$cds_lines{$tr}}+1;
		    $t_C2=$array_tr_co{$tr}[$index_coA{$tr}{$C1_incA}+1]; # gets a proper C2, not an acceptor
		    ($i,$f)=$t_C2=~/\:(.+?)\-(.+)/;
		    $OK="";
		    $OK=1 if $i>$Af && $str{$g} eq "+";
		    $OK=1 if $f<$Ai && $str{$g} eq "-";
		    $selected_tr="7=$tr" unless ($selected_tr eq "7=$ref_tr" || $selected_tr=~/[123456]\=/ || !$OK || $exons_sel>=$exons_tr);
		    $C1_accepted_index{$ev}{$tr}=$index_coA{$tr}{$C1_incA};
		}
		elsif (($tr_coA{$tr}{$C1_refA} && $str{$g} eq "-") || ($tr_coB{$tr}{$C1_refA} && $str{$g} eq "+")){ # C1do exists in transcript
		    ($temp_tr)=$selected_tr=~/[12345678]\=(.+)/; # it may be empty
		    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
		    $exons_tr=$#{$cds_lines{$tr}}+1;
		    $t_C2=$array_tr_co{$tr}[$index_coA{$tr}{$C1_refA}+1]; # gets a proper C2, not an acceptor
		    ($i,$f)=$t_C2=~/\:(.+?)\-(.+)/;
$OK="";
		    $OK=1 if $i>$Af && $str{$g} eq "+";
		    $OK=1 if $f<$Ai && $str{$g} eq "-";
#print "$ev\t$A_ref{$ev}\t$t_C2\tOK\n";
		    $selected_tr="8=$tr" unless ($selected_tr eq "8=$ref_tr" || $selected_tr=~/[1234567]\=/ || !$OK || $exons_sel>=$exons_tr);
		    $C1_accepted_index{$ev}{$tr}=$index_coA{$tr}{$C1_refA};
		}
		elsif ($tr_co{$tr}{$C2_inc}){
		    ($temp_tr)=$selected_tr=~/[123456789]\=(.+)/; # it may be empty
		    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
		    $exons_tr=$#{$cds_lines{$tr}}+1;
		    $t_C1=$array_tr_co{$tr}[$index_co{$tr}{$C2_inc}-1];
		    ($i,$f)=$t_C1=~/\:(.+?)\-(.+)/;
$OK="";
		    $OK=1 if $f<$Ai && $str{$g} eq "+";
		    $OK=1 if $i>$Af && $str{$g} eq "-";
		    $selected_tr="9=$tr" unless ($selected_tr eq "9=$ref_tr" || $selected_tr=~/[12345678]\=/ || !$OK || $exons_sel>=$exons_tr);
		    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C2_inc}-1;
		}
		elsif ($tr_co{$tr}{$C2_ref}){
		    ($temp_tr)=$selected_tr=~/[1234567890]\=(.+)/; # it may be empty
		    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
		    $exons_tr=$#{$cds_lines{$tr}}+1;
		    $t_C1=$array_tr_co{$tr}[$index_co{$tr}{$C2_ref}-1];
		    ($i,$f)=$t_C1=~/\:(.+?)\-(.+)/;
		    $OK="";
		    $OK=1 if $f<$Ai && $str{$g} eq "+";
		    $OK=1 if $i>$Af && $str{$g} eq "-";
		    $selected_tr="10=$tr" unless ($selected_tr eq "10=$ref_tr" || $selected_tr=~/[123456789]\=/ || !$OK || $exons_sel>=$exons_tr);
		    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C2_ref}-1;
		}
	    }
	    if ($selected_tr){
		$tally_non_annot_hit++;
		($rescue_type,$final_tr)=$selected_tr=~/(.+?)\=(.+)/;
		$tally_solutions{$rescue_type}++;

		    ### creating fake transcripts:
		$exN=0;
		$out_frame="";
		$started="";
		$finished="";
		$comment="";
		$CDS_seq="";
		$CDS_ini=$CDS_end=0;
		$STOP_detected="";
		$first_offset="";

		foreach $lineN (0..$#{$tr_lines{$final_tr}}){
		    $exN++;

############### exon
		    @t_line=split(/\t/,$tr_lines{$final_tr}[$lineN]);
		    $t_line[8]="gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
		    $full_line=join("\t",@t_line);
		    print O "$full_line\n";

############### CDS
		    if ($cds_lines{$final_tr}[$lineN] && !$out_frame){ # if !out_frame is either in-frame OR before the exon
			@t_line2=split(/\t/,$cds_lines{$final_tr}[$lineN]);
			    $t_line2[8]="gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; protein_id \"$final_tr"."fB$tally_non_annot_hit\"\;".
				" gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
			$full_line2=join("\t",@t_line2);
			print O "$full_line2\n";

			if (!$started){
			    $first_offset=$t_line2[7] if !$started;
			    $CDS_seq="N" if $first_offset==2;
			    $CDS_seq="NN" if $first_offset==1;
			    $CDS_seq="" if $first_offset==0; # just in case, reinitialized
			}
			$started=1; # means the CDS has already begun

			    #### prepare the translation:
			$le_CDS_ex=$cds_end{$final_tr}[$lineN]-$cds_ini{$final_tr}[$lineN]+1;
			$CDS_ini=$CDS_end+1;
			$CDS_end=$CDS_ini+$le_CDS_ex-1;
			$chr=$t_line2[0];
			$str=$t_line2[6];
			if ($str eq "+"){
			    $seq_bit=substr($seq{$chr},$cds_ini{$final_tr}[$lineN]-1,$le_CDS_ex);
			    $CDS_seq.=$seq_bit;
			}
			elsif ($str eq "-"){
			    $seq_bit=substr($seq{$chr},$cds_ini{$final_tr}[$lineN]-1,$le_CDS_ex);
			    $seq_bit=join("", reverse split (//, $seq_bit));
			    $seq_bit=~tr/ACGTacgt/TGCAtgca/;
			    $CDS_seq.=$seq_bit;    
			}
		    }
############### CDS when a novel exon was introduced and it is not 3n
		    elsif ($tr_lines{$final_tr}[$lineN] && $out_frame && !$STOP_detected){ 
			@t_line2=split(/\t/,$tr_lines{$final_tr}[$lineN]);
			$new_offset=0 if length($CDS_seq)%3==0;
			$new_offset=1 if length($CDS_seq)%3==2;
			$new_offset=2 if length($CDS_seq)%3==1;
			    
			print "$ev\t$t_line2[3]\t$new_offset\t$t_line2[7]\t$CDS_seq\n" if $ev eq "DreEX0000046";

			@t_line=split(/\t/,$tr_lines{$final_tr}[$lineN]); # uses the exon line again
			$Ai=$t_line[3]; # treated like a AS exon
			$Af=$t_line[4];

			    #### prepare the translation:
			$le_CDS_ex=$Af-$Ai+1;
			$CDS_ini=$CDS_end+1; # from upstream
			$CDS_end=$CDS_ini+$le_CDS_ex-1;
			$chr=$t_line[0]; # as it's feeding from exon line
			$str=$t_line[6]; # as it's feeding from exon line
			    
			if ($new_offset==0){
			    $extra="";
			}
			elsif ($new_offset==1){
			    ($extra)=$CDS_seq=~/.+(..)/;
    }
			elsif ($new_offset==2){
			    ($extra)=$CDS_seq=~/.+(.)/;
    }
    
			if ($str eq "+"){
			    $seq_bit=substr($seq{$chr},$Ai-1,$le_CDS_ex); # novel exon sequence
			    $seq_bit_test="$extra$seq_bit";
			    $loop=0;
			  LBL:while ($seq_bit_test=~s/(.{3})//){
			      $codon=$1;
			      $loop++; # 27/11/16 (stops in between EEJs)
			      if ($codon eq "TGA" || $codon eq "TAG" || $codon eq "TAA"){
				  $STOP_detected=1; # so it doesn't go to more CDS lines after
				  $stop_co_f=$Af-length($seq_bit_test);
				  $stop_co_i=$stop_co_f-2;
				  $co_Af_CDS=$stop_co_i-1;
				        
				  if ($loop>1){
				        $new_line2="$t_line[0]\t$t_line[1]\tCDS\t$Ai\t$co_Af_CDS\t$t_line[5]\t$t_line[6]\t$new_offset\t".
					          "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; protein_id \"$final_tr".
						  "fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
					print O "$new_line2\n";
					  $new_line3="$t_line[0]\t$t_line[1]\tstop_codon\t$stop_co_i\t$stop_co_f\t$t_line[5]\t$t_line[6]\t$t_line[7]\t".
					      "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
					print O "$new_line3\n";
				  }
				  elsif ($loop==1){
				      $bit=$t_line[3]+2-length($extra);
				        $new_line3="$t_line[0]\t$t_line[1]\tstop_codon\t$t_line[3]\t$bit\t$t_line[5]\t$t_line[6]\t$t_line[7]\t".
					    "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
				      print O "$new_line3\n";
				  }
				  last LBL;
			      }
			  }
			}
			elsif ($str eq "-"){
			    $seq_bit=substr($seq{$chr},$Ai-1,$le_CDS_ex);
			    $seq_bit=join("", reverse split (//, $seq_bit));
			    $seq_bit=~tr/ACGTacgt/TGCAtgca/;
			    $seq_bit_test="$extra$seq_bit";
			    $loop=0;
			  LBL:while ($seq_bit_test=~s/(.{3})//){
			      $codon=$1;
			      $loop++; # 27/11/16 (stops in between EEJs)
			      if ($codon eq "TGA" || $codon eq "TAG" || $codon eq "TAA"){
				  $STOP_detected=1; # so it doesn't go to more CDS lines after
				  $stop_co_i=$Ai+length($seq_bit_test); # as f has to be > i
				  $stop_co_f=$stop_co_i+2;
				  $co_Ai_CDS=$stop_co_f+1;

				  if ($loop>1){
				        $new_line2="$t_line[0]\t$t_line[1]\tCDS\t$co_Ai_CDS\t$Af\t$t_line[5]\t$t_line[6]\t$new_offset\t".
					          "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; protein_id \"$final_tr".
						  "fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
					print O "$new_line2\n";
					  $new_line3="$t_line[0]\t$t_line[1]\tstop_codon\t$stop_co_i\t$stop_co_f\t$t_line[5]\t$t_line[6]\t$t_line[7]\t".
					      "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
					print O "$new_line3\n";
				  }
				  elsif ($loop==1){
				      $bit=$t_line[4]-2+length($extra);
				        $new_line3="$t_line[0]\t$t_line[1]\tstop_codon\t$bit\t$t_line[4]\t$t_line[5]\t$t_line[6]\t$t_line[7]\t".
					    "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
				      print O "$new_line3\n";
				  }
				  last LBL;
			      }
			  } 
			}
			if  (!$STOP_detected){
			    $CDS_seq.=$seq_bit;    
			    $new_line2="$t_line[0]\t$t_line[1]\tCDS\t$Ai\t$Af\t$t_line[5]\t$t_line[6]\t$new_offset\t".
				    "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; protein_id \"$final_tr".
				    "fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
			    print O "$new_line2\n";
			    
#if ($length_A%3==0){
#    $comment="Alt_prot"; # all set
#}
#else {
#    $out_frame=1; # and moves the problem to the "next" CDS
#}
			}
		    }
		    
############### start codon
		    if ($start_lines{$final_tr}[$lineN] && !$out_frame){
			@t_line3=split(/\t/,$start_lines{$final_tr}[$lineN]);
			    $t_line3[8]="gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\;".
				" gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
			$full_line3=join("\t",@t_line3);
			print O "$full_line3\n";
		    }
############### stop codon
		    if ($stop_lines{$final_tr}[$lineN] && !$out_frame){
			@t_line4=split(/\t/,$stop_lines{$final_tr}[$lineN]);
			    $t_line4[8]="gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\;".
				" gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
			$full_line4=join("\t",@t_line4);
			print O "$full_line4\n";
			$finished=1; # means the CDS has already finished
		    }

############### Adding the exon
		    if ($lineN==$C1_accepted_index{$ev}{$final_tr}){
			$exN++;
			    #### exon
			($Ai,$Af)=$A_ref{$ev}=~/\:(.+?)\-(.+)/;
			$length_A=$Af-$Ai+1;
			    $new_line="$t_line[0]\t$t_line[1]\texon\t$Ai\t$Af\t$t_line[5]\t$t_line[6]\t$t_line[7]\t".
				"gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
			print O "$new_line\n";
			    #### evaluate effects:
			if (!$started){ # it's in the 5' UTR and no problems downstream
			    $comment="UTR_5";
			}
			elsif ($started && $finished){ # it is in the 3' UTR. It may coincide in BETWEEN exons
			    $comment="UTR_3";
			}
			elsif ($started && !$finished){ # the big fun
			    #### prepare the translation:
			    $le_CDS_ex=$Af-$Ai+1;
			    $CDS_ini=$CDS_end+1; # from upstream
			    $CDS_end=$CDS_ini+$le_CDS_ex-1;
			    $chr=$t_line[0]; # as it's feeding from exon line
			    $str=$t_line[6]; # as it's feeding from exon line
			    #### can take the offset of C2:
			    @t_line_C2=split(/\t/,$cds_lines{$final_tr}[$lineN+1]);
			    $offset_A=$t_line_C2[7];
			    if ($offset_A==0){
				$extra="";
			    }
			    elsif ($offset_A==1){
				($extra)=$CDS_seq=~/.+(..)/;
			    }
			    elsif ($offset_A==2){
				($extra)=$CDS_seq=~/.+(.)/;
			    }
			    
			    if ($str eq "+"){
				$seq_bit=substr($seq{$chr},$Ai-1,$le_CDS_ex); # novel exon sequence
				$seq_bit_test="$extra$seq_bit";
				$loop=0;

			      LBL:while ($seq_bit_test=~s/(.{3})//){
				  $codon=$1;
				  $loop++; # 27/11/16 (stops in between EEJs)
				        
				  if ($codon eq "TGA" || $codon eq "TAG" || $codon eq "TAA"){
				      $STOP_detected=1; # so it doesn't go to more CDS lines after
				      $stop_co_f=$Af-length($seq_bit_test);
				      $stop_co_i=$stop_co_f-2;
				      $co_Af_CDS=$stop_co_i-1;

				      if ($loop>1){
					        $new_line2="$t_line[0]\t$t_line[1]\tCDS\t$Ai\t$co_Af_CDS\t$t_line[5]\t$t_line[6]\t$offset_A\t".
						    "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
						print O "$new_line2\n";
						      $new_line3="$t_line[0]\t$t_line[1]\tstop_codon\t$stop_co_i\t$stop_co_f\t$t_line[5]\t$t_line[6]\t$t_line[7]\t".
							    "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; protein_id \"$final_tr".
							    "fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
						print O "$new_line3\n";
				      }
				      elsif ($loop==1){
					  $bit=$t_line[3]+2-length($extra);
					        $new_line3="$t_line[0]\t$t_line[1]\tstop_codon\t$t_line[3]\t$bit\t$t_line[5]\t$t_line[6]\t$t_line[7]\t".
						      "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; protein_id \"$final_tr".
						      "fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
					  print O "$new_line3\n";      
				      }

				      $comment="In_frame_STOP";
				      $out_frame=1;
				      last LBL;
				  }
			      }
			    }
			    elsif ($str eq "-"){
				$seq_bit=substr($seq{$chr},$Ai-1,$le_CDS_ex);
				$seq_bit=join("", reverse split (//, $seq_bit));
				$seq_bit=~tr/ACGTacgt/TGCAtgca/;
				$seq_bit_test="$extra$seq_bit";
				$loop=0;
			      LBL:while ($seq_bit_test=~s/(.{3})//){
      $codon=$1;
      $loop++; # 27/11/16 (stops in between EEJs)
      if ($codon eq "TGA" || $codon eq "TAG" || $codon eq "TAA"){
	  $STOP_detected=1; # so it doesn't go to more CDS lines after
	  $stop_co_i=$Ai+length($seq_bit_test); # as f has to be > i
	  $stop_co_f=$stop_co_i+2;
	  $co_Ai_CDS=$stop_co_f+1;

	  if ($loop>1){
	            $new_line2="$t_line[0]\t$t_line[1]\tCDS\t$co_Ai_CDS\t$Af\t$t_line[5]\t$t_line[6]\t$offset_A\t".
			  "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; protein_id \"$final_tr".
			  "fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
		    print O "$new_line2\n";
		          $new_line3="$t_line[0]\t$t_line[1]\tstop_codon\t$stop_co_i\t$stop_co_f\t$t_line[5]\t$t_line[6]\t$t_line[7]\t".
			      "gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
		    print O "$new_line3\n";
	  }
	  elsif ($loop==1){
	      $bit=$t_line[4]-2+length($extra);
	            $new_line3="$t_line[0]\t$t_line[1]\tstop_codon\t$bit\t$t_line[4]\t$t_line[5]\t$t_line[6]\t$t_line[7]\t".
			"gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
	      print O "$new_line3\n";
	  }
	  $comment="In_frame_STOP";
	  $out_frame=1;
	  last LBL;
      }
			      } 
			    }
			    if  (!$STOP_detected){
				$CDS_seq.=$seq_bit;
				    $new_line2="$t_line[0]\t$t_line[1]\tCDS\t$Ai\t$Af\t$t_line[5]\t$t_line[6]\t$offset_A\t".
					"gene_id \"$g\"\; transcript_id \"$final_tr"."fB$tally_non_annot_hit\"\; protein_id \"$final_tr".
					"fB$tally_non_annot_hit\"\; gene_name \"$gene_name{$g}\"\; exon_number \"$exN\"\;";
				print O "$new_line2\n";
				    
				if ($length_A%3==0){
				    $comment="Alt_prot"; # all set
				}
				else {
				    $out_frame=1; # and moves the problem to the "next" CDS
				    $comment="Frame_shifting";
				}
			    }
			}
			print LOG "$ev\t$g\t$exN\t$final_tr"."fB$tally_non_annot_hit\t$comment\t$rescue_type\n";
		    }
		}
	    }
	    elsif (!$selected_tr && $partial{$ev}){
		$tally_partial++;
	    }
	    else {
		$tally_non_annot_no_hit++;
#    print "$g\t$ev\t$A_ref{$ev}\n";
	    }
	}
    }

print "Fully Annotated\t$tally_annotated\n".
    "Rescued skipping\t$tally_non_annot_hit\n".
    "Partially Annotated\t$tally_partial\n".
    "Not rescued\t$tally_non_annot_no_hit\n";
print LOG2 "Fully Annotated\t$tally_annotated\n".
    "Rescued skipping\t$tally_non_annot_hit\n".
    "Partially Annotated\t$tally_partial\n".
    "Not rescued\t$tally_non_annot_no_hit\n";

    print "\nSolutions:\n";
    print LOG2 "\nSolutions:\n";
    foreach $type (sort {$a<=>$b} keys %tally_solutions){
	print "Type $type\t$tally_solutions{$type}\n";
	print LOG2 "Type $type\t$tally_solutions{$type}\n";
    }

##################################
###Finish IntroduceEXON_to_GTF####

######Get exint file with added exons######
    my $tmp_exint_output = "$exons_db_folder/$species/tmp_$species.exint";
    my $fakeGTF="$exons_db_folder/$species/FakeTranscripts-$species-vB.gtf";
    my (%ntseq2, %rseq2);
#1) Read FakeGTF
    open(FAKEGTF, $fakeGTF) || die "cannot open $fakeGTF";
    while (<FAKEGTF>){
	chomp($_); 
	@line=split(/\t/,$_);
	if ($seq{$line[0]}){
	    if ($line[2] eq "CDS"){
		@l1=split(/\"/,$line[8]);
		$gid=$l1[1]; ##change here
		    #print "GID: $gid\n";
		if ($line[8]=~/protein_id/){
		    $tmp="protein_id";
		} 
		else { $tmp="transcript_id"; }
		@l2=split(/$tmp/,$line[8]);
		@l3=split(/\"/,$l2[1]);
		$prot=$l3[1];
		if (!$pid{$prot}){
		    $m=0;
		    if ($line[7] ne "."){  $phase=$line[7]; $fex{$prot}=$phase; $m=$phase; } ##
		    else { 
			$phase=0; 
			$fex{$prot}=$phase;
		    }
		    $nuc=0;
		    $size=$line[4]-$line[3]+1;
		    $tmpseq=substr($seq{$line[0]},($line[3]-1),($size));
		    $size=$size-$m;
		    if ($line[6] eq "-"){  $tmpseq=~tr/ATCG/TAGC/; $tmpseq=reverse($tmpseq); }
		    $ntseq2{$prot}=$tmpseq;
		    $rseq2{$prot}=$tmpseq;
		    $res=$size%3;
		    $nuc=$nuc+$size;
		    $aa=int(1+$nuc/3); 
		    $pid{$prot}=1;
		    $header{$prot}="$prot|$gid ";
		}
		else {
		    $header{$prot}.=" $aa.";
		    if ($line[7] ne "."){  $phase=$line[7]; }
		    else {  
			if ($res==0){ $phase=0; } elsif ($res==1 ){ $phase=2; } elsif ($res==2) { $phase=1; } 
		    }
		    $size=($line[4]-$line[3])+1;
		    $tmpseq=substr($seq{$line[0]},($line[3]-1),($size));
		    if ($line[6] eq "-"){  $tmpseq=~tr/ATCG/TAGC/; $tmpseq=reverse($tmpseq); }
		    $ntseq2{$prot}.=$tmpseq;
		    $rseq2{$prot}.="|".$tmpseq;
		    $nuc=$nuc+$size;
		    $size=$size+$res;
		    $res=$size%3;
		    $aa=int(1+$nuc/3); ##adding the residues that were left
		    $header{$prot}.="$phase";
		}
	    }
	}
    }
    close (GTF);
    
#2) Get exint file for fake transcripts
    my @keys2=keys(%ntseq2);
    verbPrint("Generating exint files with added exons for $species\n");
    open (EXINT_OUT, ">$tmp_exint_output") || die "Cannot open $exint_output (loop 1)\n";
    foreach $el(@keys2){ 
	$nseq=$ntseq2{$el};
	$name=$header{$el};
	$p=$fex{$el}; ##checking phase of the first exon
	$protein="";
	$c=0;
	for($l=$p; $l<length($nseq); $l=$l+3){
	    $triplet=substr($nseq,$l,3);
	    if ($c==0){
		if ($codons{$triplet}){
		    if ($l>=(length($nseq)-3) && $codons{$triplet} eq "stop"){
			$c=1;
		    }
		    elsif ($codons{$triplet} eq "stop"){ $c=1; } #$protein.=$codons{$triplet}; } #$c=1; }
		    else {  $protein.=$codons{$triplet};   }
		}
		else {
		    $protein.="X";
		}
	    }
	}
	print EXINT_OUT ">$name\n$protein\n";
    }
    close (EXINT_OUT);

###Joining 1. GTF files (Annot+Fake)
    my $tmp="$exons_db_folder/$species/tmp.gtf";
    `cat $gtf_file $fakeGTF > $tmp`;
    my $ntmp=$exons_db_folder."/".$species."/".$species."_annot_fake.gtf";
    `mv $tmp $ntmp`;

###Joining 2. Exint files
    my $tmp="$exons_db_folder/$species/tmp.exint";
    `cat $exint_output $tmp_exint_output > $tmp`;
    `mv $tmp $exint_output`;


}##if added exons ne "NA" --> if additional exons are provided by the user 

### loops for each species (loop 2: get_trs_gtf.pl)
#foreach my $species (@SPECIES){
my (@line, @coords1, @coords2, @gene, @l1, @l2, @l3);
my ($l, $grep, $gid, $int, $r, $prot, $tmp, $trid, $n);
my (%tr1, %tr2, %tpid);
my %prot;
my (%coords1, %coords2);
my (%strand1, %strand2);
my %intron;
my (%chr1, %chr2);
if ($gtf_file =~ /.gz$/) {
    open(INFILEONE, "gunzip -c $gtf_file |") || die "cannot open pipe to $gtf_file";
}
else {
    open(INFILEONE, $gtf_file) || die "t open $gtf_file";
}    
while (<INFILEONE>){ 
    if($_){ 
	chomp($_);
	@line=split(/\t/,$_);
	if ($line[2] eq "exon"){
	    @l1=split(/\"/,$line[8]);
	    $gid=$l1[1];
	    $tmp="transcript_id";
	    @l2=split(/$tmp/,$line[8]);
	    @l3=split(/\"/,$l2[1]);
	    $trid=$l3[1];
	    $strand1{$trid."|".$gid}=$line[6];
	    $chr1{$trid."|".$gid}=$line[0];
	    if (!$tr1{$trid."|".$gid}){
		$tr1{$trid."|".$gid}=$line[3].",".$line[4];
	    }
	    else {
		$tr1{$trid."|".$gid}.="-".$line[3].",".$line[4];
	    }
	}
	elsif ($line[2] eq "CDS"){
	    @l1=split(/\"/,$line[8]);
	    if ($line[8]=~/protein_id/) { $tmp="protein_id"; }
	    else { $tmp="transcript_id"; }
	    @l2=split(/$tmp/,$line[8]);
	    @l3=split(/\"/,$l2[1]);
	    $prot=$l3[1];
	    $tmp="transcript_id";
	    @l2=split(/$tmp/,$line[8]);
	    @l3=split(/\"/,$l2[1]);
	    $trid=$l3[1];
	    $tpid{$trid."\t".$prot}=1;
	    $strand2{$prot."|".$gid}=$line[6];
	    $chr2{$prot."|".$gid}=$line[0];
	    if (!$tr2{$prot."|".$gid}){ ##modifying when the first exon of the CDS is not in phase 0
		if ($line[7] == 1){
		    if ($line[6] eq "+"){ 
			$n=$line[3]+1;
			$tr2{$prot."|".$gid}=$n.",".$line[4];
		    } else { 
			$n=$line[4]-1;
			$tr2{$prot."|".$gid}=$line[3].",".$n;
		    }
		}
		elsif ($line[7] == 2) {
		    if ($line[6] eq "+"){ 
			$n=$line[3]+2;
			$tr2{$prot."|".$gid}=$n.",".$line[4];
		    } else { 
			$n=$line[4]-2;
			$tr2{$prot."|".$gid}=$line[3].",".$n;
		    }
		}
		else {
		    $tr2{$prot."|".$gid}=$line[3].",".$line[4];
		}
	    }
	    else {
		$tr2{$prot."|".$gid}.="-".$line[3].",".$line[4];
	    }
	}
    }
}
close (INFILEONE);
### Generates the outputs for Loop 2
# -o1 Hsa_tr_coords.txt -o2 Hsa_tr_coords_CDS.txt -o3 Hsa_trid_protid.txt -o4 Hsa_annot_exons_prot_ids.txt
my $out1 = "$exons_db_folder/$species/$species"."_tr_coords.txt";
my $out2 = "$exons_db_folder/$species/$species"."_tr_coords_CDS.txt";
my $out3 = "$exons_db_folder/$species/$species"."_trid_protid.txt";
my $out4 = "$exons_db_folder/$species/$species"."_annot_exons_prot_ids.txt";
my @keys=keys(%tr1);
my $element;
verbPrint("Generating $species"."_tr_coords.txt\n");
open (OUTFILE1, ">$out1") || die "Cannot open $out1 (loop 2)\n";
foreach $element (@keys){
    if ($strand1{$element} eq "+"){
	print OUTFILE1 "$element\t$chr1{$element}\t$strand1{$element}\t$tr1{$element}\n";
    }
    else {
	@coords1=split(/\-/,$tr1{$element});
	print OUTFILE1 "$element\t$chr1{$element}\t$strand1{$element}\t";
	for ($l=scalar(@coords1)-1; $l>=1; $l--){
	    print OUTFILE1  "$coords1[$l]-";
	}
	print OUTFILE1 "$coords1[0]\n";
    }
}
close OUTFILE1;
    
my @keys2=keys(%tr2);
my $element2;
my %print;
verbPrint("Generating $species"."_tr_coords_CDS.txt and $species"."_trid_protid.txt\n");
open (OUTFILE2, ">$out2") || die "Cannot open $out2 (loop 2)\n";
open (OUTFILE4, ">$out4")  || die "Cannot open $out4 (loop 2)\n";
foreach $element2 (@keys2){
    if ($strand2{$element2} eq "+"){
	print OUTFILE2 "$element2\t$chr2{$element2}\t$strand2{$element2}\t$tr2{$element2}\n";
	if (!$print{$element2}){
	    print OUTFILE4 "$element2\n";
	    $print{$element2}=1;
	}
    }
    else {
	@coords2=split(/\-/,$tr2{$element2});
	print OUTFILE2 "$element2\t$chr2{$element2}\t$strand2{$element2}\t";
	if (!$print{$element2}){
	    print OUTFILE4 "$element2\n";
	    $print{$element2}=1;
	}
	for ($l=scalar(@coords2)-1; $l>=1; $l--){
	    print OUTFILE2  "$coords2[$l]-";
	}
	print OUTFILE2 "$coords2[0]\n";
    }
}
close OUTFILE2;
close OUTFILE4;

my @keys3=keys(%tpid);
@keys=sort(@keys3);
my $tp;
verbPrint("Generating $species"."_annot_exons_prot_ids.txt\n");
open (OUTFILE3, ">$out3") || die "Cannot open $out3 (loop 2)\n";
foreach $tp (@keys3){
    print OUTFILE3 "$tp\n";
}
close OUTFILE3;
### Only if VASTDB files provided
if ($add_exons ne "NA"){
    ### loops for each species (Loop 3: get_prot_isof_exon_2.pl)
    #    foreach my $species (@SPECIES){  
        # Declaration of variables
    my (@line);
    my %pid;
    my %lprot;
    my ($ex, $l, $intron, $ex2);
    my (@t1,@t2, @c);
    my %psize;
    my $pr;   
    # Skips the step for the species if NA was provided
    next if $add_exons eq "NA";    
    # Declaration of input files
    my $i1 = "$exons_db_folder/$species/$species"."_prot_sizes.txt";
    my $i2 = "$exons_db_folder/$species/$species"."_tr_coords_CDS.txt";
    my $i3 = $add_exons;    
    # Declaration of output files
    my $out = "$exons_db_folder/$species/Ref_protein_exons_".$species."_1.txt";
    my $out2 = "$exons_db_folder/$species/missing_exons.txt";
    open (OUT,">$out") || die "Cannot open $out (loop 3)\n";
    open (MISS,">$out2") || die "Cannot open $out2 (loop 3)\n";    
    verbPrint("Generating Ref_protein_exons_".$species."_1.txt and missing_exons.txt from $add_exons file\n");
    open (INONE, $i1) || die "Cannot open $i1 (loop 3)\n";
    # Format: BL21561_evm0BL21561476
    while (<INONE>){ 
	chomp($_);
	@line=split(/\t/,$_);
	$pr=$line[0]."|".$line[1];
	$psize{$pr}=$line[2];
    }
    close (INONE);    
    my ($tmp, $k, $r, $s);
    my %coords;
    open (INTWO, $i2) || die "Cannot open $i2 (loop 3)\n"; 
    #Format: BL20899_cuf1|BL20899Sc0000317-255783,255866-258240,258374-259202,259294
    while (<INTWO>){ 
	chomp($_);
	@line=split(/\t/,$_);
	$s=$line[3]; $s=~s/\,/\|/g; $s=~s/\-/\,/g; $s=~s/\|/\-/g;
	$coords{$line[0]}=$line[1].":".$line[3].":".$line[2];    
	if ($line[3]=~/\-/){
	    @c=split(/\-/,$_);
	    for ($l=1; $l<scalar(@c)-1; $l++){
		my $r=$l-1; my $k=$l+1;
		@t1=split(/\,/,$c[$r]); @t2=split(/\,/,$c[$k]);
		$tmp=$c[$l];
		$tmp=~s/\,/\-/;
		$ex2=$line[1].":".$tmp;
		if ($line[2] eq "+"){ 
		        # WARNING :This line emited a nasty warning, but it's probably OK (the + value is ignored)
		    $t1[1]= "NA" if !defined $t1[1];
		    $ex=$line[1].":".$t1[1].",".$tmp.",".$t2[0].":".$line[2];
		}
		else {
		        # WARNING :This line emited a nasty warning, but it's probably OK (the + value is ignored)
		    $t1[1]= "NA" if !defined $t1[1];
		    $ex=$line[1].":".$t2[0].",".$tmp.",".$t1[1].":".$line[2]; 
		}

		if (!$pid{$ex}){
		    $pid{$ex}=$line[0];
		    $pid{$ex2}=$line[0];
		    $lprot{$ex}=$psize{$line[0]};
		    $lprot{$ex2}=$psize{$line[0]};
		}
		else {  
		    if ($psize{$line[0]} > $lprot{$ex}){
			$pid{$ex}=$line[0];
			$pid{$ex2}=$line[0];
			$lprot{$ex}=$psize{$line[0]};
			$lprot{$ex2}=$psize{$line[0]};
		    }
		}
	    }
	    @c=split(/\,/,$_);
	    for ($l=1; $l<scalar(@c)-1; $l++){
		@t1=split(/\-/,$c[$l]);
		if ($line[2] eq "+"){
		    $intron=$line[1].":".$t1[0].",".$t1[1].":".$line[2];
		}
		else {  
		    $intron=$line[1].":".$t1[1].",".$t1[0].":".$line[2]; 
		}
		if (!$pid{$intron}){
		    $pid{$intron}=$line[0];
		    $lprot{$intron}=$psize{$line[0]};
		}
		else {  
		    if ($psize{$line[0]} > $lprot{$intron}){
			$pid{$intron}=$line[0];
			$lprot{$intron}=$psize{$line[0]};
		    }
		}
	    }
	}
    }
    close (INTWO);    
    my ($id,$id4, $trcoords);
    my @ex;
    open (INTHREE, $i3) || die "Cannot open $i3 (loop 3)\n"; 
# Format: FAM13ABlaEX0015343Sc0000095:756974-75707299Sc0000095:757428,756974-757072,756243S
#         Sc0000095:757428,756974-757072,756243:-99
#         Sc0000095:757428-757560Sc0000095:756974-757072Sc0000095:756125-756243
#Format: GeneID/GeneNameExonIDExon_coordsExon_coords:C1,A,C2:strand
#BL005147BlaEX0015343Sc0000095:756974-75707299 Sc0000095:757428,756974-757072,756243:- 
#ENSG00000029534HsaEX0004110chr8:41758036-41758137chr8:41797512-41797622chr8:41733971-41734069chr8:41797512,41758036-41758137,41734069:-
#From VastDB: Columns 1,2,3 & 7.
    while (<INTHREE>){ 
	chomp($_);
	@line=split(/\t/,$_);
	if ($pid{$line[3]}){
	    $trcoords=$coords{$pid{$line[5]}};
	    print OUT "$line[0]\t$line[1]\t$line[2]\t$pid{$line[5]}\tAnnotated_exon_C1,A,C2\t$trcoords\n";
	}
	elsif ($pid{$line[2]}) {
	    $trcoords=$coords{$pid{$line[2]}};
	    print OUT "$line[0]\t$line[1]\t$line[2]\t$pid{$line[5]}\tAnnotated_exon_A\t$trcoords\n";
	}
	else {
	    @t1=split(/\,/,$line[5]);
	    $id=$t1[0].",".$t1[2]; ##searching then only the intron, the exon might be not annotated
	    if ($pid{$id}){
		$trcoords=$coords{$pid{$id}};
		print OUT "$line[0]\t$line[1]\t$line[2]\t$pid{$id}\tNot_annotated_exon|intron:$id\t$trcoords\n"; 
	    }
	    else { print MISS "$_\n"; }
	}
    }
    close (INTHREE);
    close OUT;
    close MISS;
### loops for each species (Loop 4: get_prot_isof_exon_3.pl) only if VASTDB refs
#   -sp Hsa 
#   -i1 Hsa_prot_sizes.txt -i2 Hsa_trid_protid.txt -i3 Hsa_tr_coords.txt -i4 Ref_protein_exons_Hsa_1.txt -i5 missing_exons.txt 
#   -o1 Ref_protein_exons_Hsa_2.txt -o2 missing_exons_2.txt -o3 Final_exons_ref_proteins_Hsa.txt -o4 Hsa_vast_exons_prot_ids.txt
#    foreach my $species (@SPECIES){
    # Declaration of variables
    my (@line);
    my %pid;
    my %lprot;
    my ($ex, $l, $intron);
    my (@t1,@t2, @c);
    my %gene;
    my %psize;
    my $pr;
# Skips the step for the species if NA was provided
    next if $add_exons eq "NA";
    # Declaration of input files
    my $i1 = "$exons_db_folder/$species/$species"."_prot_sizes.txt";
    my $i2 = "$exons_db_folder/$species/$species"."_trid_protid.txt";
    my $i3 = "$exons_db_folder/$species/$species"."_tr_coords.txt";
    my $i4 = "$exons_db_folder/$species/Ref_protein_exons_$species"."_1.txt";
    my $i5 = "$exons_db_folder/$species/missing_exons.txt";
    die "Cannot find $i4 (loop4)\n" unless (-e $i4); # sanity check
    
    # Declaration of output
    #   -o1 Ref_protein_exons_Hsa_2.txt -o2 missing_exons_2.txt -o3 Final_exons_ref_proteins_Hsa.txt -o4 Hsa_vast_exons_prot_ids.txt
    
    my $out = "$exons_db_folder/$species/Ref_protein_exons_$species"."_2.txt";
    my $out2 = "$exons_db_folder/$species/missing_exons_2.txt";
    my $out3 = "$exons_db_folder/$species/Final_exons_ref_proteins_$species.txt";
    my $out4 = "$exons_db_folder/$species/$species"."_add_exons_prot_ids.txt"; ##add_exons file prot ids
    my $mex = "$exons_db_folder/$species/$species"."_not_annotated_exons.txt";
    
    open (OUT,">$out") || die "Cannot open $out (loop4)\n";
    open (MISS,">$out2") || die "Cannot open $out2 (loop4)\n";
    open (OUTF,">$out4") || die "Cannot open $out4 (loop 4)\n"; # out3 is done by a cat
    open (MEX, ">$mex") || die "Cannot open $mex (loop 4)\n";
    
    ### Parsing first input
    open (INONE, $i1) || die "Cannot open $i1 (loop 4)\n";
    # Format: BL21561_evm0BL21561476
    while (<INONE>){ 
	chomp($_);
	@line=split(/\t/,$_);
	$pr=$line[0];
	$psize{$pr}=$line[2];
	$gene{$pr}=$line[1];
    }
    close (INONE);
    
    ### Parsing second input
    open (INTWO, $i2) || die "Cannot open $i2 (loop 4)\n";
    my %trid;
    while (<INTWO>){ 
	chomp($_);
	@line=split(/\t/,$_);
	$trid{$line[0]}=$line[1];
    }
    close (INTWO);
    
    ### Parsing third input
    my ($tmp, $k, $r, $s);
    my @n;
    my %coords;
    open (INTHR, $i3) || die "Cannot open $i3 (loop 4)\n";
    # Format: BL20899_cuf1|BL20899Sc0000317-255783,255866-258240,258374-259202,259294
    while (<INTHR>){ 
	chomp($_);
	@line=split(/\t/,$_);
	$s=$line[3]; $s=~s/\,/\|/g; $s=~s/\-/\,/g; $s=~s/\|/\-/g;
	@n=split(/\|/,$line[0]);
	    ### This gave WARNINGS: basically, nc transcripts do not have prot?
#    $trid{$n[0]}="NA" if (!$trid{$n[0]}); # this solution adds lines in the output for nc exons
	$coords{$trid{$n[0]}}=$line[1].":".$line[3].":".$line[2] if (defined $trid{$n[0]});
	if ($line[3]=~/\-/){
	    @c=split(/\-/,$_);
	    for ($l=1; $l<scalar(@c)-1; $l++){
		my $r=$l-1; my $k=$l+1;
		@t1=split(/\,/,$c[$r]); @t2=split(/\,/,$c[$k]);
		$tmp=$c[$l];
		$tmp=~s/\,/\-/;
		if ($line[2] eq "+"){
		        # WARNINGS: same issues as above
		    $t1[1]= "NA" if (!defined $t1[1]);
		    $ex=$line[1].":".$t1[1].",".$tmp.",".$t2[0].":".$line[2];
		}
		else {  
		        # WARNINGS: same issues as above
		    $t1[1]= "NA" if (!defined $t1[1]);
		    $ex=$line[1].":".$t2[0].",".$tmp.",".$t1[1].":".$line[2]; 
		}
		if (!$pid{$ex}){
		    $pid{$ex}=$trid{$n[0]};
		        # This gave warnings (likely for nc transcripts)
		    $lprot{$ex}=$psize{$trid{$n[0]}} if (defined $trid{$n[0]});
		}
		else {  
		    if (defined $trid{$n[0]}){
			if ($psize{$trid{$n[0]}} > $lprot{$ex}){
			    $pid{$ex}=$trid{$n[0]};
			    $lprot{$ex}=$psize{$trid{$n[0]}};
			}
		    }
		}
	    }
	    @c=split(/\,/,$_);
	    for ($l=1; $l<scalar(@c)-1; $l++){
		@t1=split(/\-/,$c[$l]);
		if ($line[2] eq "+"){
		        # WARNINGS: same issues as above
		    $t1[1]= "NA" if (!defined $t1[1]);
		    $intron=$line[1].":".$t1[0].",".$t1[1].":".$line[2];
		}
		else { 
		        # WARNINGS: same issues as above
		    $t1[1]= "NA" if (!defined $t1[1]);
		    $intron=$line[1].":".$t1[1].",".$t1[0].":".$line[2]; 
		}
		if (!$pid{$intron}){
		    $pid{$intron}=$trid{$n[0]};
		        # This gave warnings (likely for nc transcripts)
		    $lprot{$intron}=$psize{$trid{$n[0]}} if (defined $trid{$n[0]});
		}
		else {  
		    if (defined $trid{$n[0]}){ # This gave warnings (likely for nc transcripts) 
			if ($psize{$trid{$n[0]}} > $lprot{$intron}){
			    $pid{$intron}=$trid{$n[0]};
			    $lprot{$intron}=$psize{$trid{$n[0]}} if (defined $trid{$n[0]});
			}
		    }
		}
	    }
	}
    }
    close (INTHR);
    
### Parsing fifth input
    verbPrint("Generating Ref_protein_exons_$species"."_2.txt\n");
    verbPrint("Generating missing_exons_2.txt from $add_exons file");
    my ($id,$g, $trcoords);
    my @ex;
    open (INFIVE, $i5) || die "Cannot open $i5 (loop 4)\n";
# Format: FAM13ABlaEX0015343Sc0000095:756974-75707299Sc0000095:757428,756974-757072,756243S
#         Sc0000095:757428,756974-757072,756243:-99
#         Sc0000095:757428-757560Sc0000095:756974-757072Sc0000095:756125-756243
##BL005147BlaEX0015343Sc0000095:756974-75707299 Sc0000095:757428,756974-757072,756243:-
    while (<INFIVE>){ 
	chomp($_);
	@line=split(/\t/,$_);
	if ($pid{$line[3]}){
	    $g=$gene{$pid{$line[3]}};
	    $trcoords=$coords{$pid{$line[3]}};
	    print OUT "$line[0]\t$line[1]\t$line[2]\t$pid{$line[5]}|$g\tAnnotated_exon_C1,A,C2*\t$trcoords\n";
	}
	else {
	    print MISS "$_\n";
	}
    }
    close (INFIVE);
    close OUT;
    close MISS;

        ##Getting final file of exons
    verbPrint("Generating Final_exons_ref_proteins_$species.txt");
    system ("cat $out $i4 > $out3"); # uses forth input here
    
        ##Getting file with vast exons and protein ids
    verbPrint("Generating $species"."_add_exons_prot_ids.txt");
    verbPrint("Generating $species"."_not_annotated_exons.txt");
    
    my %print;
    open (IN, $out3) || die "Cannot open $out3 (loop 4)\n";
    while (<IN>){ 
	chomp($_);
	@line=split(/\t/,$_);
	if (!$print{$line[3]}){
	    print OUTF "$line[3]\n";
	    $print{$line[3]}=1;
	}
	if ($line[4]=~/Not_annotated_exon/){
	    print MEX "$line[0]\t$line[1]\t$line[2]\t$line[3]\n";
	}
    }
    close (MEX);
    close (OUTF);
    close IN;
}


### loops for each species (Loop 5: get_aa_pos_exon.pl)
#   -i Hsa_annot_exons_prot_ids.#   -out Hsa_protein_ids_exons_pos.txt
#foreach my $species (@SPECIES){  
    #Declaration of variables
my (%pid, %protein);
my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m, $ex);
my (%header, %ntseq);
#my $b;
my (@l, @l1, @line, @l2, @l3, @s);
    
    # Declaration of input files
my $p = "$exons_db_folder/$species/$species"."_annot_exons_prot_ids.txt";
my $v = "$exons_db_folder/$species/$species"."_add_exons_prot_ids.txt" if (-e "$exons_db_folder/$species/$species"."_add_exons_prot_ids.txt");
    #my $gtf_file = "$gtf_folder/$species"."_annot.gtf";
    
    # Declaration of output files
my $outfile = "$exons_db_folder/$species/$species"."_protein_ids_exons_pos.txt";
open (OUT, ">$outfile") || die "Cannot open $outfile (loop 5)\n";
    
verbPrint("Generating $species"."_protein_ids_exons_pos.txt\n");
    
    ### open main input file
open (PROTS, "$p") || die "Cannot open $p (loop 5)\n";
while (<PROTS>){
    chomp($_);
    $protein{$_}=1;
}
close (PROTS);
### if VASTDB refs were provided
if ($v){
    open (VAST, "$v");
    while (<VAST>){
	chomp($_);
	$protein{$_}=1;
    }
    close (VAST);
}
### Parsing GTF file
my %rseq;
my ($res, $size, $new, $id, $e);
my %fex; ##saving the phase of the first exon
if ($gtf_file =~ /.gz$/) {
    open(GTF, "gunzip -c $gtf_file |") || die "t open pipe to $gtf_file";
}
else {
    open(GTF, $gtf_file) || die "t open $gtf_file";
}
while (<GTF>){
    chomp($_); 
    @line=split(/\t/,$_);
    if ($line[2] eq "CDS"){
	@l1=split(/\"/,$line[8]);
	$gid=$l1[1]; ##change here
	if ($line[8]=~/protein_id/){ $tmp="protein_id"; }
	else { $tmp="transcript_id";  }
	@l2=split(/$tmp/,$line[8]);
	@l3=split(/\"/,$l2[1]);
	$prot=$l3[1];
	$id=$prot."|".$gid;
	$tmp="exon_number";
	@l2=split(/$tmp/,$line[8]);
	@l3=split(/\"/,$l2[1]);
	$ex=$l3[1];
	if ($protein{$id}){
	    if (!$pid{$prot}){
		$e=1;
		$pid{$prot}=1;
		$m=0;
		if ($line[7] ne "."){  $phase=$line[7]; $fex{$prot}=$phase; $m=$phase; } ##
		else { 
		    $phase=0; 
		    $fex{$prot}=$phase;
		}
		$nuc=0;
		$size=$line[4]-$line[3]+1;
		$size=$size-$m;
		$res=$size%3;
		$nuc=$nuc+$size;
		$aa=int(1+$nuc/3); 
		print OUT "$prot|$gid\texon_$e\t1-$aa\t$line[0]\t$line[3]-$line[4]\t$line[6]\n"; ###exon numbering according to CDS not to transcript, otherwise $ex variable shoud be use
	    }
	    else {
		$e++;
		$new=$aa+1;
		if ($line[7] ne "."){  $phase=$line[7]; }
		else {  
		    if ($res==0){ $phase=0; } elsif ($res==1 ){ $phase=2; } elsif ($res==2) { $phase=1; } 
		}
		$size=($line[4]-$line[3])+1;
		$nuc=$nuc+$size;
		$size=$size+$res;
		$res=$size%3;
		$aa=int(1+$nuc/3);
		print OUT "$prot|$gid\texon_$e\t$new-$aa\t$line[0]\t$line[3]-$line[4]\t$line[6]\n"; 
	    }
	}
    }
}
close (GTF);
close OUT;
#}
### loops for each species (Loop 6: get_genomic_coords_intron.pl
#   -p Hsa_annot_exons_prot_ids.txt -v Hsa_vast_exons_prot_ids.txt -t Hsa_trid_protid.txt 
#   -c1 Hsa_tr_coords.txt -c2 Hsa_tr_coords_CDS.txt 
#   -o1 Hsa_protein_ids_intron_pos.txt -o2 Hsa_protein_ids_intron_pos_CDS.txt
#foreach my $species (@SPECIES){  
    # Declaration of variables
my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
my (%pid, %tid, %ntseq);
#my $b=0;
my (@l, @l1, @line, @l2, @l3, @coords);
my ($l, $size, $tmpseq, $c);
# Declaration of input files
my $p = "$exons_db_folder/$species/$species"."_annot_exons_prot_ids.txt";
my $v = "$exons_db_folder/$species/$species"."_add_exons_prot_ids.txt" if (-e "$exons_db_folder/$species/$species"."_add_exons_prot_ids.txt");
my $t = "$exons_db_folder/$species/$species"."_trid_protid.txt";
my $i1 = "$exons_db_folder/$species/$species"."_tr_coords.txt";
my $i2 = "$exons_db_folder/$species/$species"."_tr_coords_CDS.txt";
# Declaration of output files
my $out1 = "$exons_db_folder/$species/$species"."_protein_ids_intron_pos.txt";
my $out2 = "$exons_db_folder/$species/$species"."_protein_ids_intron_pos_CDS.txt";
open (ONE, ">$out1") || die "Cannot create $out1 (loop 6)\n";
open (TWO, ">$out2") || die "Cannot create $out2 (loop 6)\n";
# parses input 1
open (PROTS,"$p") || die "Cannot open $p (loop 6)\n";
while (<PROTS>){
    chomp($_); 
    $pid{$_}=1;
}
close (PROTS);
# parses input 2
if ($v){
    open (VAST,"$v") || die "Cannot open $v (loop 6)\n";
    while (<VAST>){
	chomp($_); 
	$pid{$_}=1;
    }
    close (VAST);
}
# parses input 3
open (TRS,"$t") || die "Cannot open $t (loop 6)\n";
# Format: ENSBTAT00000000005ENSBTAP00000000005
while (<TRS>){
    chomp($_); 
    @line=split(/\t/,$_);
    $tid{$line[0]}=$line[1];
}
close (TRS);
verbPrint("Generating $species"."_protein_ids_intron_pos.txt\n");
# parses input 4
open (FILE,"$i1") || die "Cannot open $i1 (loop 6)\n";
while (<FILE>){
    chomp($_); 
    @line=split(/\t/,$_);
    @l=split(/\|/,$line[0]);
    if ($tid{$l[0]}){ 
	$tmp=$tid{$l[0]}."|".$l[1];
	if ($pid{$tmp}){
	    if ($line[3]=~/\-/){
		@coords=split(/\,/,$line[3]);
		if ($line[2] eq "+"){
		    $c=0;
		    for ($l=1; $l<scalar(@coords)-1; $l++){
			$c++;
			print ONE "$tmp\tintron\t$c\t$line[1]\t$coords[$l]\t$line[2]\n";
		    }
		}
		elsif ($line[2] eq "-"){
		    $c=0;
		    for ($l=scalar(@coords)-2; $l>=1; $l--){
			$c++;
			print ONE "$tmp\tintron\t$c\t$line[1]\t$coords[$l]\t$line[2]\n";
		    }
		}
	    }
	    else {
		print ONE "$tmp\tintronless\tNA\t$line[1]\tNA\t$line[2]\n";
	    }
	}
    }
}
close (FILE);
close (ONE);
verbPrint ("Generating $species"."_protein_ids_intron_pos.txt\n");
# parsing input 5
open (FILE2,"$i2") || die "Cannot open $i2 (loop 6)\n";
while (<FILE2>){
    chomp($_); 
    @line=split(/\t/,$_);
    $tmp=$line[0];
    if ($pid{$tmp}){
	if ($line[3]=~/\-/){
	    @coords=split(/\,/,$line[3]);
	    if ($line[2] eq "+"){
		$c=0;
		for ($l=1; $l<scalar(@coords)-1; $l++){
		    $c++;
		    print TWO "$tmp\t$line[1]\t$line[2]\tIntron_$c\t$coords[$l]\n";
		}
	    }
	    elsif ($line[2] eq "-"){
		$c=0;
		for ($l=scalar(@coords)-2; $l>=1; $l--){
		    $c++;
		    print TWO "$tmp\t$line[1]\t$line[2]\tIntron_$c\t$coords[$l]\n";
		}
	    }
	}
	else {
	    print TWO "$tmp\t$line[1]\t$line[2]\tIntronless\n";
	}
    }
}
close (FILE2);
close (TWO);
#}
### loops for each species (Loop 7: get_ex_size_phase.pl) => try to intergrate in exint?
# -GTF Hsa_annot.gtf -out Hsa_annot_exon_sizes.txt
#foreach my $species (@SPECIES){  
    # Declaration of variables (same as loop 1)
my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
my (%pid, %header, %ntseq);
#my $b=0;
my (@l, @l1, @line, @l2, @l3, @s);
my %rseq;
my ($res, $size, $tmpseq);
my %fex; ##saving the phase of the first exon

my $outfile = "$exons_db_folder/$species/$species"."_annot_exon_sizes.txt";
open (OUT, ">$outfile") || die "Cannot create $outfile (loop 7)\n";
    
verbPrint("Generating $species"."_annot_exon_sizes.txt\n");

    ### Parsing GTF file
    #my $gtf_file = "$gtf_folder/$species"."_annot.gtf";
if ($gtf_file =~ /.gz$/) {
    open(GTF, "gunzip -c $gtf_file |") || die "t open pipe to $gtf_file";
}
else {
    open(GTF, $gtf_file) || die "t open $gtf_file";
}
while (<GTF>){
    chomp($_); 
    @line=split(/\t/,$_);
    if ($line[2] eq "CDS"){
	@l1=split(/\"/,$line[8]);
	$gid=$l1[1]; ##change here 
	if ($line[8]=~/protein_id/){
	    $tmp="protein_id";
	} 
	else { $tmp="transcript_id"; }
	@l2=split(/$tmp/,$line[8]);
	@l3=split(/\"/,$l2[1]);
	$prot=$l3[1];
	if (!$pid{$prot}){
	    $m=0;
	    if ($line[7] ne "."){  $phase=$line[7]; $fex{$prot}=$phase; $m=$phase; } ##
	    else { 
		$phase=0; 
		$fex{$prot}=$phase;
	    }
	    $nuc=0;
	    $size=$line[4]-$line[3]+1;
	    $size=$size-$m;
	    $res=$size%3;
	    $nuc=$nuc+$size;
	    $pid{$prot}=1;
	    $header{$prot}="$prot|$gid $size.$phase";
	}
	else {
	    if ($line[7] ne "."){  $phase=$line[7]; }
	    else {  
		if ($res==0){ $phase=0; } elsif ($res==1 ){ $phase=2; } elsif ($res==2) { $phase=1; } 
	    }
	    $size=($line[4]-$line[3])+1;
	    $nuc=$nuc+$size;
	    $header{$prot}.=" $size.$phase";
	}
    }
}
close (GTF);
my @keys=keys(%header);
my $el;
my $p;
my ($l, $triplet, $protein, $nseq, $name);
my $c=0;
foreach $el(@keys){ 
    $name=$header{$el};
    print OUT "$name\n";
}
close (OUT);
verbPrint ">>> All annotation files generated for $species\n";
#A. Removing multi-skipping exon transcripts
my $tint = "$exons_db_folder/$species/$species"."_protein_ids_intron_pos_CDS.txt";
my $oint = "$exons_db_folder/$species/$species"."_CDS_introns.txt";
open (TMP, "$tint");
open (TOUT, ">$oint");
my %ints;
while (<TMP>){
    chomp($_);
    $_=~s/\|/\t/;
    my @ln=split(/\t/,$_);
    if ($ln[4] ne "Intronless"){
	my $tid=$ln[1]."\t".$ln[5];
	$ints{$tid}=1;
    }
}
my $mi;
my @ki=sort(keys(%ints));
foreach $mi (@ki){
    print TOUT "$mi\n";
}
close (TOUT);
#`cat $tint | grep -v Intronless | perl -n -e '$_=~s/\|/\t/; print "$_";' | cut -f2,6 | sort -k1 | uniq > $oint`;
my $ovint= "$exons_db_folder/$species/$species"."_overlap_CDS_introns.txt";
##A.1 Getting overlapping introns
verbPrint ("Getting overlapping introns\n"); 
open (OUT, ">$ovint");
open (INFILE, $oint);
my (@pos);
my (%junctions,%ov_juncs);
my ($sum_junc,$flag,$pos5p,$pos3p,$count,$bandera)=(0,0,0,0,0,0);
my ($gene,$id);
while (<INFILE>){ #AT1G010103913,3996196intronic_reads
    if($_){
	chomp($_);
	@line=split(/\t/,$_); 
	$line[0]=~s/\s+//; ###saving_gene_id
	$line[1]=~s/\s+//;
	@pos=split(/\-/,$line[1]);
	if(!$junctions{$line[0]}){
	    if ($pos5p!=0) {
		print OUT "$ov_juncs{$id}\n";
	    }
	    $count++;
	    $id="OV_INT_".$species."_".$count;
	    $ov_juncs{$id}=$id."\t".$_;
	    $gene=$line[0];
	    $junctions{$line[0]} = $line[1];
	    $sum_junc=1;
	    $flag=1;
	    $pos5p=$pos[0];
	    $pos3p=$pos[1]; 
	}
	else {
	    if ($pos[0]>=$pos5p && $pos[0]<=$pos3p){ ##the junction overlaps
		if ($pos[1]>$pos3p){
		    $pos3p=$pos[1];
		}
		$sum_junc++;
		$flag=0;
		$ov_juncs{$id}.="\n".$id."\t".$_;
	    }
	    else{ ### the junction does not overlap, print all the information from previous junctions
		$sum_junc=0;
		$pos5p=$pos[0];
		$pos3p=$pos[1];
		$sum_junc=1;
		$flag=0;
		$bandera=0;
		print OUT "$ov_juncs{$id}\n";
		$count++;
		$id="OV_INT_".$species."_".$count;
		$ov_juncs{$id}=$id."\t".$_;
	    }
	}
    }
}
close (INFILE);
print OUT "$ov_juncs{$id}\n";
close (OUT);
##A.2 Getting multi-exon introns
verbPrint ("Getting multi-exon introns\n"); 
my (@l,@k,@ints,@i1,@i2);
my (%g,%cov, %ov,%ovint);
my ($el, $r, $s,$id);
open (INONE,"$ovint");
while (<INONE>){
    chomp($_);
    @l=split(/\t/,$_);
    $g{$l[0]}=$l[1];
    if (!$ov{$l[0]}){
	$cov{$l[0]}++;
	$ov{$l[0]}=$l[2];
    }else {
	$cov{$l[0]}++;
	$ov{$l[0]}.=",".$l[2];
    }
}
my @k=keys(%cov);
foreach $el (@k){
    if ($cov{$el}>=4){
	@ints=split(/\,/,$ov{$el});
	for ($r=0; $r<scalar(@ints); $r++){
	    for ($s=0; $s<scalar(@ints); $s++){
		if ($r!=$s){
		    @i1=split(/\-/,$ints[$r]);
		    @i2=split(/\-/,$ints[$s]);
		    if ($i2[0]>$i1[0] && $i2[1]<$i1[1]){
			$id=$el."\t".$g{$el}."\t".$ints[$r];
			$ovint{$id}++;
		    }
		}
	    }
	}

    }
}
my $omul="$exons_db_folder/$species/$species"."_multex_introns.tab";
open (OUT, ">$omul");
my @k=sort(keys(%ovint));
foreach $el (@k){
    if ($ovint{$el}>=3){
	print OUT "$el\t$ovint{$el}\n";
    }
}

##A.3 Removing multiexon skipping transcripts
verbPrint ("Removing multiexon skipping transcripts\n"); 
my (%ints);
open (INONE,"$omul");
while (<INONE>){
    chomp($_);
    @l=split(/\t/,$_);
    if (!$ints{$l[1]}){
	$ints{$l[1]}=$l[2];
    }else {
	$ints{$l[1]}.="|".$l[2];
    }
}
my $ti="$exons_db_folder/$species/$species"."_tr_coords_CDS.txt";
my (@n,@in,@ex,@e1,@i1);
my (%trs,%mex);
my ($s,$idint,$cex,$r);
open (INTWO,"$ti");
while (<INTWO>){
    chomp($_);
    @l=split(/\t/,$_);
    @n=split(/\|/,$l[0]);
    @in=();
    if ($ints{$n[1]}){
	if ($l[3]=~/\-/){
	    @ex=split(/\-/,$l[3]);
	    if ($ints{$n[1]}=~/\|/){
		@in=split(/\|/,$ints{$n[1]});
	    }else { push (@in, $ints{$n[1]}); }
	    for ($s=0; $s<scalar(@in); $s++){
		@i1=split(/\-/,$in[$s]);
		$idint=$n[1]."\t".$in[$s];
		$cex=0;
		if ($l[3]=~/$in[$s]/){ 
		    if (!$trs{$idint}){  $trs{$idint}=$l[0]; }
		    else { $trs{$idint}.=",".$l[0]  }
		}
		for ($r=0; $r<scalar(@ex); $r++){
		    @e1=split(/\,/,$ex[$r]);
		    if ($e1[0]>$i1[0] && $e1[1]<$i1[1]){
			$cex++;
		    }
		}
		#print "$idint\t#$cex#\n";
		if ($cex>=3){
		    $mex{$idint}=1;
		}
	    }
	}
    }
}

my %rmtr;
my @k=keys(%mex);
my @t;
foreach $el(@k){
    if ($trs{$el}){
	#print "#$el#\t$trs{$el}\n"; 
	if ($trs{$el}=~/\,/){
	    @t=split(/\,/,$trs{$el});
	    foreach $m(@t){
		$rmtr{$t[$m]}=1;
	    }
	}else { $rmtr{$trs{$el}}=1;  }
    }
}

my $to="$exons_db_folder/$species/$species"."_multiexsk_trs.txt"; 
print "$to\n";
open (OUT, ">$to");
open (INTWO,"$ti");
while (<INTWO>){
    chomp($_);
    @l=split(/\t/,$_);
    if(!$rmtr{$l[0]}){
	print OUT "$_\n";
    }
}
my $nrm=scalar(keys(%rmtr));
verbPrint ("Number of transcripts to be removed: $nrm\n");
verbPrint ("Replacing old file of transcripts for the filtered one\n");
##replacing old file of transcripts for the filtered one
`mv $to $ti`;
verbPrint ("Annotations for $species finished!!!..."); 