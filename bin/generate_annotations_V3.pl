#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path cwd);

#Declaration of variables
my $gtf_file;
my $genome_file;
my $exons_db_folder="./";
my $species;
my $verboseFlag=1;
my $help;
my $add_exons;
my $do_all_steps=1;

#Arguments
Getopt::Long::Configure("no_auto_abbrev");
GetOptions("GTF=s" => \$gtf_file,
	   "G=s" => \$genome_file,
	   "EX_DB=s" => \$exons_db_folder,
	   "sp=s" => \$species,
	   "add_exons=s" => \$add_exons, #before $vastdb_refs
	   "verbose=s" => \$verboseFlag,
	   "h" => \$help,
	   "help" => \$help
    );


sub verbPrint {
    my $verbMsg = shift;
    unless ($verboseFlag == 0 || $verboseFlag eq "F" || $verboseFlag eq "FALSE") {
	chomp($verbMsg);
	print STDERR "[ExOrthist annotation]: $verbMsg\n";
    }
}
if (!defined $genome_file || !defined $gtf_file || !defined $species || defined $help){
    die "\nUsage: generate_annotations.pl -GTF path_to_gtfs/ -G path_to_genomes/ -sp Sp1 [-EX_DB path_to_EXONS_DB/ -add_exons exon_file ]

Script that creates all annotation files needed for the second module of the pipeline

OMPULSORY
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
system "mkdir $exons_db_folder/$species" unless (-e "$exons_db_folder/$species");
verbPrint("EXONS_DB path set to $full_path_exons_db\n");

### Future GTF file with original + fake transcripts (if extra exons are provided)
my $merged_GTF = $exons_db_folder."/".$species."/".$species."_annot_fake.gtf";

### Loading sequence info
my (%seq);
my $chr;
my @g;
if ($genome_file =~ /.gz$/) {
    open(GENOME, "gunzip -c $genome_file |") || die "It cannot open pipe to $genome_file";
}
else {
    open(GENOME, $genome_file) || die "It cannot open $genome_file";
}
verbPrint("Parsing gDNA file for $species\n");
while (<GENOME>){
    chomp($_);
    if ($_=~/\>/){ $_=~s/\>//; @g=split(/\s+/,$_); $chr=$g[0]; }
    else { $seq{$chr}.=$_;  }    
}
close (GENOME);

### PART 1: get_ref_proteins.pl and getting exint file
if ($do_all_steps){
    verbPrint("Generating exint file for $species\n");
    my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
    my (%pid, %header, %ntseq);
    my (@l, @l1, @line, @l2, @l3, @s);
    my %rseq;
    my ($res, $size, $tmpseq);
    my %fex; ##saving the phase of the first exon
    
#### Opening and reading GTF file
    if ($gtf_file =~ /.gz$/) {
	open(GTF, "gunzip -c $gtf_file |") || die "It cannot open pipe to $gtf_file";
    }else {
	open(GTF, $gtf_file) || die "It cannot open $gtf_file";
    }
    while (<GTF>){
	chomp($_); 
	@line=split(/\t/,$_);
	if ($seq{$line[0]}){
	    if ($line[2] eq "CDS"){
		@l1=split(/\"/,$line[8]);
		$gid=$l1[1]; 
		if ($line[8]=~/protein_id/){
		    $tmp="protein_id";
		} 
		else {$tmp="transcript_id"; }
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
    my $p;
    my ($l, $triplet, $protein, $nseq, $name);
    my $c=0;
### output files (former -F 1 run)
    my $exint_output = "$exons_db_folder/$species/$species.exint";
    open (EXINT_OUT, ">$exint_output") || die "Cannot open $exint_output (loop 1)\n";
    foreach my $el (@keys){ 
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
		    elsif ($codons{$triplet} eq "stop"){ $c=1; } 
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

    ### Creates files with the longest protein per gene (exint not maintained, just list)
    my $longest_prot_output = "$exons_db_folder/$species/$species"."_ref_proteins.txt";
#    my $longest_exint_output = "$exons_db_folder/$species/$species"."_ref_proteins.exint";
    open (L_P_OUTPUT, ">$longest_prot_output") || die "Cannot open $longest_prot_output (loop 1)\n";
#    print L_P_OUTPUT "GeneID\tProteinID\n"; # maintain without header, as the original
#    open (L_E_OUTPUT, ">$longest_exint_output") || die "Cannot open $longest_exint_output (loop 1)\n";
    foreach my $temp_exint (sort keys %longest_prot){
	print L_P_OUTPUT "$temp_exint\t$longest_prot{$temp_exint}\n";
#	print L_E_OUTPUT "$exint_data{$longest_prot{$temp_exint}}\n";
    }
    close L_P_OUTPUT;
#    close L_E_OUTPUT;
    
##### Introducing missing exons in GTF if extra exons are provided
    if (defined $add_exons){ 
	verbPrint("Introducing additional exons for $species\n");
	#ENSG00000029534  HsaEX0004110  chr8:41758036-41758137  chr8:41797512-41797622  chr8:41733971-41734069  chr8:41797512,41758036-41758137,41734069:-
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
	    $A_ref{$ev}=$t[2];
	    $C1_ref{$ev}=$t[3];
	    $C2_ref{$ev}=$t[4];
	    $C1_inc{$ev}=$t[3]; ##same as C1 reference
	    $C2_inc{$ev}=$t[4]; ##same as C2 reference 	
	    $co_ev{$t[2]}=$ev;
	    ($chr,$i,$f)=$t[2]=~/(.+?)\:(.+?)\-(.+)/;
	    $coA="$chr:$i";
	    $coB="$chr:$f";
	    $coA_ev{$coA}=$ev;
	    $coB_ev{$coB}=$ev;
	}
	close EXONS;
	
	my ($tr, $co,$exon_number,$g_2,$tr_2,);
	my (%gene_name,%annotated,%partial,%array_tr_co,%array_tr_coA,%array_tr_coB,%index_coA,%index_coB,%TR,%index_co,%tr_co,%tr_coA,%tr_coB,%str,%done);
	my (%cds_lines,%cds_ini,%cds_end,%offset,%start_lines,%stop_lines,%tr_lines);
	
#### Starts checking for annotation
	if ($gtf_file =~ /.gz$/) {
	    open(GTF, "gunzip -c $gtf_file |") || die "It cannot open pipe to $gtf_file";
	}else {
	    open(GTF, $gtf_file) || die "It cannot open $gtf_file";
	}
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
	# A lot variable declarations to avoid issues from merging scripts
	my ($tally_annotated,$C1_inc,$C2_inc,$C1_ref,$C2_ref,$C1i,$C1r,$C1_incA,$C1_refA,$C1f,$C2f,$C2_incA,$C2_refA,$C2i,$C2r,$Ai,$Af,$coA_A,$coB_A,$temp_tr,$exons_tr,$selected_tr,$ref_tr,$OK,$t_C2,$exons_sel,$t_C1,$CDS_end,$CDS_ini,$CDS_seq,$STOP_detected,$bit,$co_Af_CDS,$co_Ai_CDS,$codon,$comment,$exN,$extra,$final_tr,$finished,$first_offset,$full_line,$full_line2,$full_line3,$full_line4,$le_CDS_ex,$length_A,$lineN,$loop,$new_line,$new_line2,$new_line3,$new_offset,$offset_A,$out_frame,$rescue_type,$seq_bit,$seq_bit_test,$started,$stop_co_f,$stop_co_i,$str,$tally_non_annot_hit,$tally_non_annot_no_hit,$tally_partial,$type);
	my (@t_line,@t_line2,@t_line3,@t_line4,@t_line_C2);
	my (%C1_accepted_index,%tally_solutions);
	foreach $ev (sort keys %A_ref){
	    $g=$Ev_Gene{$ev};
	    if ($annotated{$ev}){
		$tally_annotated++;
	    }
	    else {
		$selected_tr="";
# 		my $ref_tr; # not provided 
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
		
		foreach $tr (sort keys %{$TR{$g}}){
		    next if $tr_coA{$tr}{$coA_A} || $tr_coB{$tr}{$coB_A}; # if the transcript contains any junction of the A exon
		    
		    if ($tr_co{$tr}{$C1_inc} && $tr_co{$tr}{$C2_inc} && $index_co{$tr}{$C1_inc}==$index_co{$tr}{$C2_inc}-1){
			($temp_tr)=$selected_tr=~/1\=(.+)/; # it may be empty
			if (defined $temp_tr){ 
			    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
			}else {$exons_sel=1;}
			$exons_tr=$#{$cds_lines{$tr}}+1;
			if (defined $ref_tr){
			    $selected_tr="1=$tr" unless ($selected_tr eq "1=$ref_tr" || $exons_sel>=$exons_tr);
			} else {$selected_tr="1=$tr" unless ($exons_sel>=$exons_tr);}
			$C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_inc};
		    }
		    elsif ($tr_co{$tr}{$C1_inc} && $tr_co{$tr}{$C2_ref} && $index_co{$tr}{$C1_inc}==$index_co{$tr}{$C2_ref}-1){
			($temp_tr)=$selected_tr=~/[12]\=(.+)/; # it may be empty
			if (defined $temp_tr){ 
			    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
			}else {$exons_sel=1;}
			$exons_tr=$#{$cds_lines{$tr}}+1;
			if (defined $ref_tr){
			    $selected_tr="2=$tr" unless ($selected_tr eq "2=$ref_tr" || $selected_tr=~/1\=/ || $exons_sel>=$exons_tr);
			} else {$selected_tr="2=$tr" unless ($selected_tr=~/1\=/ || $exons_sel>=$exons_tr);}
			$C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_inc};
		    }
		    elsif ($tr_co{$tr}{$C1_ref} && $tr_co{$tr}{$C2_inc} && $index_co{$tr}{$C1_ref}==$index_co{$tr}{$C2_inc}-1){
			($temp_tr)=$selected_tr=~/[123]\=(.+)/; # it may be empty
			if (defined $temp_tr){
			    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
			}else {$exons_sel=1;}
			$exons_tr=$#{$cds_lines{$tr}}+1;
			if (defined $ref_tr){
			    $selected_tr="3=$tr" unless ($selected_tr eq "3=$ref_tr" || $selected_tr=~/[12]\=/ || $exons_sel>=$exons_tr);
			} else {$selected_tr="3=$tr" unless ($selected_tr=~/[12]\=/ || $exons_sel>=$exons_tr);}
			$C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_ref};
		    }
		    elsif ($tr_co{$tr}{$C1_ref} && $tr_co{$tr}{$C2_ref} && $index_co{$tr}{$C1_ref}==$index_co{$tr}{$C2_ref}-1){
			($temp_tr)=$selected_tr=~/[1234]\=(.+)/; # it may be empty
			if (defined $temp_tr){
			    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
			}else {$exons_sel=1;}
			$exons_tr=$#{$cds_lines{$tr}}+1;
			if (defined $ref_tr){
			    $selected_tr="4=$tr" unless ($selected_tr eq "4=$ref_tr" || $selected_tr=~/[123]\=/ || $exons_sel>=$exons_tr);
			} else {$selected_tr="4=$tr" unless ($selected_tr=~/[123]\=/ || $exons_sel>=$exons_tr);}
			$C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_ref};
		    }
		    elsif ($tr_co{$tr}{$C1_inc}){
			($temp_tr)=$selected_tr=~/[12345]\=(.+)/; # it may be empty
			if (defined $temp_tr){
			    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
			}else {$exons_sel=1;}
			$exons_tr=$#{$cds_lines{$tr}}+1;
			if (defined $index_co{$tr}{$C1_inc}){
			    $t_C2=$array_tr_co{$tr}[$index_co{$tr}{$C1_inc}+1];
			    ($i,$f)=$t_C2=~/\:(.+?)\-(.+)/;
			    $OK="";
			    $OK=1 if $i>$Af && $str{$g} eq "+";
			    $OK=1 if $f<$Ai && $str{$g} eq "-";
			    if (defined $ref_tr){
				$selected_tr="5=$tr" unless ($selected_tr eq "5=$ref_tr" || $selected_tr=~/[1234]\=/ || $exons_sel>=$exons_tr || !$OK);
			    } else {$selected_tr="5=$tr" unless ($selected_tr=~/[1234]\=/ || $exons_sel>=$exons_tr || !$OK);}
			    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_inc};
			}
		    }
		    elsif ($tr_co{$tr}{$C1_ref}){
			($temp_tr)=$selected_tr=~/[123456]\=(.+)/; # it may be empty
			if (defined $temp_tr){
			    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
			}else {$exons_sel=1;}
			$exons_tr=$#{$cds_lines{$tr}}+1;
			if (defined $index_co{$tr}{$C1_ref}){
			    $t_C2=$array_tr_co{$tr}[$index_co{$tr}{$C1_ref}+1];
			    ($i,$f)=$t_C2=~/\:(.+?)\-(.+)/;
			    $OK="";
			    $OK=1 if $i>$Af && $str{$g} eq "+";
			    $OK=1 if $f<$Ai && $str{$g} eq "-";
			    if (defined $ref_tr){
				$selected_tr="6=$tr" unless ($selected_tr eq "6=$ref_tr" || $selected_tr=~/[12345]\=/ || $exons_sel>=$exons_tr || !$OK);
			    } else {$selected_tr="6=$tr" unless ($selected_tr=~/[12345]\=/ || $exons_sel>=$exons_tr || !$OK);}
			    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C1_ref};
			}
		    }
		    elsif (($tr_coA{$tr}{$C1_incA} && $str{$g} eq "-") || ($tr_coB{$tr}{$C1_incA} && $str{$g} eq "+") ){ # C1do exists in transcript
			($temp_tr)=$selected_tr=~/[1234567]\=(.+)/; # it may be empty
			if (defined $temp_tr){
			    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
			}else {$exons_sel=1;}
			$exons_tr=$#{$cds_lines{$tr}}+1;
			if (defined $index_coA{$tr}{$C1_incA}){
			    $t_C2=$array_tr_co{$tr}[$index_coA{$tr}{$C1_incA}+1]; # gets a proper C2, not an acceptor
			    ($i,$f)=$t_C2=~/\:(.+?)\-(.+)/;
			    $OK="";
			    $OK=1 if $i>$Af && $str{$g} eq "+";
			    $OK=1 if $f<$Ai && $str{$g} eq "-";
			    if (defined $ref_tr){
				$selected_tr="7=$tr" unless ($selected_tr eq "7=$ref_tr" || $selected_tr=~/[123456]\=/ || !$OK || $exons_sel>=$exons_tr);
			    } else {$selected_tr="7=$tr" unless ($selected_tr=~/[123456]\=/ || !$OK || $exons_sel>=$exons_tr);}
			    $C1_accepted_index{$ev}{$tr}=$index_coA{$tr}{$C1_incA};
			}
		    }
		    elsif (($tr_coA{$tr}{$C1_refA} && $str{$g} eq "-") || ($tr_coB{$tr}{$C1_refA} && $str{$g} eq "+")){ # C1do exists in transcript
			($temp_tr)=$selected_tr=~/[12345678]\=(.+)/; # it may be empty
			if (defined $temp_tr){
			    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
			}else {$exons_sel=1;}
			$exons_tr=$#{$cds_lines{$tr}}+1;
			if (defined $index_coA{$tr}{$C1_refA}){
			    $t_C2=$array_tr_co{$tr}[$index_coA{$tr}{$C1_refA}+1]; # gets a proper C2, not an acceptor
			    ($i,$f)=$t_C2=~/\:(.+?)\-(.+)/;
			    $OK="";
			    $OK=1 if $i>$Af && $str{$g} eq "+";
			    $OK=1 if $f<$Ai && $str{$g} eq "-";
			    if (defined $ref_tr){
				$selected_tr="8=$tr" unless ($selected_tr eq "8=$ref_tr" || $selected_tr=~/[1234567]\=/ || !$OK || $exons_sel>=$exons_tr);
			    } else {$selected_tr="8=$tr" unless ($selected_tr=~/[1234567]\=/ || !$OK || $exons_sel>=$exons_tr);}
			    $C1_accepted_index{$ev}{$tr}=$index_coA{$tr}{$C1_refA};
			}
		    }
		    elsif ($tr_co{$tr}{$C2_inc}){
			($temp_tr)=$selected_tr=~/[123456789]\=(.+)/; # it may be empty
			if (defined $temp_tr){
			    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
			}else {$exons_sel=1;}
			$exons_tr=$#{$cds_lines{$tr}}+1;
			if (defined $index_co{$tr}{$C2_inc}){
			    $t_C1=$array_tr_co{$tr}[$index_co{$tr}{$C2_inc}-1];
			    ($i,$f)=$t_C1=~/\:(.+?)\-(.+)/;
			    $OK="";
			    $OK=1 if $f<$Ai && $str{$g} eq "+";
			    $OK=1 if $i>$Af && $str{$g} eq "-";
			    if (defined $ref_tr){
				$selected_tr="9=$tr" unless ($selected_tr eq "9=$ref_tr" || $selected_tr=~/[12345678]\=/ || !$OK || $exons_sel>=$exons_tr);
			    } else {$selected_tr="9=$tr" unless ($selected_tr=~/[12345678]\=/ || !$OK || $exons_sel>=$exons_tr);}
			    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C2_inc}-1;
			}
		    }
		    elsif ($tr_co{$tr}{$C2_ref}){
			($temp_tr)=$selected_tr=~/[1234567890]\=(.+)/; # it may be empty
			if (defined $temp_tr){
			    $exons_sel=$#{$cds_lines{$temp_tr}}+1; # it may be empty
			}else {$exons_sel=1;}
			$exons_tr=$#{$cds_lines{$tr}}+1;
			if (defined $index_co{$tr}{$C2_ref}){
			    $t_C1=$array_tr_co{$tr}[$index_co{$tr}{$C2_ref}-1];
			    ($i,$f)=$t_C1=~/\:(.+?)\-(.+)/;
			    $OK="";
			    $OK=1 if $f<$Ai && $str{$g} eq "+";
			    $OK=1 if $i>$Af && $str{$g} eq "-";
			    if (defined $ref_tr){
				$selected_tr="10=$tr" unless ($selected_tr eq "10=$ref_tr" || $selected_tr=~/[123456789]\=/ || !$OK || $exons_sel>=$exons_tr);
			    } else {$selected_tr="10=$tr" unless ($selected_tr=~/[123456789]\=/ || !$OK || $exons_sel>=$exons_tr);}
			    $C1_accepted_index{$ev}{$tr}=$index_co{$tr}{$C2_ref}-1;
			}
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
				($extra)=$CDS_seq=~/(..)$/;
				$extra = "NN" if !defined $extra;
			    }
			    elsif ($new_offset==2){
				($extra)=$CDS_seq=~/(.)$/;
				$extra = "N" if !defined $extra;
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
				    ($extra)=$CDS_seq=~/(..)$/;
				    $extra="NN" if !defined $extra;
				}
				elsif ($offset_A==2){
				    ($extra)=$CDS_seq=~/(.)$/;
				    $extra="N" if !defined $extra;
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
		}
	    }
	}
	verbPrint("   Fully Annotated\t$tally_annotated\n");
	verbPrint("   Rescued skipping\t$tally_non_annot_hit\n");
	verbPrint("   Partially Annotated\t$tally_partial\n");
	verbPrint("   Not rescued\t$tally_non_annot_no_hit\n");
	print LOG2 "Fully Annotated\t$tally_annotated\n".
	    "Rescued skipping\t$tally_non_annot_hit\n".
	    "Partially Annotated\t$tally_partial\n".
	    "Not rescued\t$tally_non_annot_no_hit\n";
	
	verbPrint("   Solutions:\n");
	print LOG2 "\nSolutions:\n";
	foreach $type (sort {$a<=>$b} keys %tally_solutions){
	    verbPrint("   Type $type\t$tally_solutions{$type}\n");
	    print LOG2 "Type $type\t$tally_solutions{$type}\n";
	}
	close LOG2;
	close LOG;
	close O;
	
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
	close (FAKEGTF);
	
#2) Get exint file for fake transcripts
	my @keys2=keys(%ntseq2);
	verbPrint("Generating exint files with added exons for $species\n");
	open (EXINT_OUT, ">$tmp_exint_output") || die "Cannot open $exint_output (loop 1)\n";
	foreach my $el(@keys2){ 
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
	if ($gtf_file=~/\.gz$/){ # if GTF is compressed
	    `gunzip -c $gtf_file > $merged_GTF`;
	    `cat $fakeGTF >> $merged_GTF`;
	}
	else {
	    `cat $gtf_file $fakeGTF > $merged_GTF`;
	}
	
###Joining 2. Exint files
	`cat $tmp_exint_output >> $exint_output`;
	`rm $tmp_exint_output`;
    }
    else { # if no extra exons added, it creates a copy of the original GTF
	if ($gtf_file=~/\.gz$/){ # if GTF is compressed
	    `gunzip -c $gtf_file > $merged_GTF`;
	}
	else {
	    `cp $gtf_file $merged_GTF`;
	}
    }
}

### PART 2: get_trs_gtf.pl
if ($do_all_steps){ # done to re-declare variables
    my (@line, @coords1, @coords2, @gene, @l1, @l2, @l3);
    my ($l, $grep, $gid, $int, $r, $prot, $tmp, $trid, $n);
    my (%tr1, %tr2, %tpid);
    my %prot;
    my (%coords1, %coords2);
    my (%strand1, %strand2);
    my %intron;
    my (%chr1, %chr2);
    
    open (INFILEONE, $merged_GTF) || die "It cannot open $merged_GTF\n";
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
}

### PART 3: get_aa_pos_exon.pl
if ($do_all_steps){ # done to re-declare variables
#   -i Hsa_annot_exons_prot_ids.#   -out Hsa_protein_ids_exons_pos.txt
#Declaration of variables
    my (%pid, %protein);
    my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m, $ex);
    my (%header, %ntseq);
    my (@l, @l1, @line, @l2, @l3, @s);
    
# Declaration of input files
    my $p = "$exons_db_folder/$species/$species"."_annot_exons_prot_ids.txt";
    my $v = "$exons_db_folder/$species/$species"."_add_exons_prot_ids.txt" if (-e "$exons_db_folder/$species/$species"."_add_exons_prot_ids.txt");
    
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
### if VASTDB refs were provided => DEPRECATED
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
    open (GTF, $merged_GTF) || die "It cannot open $merged_GTF\n";
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
}

### PART 4: get_genomic_coaords_intron.pl
if ($do_all_steps){ # done to re-declare variables
#   -p Hsa_annot_exons_prot_ids.txt -v Hsa_vast_exons_prot_ids.txt -t Hsa_trid_protid.txt 
#   -c1 Hsa_tr_coords.txt -c2 Hsa_tr_coords_CDS.txt 
#   -o1 Hsa_protein_ids_intron_pos.txt -o2 Hsa_protein_ids_intron_pos_CDS.txt
# Declaration of variables
    my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
    my (%pid, %tid, %ntseq);
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
# parses input 2 => DEPRECATED
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
# Format: ENSBTAT00000000005  ENSBTAP00000000005
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
}

### PART 5: get_ex_size_phase.pl
if ($do_all_steps){
# -GTF Hsa_annot.gtf -out Hsa_annot_exon_sizes.txt
# Declaration of variables (same as loop 1)
    my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
    my (%pid, %header, %ntseq);
    my (@l, @l1, @line, @l2, @l3, @s);
    my %rseq;
    my ($res, $size, $tmpseq);
    my %fex; ##saving the phase of the first exon
    
    my $outfile = "$exons_db_folder/$species/$species"."_annot_exon_sizes.txt";
    open (OUT, ">$outfile") || die "Cannot create $outfile (loop 7)\n";
    
    verbPrint("Generating $species"."_annot_exon_sizes.txt\n");
### Parsing GTF file
    open(GTF, $merged_GTF) || die "It cannot open $merged_GTF\n";
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
    my $p;
    my ($l, $triplet, $protein, $nseq, $name);
    my $c=0;
    foreach my $el (@keys){ 
	$name=$header{$el};
	print OUT "$name\n";
    }
    close (OUT);
}

### PART 6: Removing multi-skipping exon transcripts
if ($do_all_steps){
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
    close TMP;

    my $mi;
    my @ki=sort(keys(%ints));
    foreach $mi (@ki){
	print TOUT "$mi\n";
    }
    close (TOUT);
    
    ## SUBPART 6.1 Getting overlapping introns
    verbPrint ("Getting overlapping introns\n"); 
    my $ovint= "$exons_db_folder/$species/$species"."_overlap_CDS_introns.txt";
    open (OUT, ">$ovint");
    open (INFILE, $oint);
    my (@pos,@line);
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
}

### SUBPART 6.2 Getting multi-exon introns
if ($do_all_steps){
    verbPrint ("Getting multi-exon introns\n"); 
    my (@l,@ints,@i1,@i2);
    my (%g,%cov, %ov,%ovint);
    my ($r, $s,$id);
    my $ovint= "$exons_db_folder/$species/$species"."_overlap_CDS_introns.txt";
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
    foreach my $el (@k){
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
    # Creates output
    my $omul="$exons_db_folder/$species/$species"."_multex_introns.tab";
    open (OUT, ">$omul");
    my @k2=sort(keys(%ovint));
    foreach my $el (@k2){
	if ($ovint{$el}>=3){
	    print OUT "$el\t$ovint{$el}\n";
	}
    }
    close OUT;
}

### SUBPART 6.3: Removing multiexon skipping transcripts
if ($do_all_steps){
    verbPrint ("Removing multiexon skipping transcripts\n"); 
    my (%ints);
    my $omul="$exons_db_folder/$species/$species"."_multex_introns.tab";
    open (INONE,"$omul");
    while (<INONE>){
	chomp($_);
	my @l=split(/\t/,$_);
	if (!$ints{$l[1]}){
	    $ints{$l[1]}=$l[2];
	}else {
	    $ints{$l[1]}.="|".$l[2];
	}
    }
    close INONE;
    
    my $f_transc_coords_CDS = "$exons_db_folder/$species/$species"."_tr_coords_CDS.txt"; # for $ti
    my (@in,@ex,@e1,@i1);
    my (%trs,%mex);
    my ($s,$idint,$cex,$r);
    open (INTWO, $f_transc_coords_CDS) || die "It cannot open file with tr_coords_CDS ($f_transc_coords_CDS)\n";
    while (<INTWO>){
	chomp($_);
	my @l=split(/\t/,$_);
	my @n=split(/\|/,$l[0]);
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
		    if ($cex>=3){
			$mex{$idint}=1;
		    }
		}
	    }
	}
    }
    close (INTWO);
    
    my %rmtr;
    my @k=keys(%mex);
    my @t;
    foreach my $el (@k){
	# $el => FBgn0003429  24698269-24703021
	# $trs{$el} => FBpp0100093|FBgn0003429
	if ($trs{$el}){
	    if ($trs{$el}=~/\,/){ # not needed...
		@t=split(/\,/,$trs{$el});
		foreach my $tmp_m (@t){
		    $rmtr{$tmp_m}=1;
		}
	    }
	    else { 
		$rmtr{$trs{$el}}=1;  
	    }
	}
    }
    
    ##### Removes multi-skipping transcripts from files
    # Sp_multiexsk_trs.txt
    my $temp_output = "$exons_db_folder/$species/$species"."_multiexsk_trs.txt"; 
    open (OUT, ">$temp_output") || die "It cannot open the temporary output without multi-sk transcripts\n";
    open (IN, $f_transc_coords_CDS) || die "It cannot open file with tr_coords_CDS ($f_transc_coords_CDS)\n";
    while (<IN>){
	chomp($_);
	my @l=split(/\t/,$_);
	if (!defined $rmtr{$l[0]}){
	    print OUT "$_\n";
	}
    }
    close IN;
    close OUT;
    system "mv $temp_output $f_transc_coords_CDS";
    # Sp_annot_exons_prot_ids.txt
    my $input_file = "$exons_db_folder/$species/$species"."_annot_exons_prot_ids.txt";
    open (OUT, ">$input_file-temp") || die "It cannot open the temporary output without multi-sk transcripts\n";
    open (IN, $input_file) || die "It cannot open file with annot_exons_prot_ids ($input_file)\n";
    while (<IN>){ #Format: FBpp0309061|FBgn0025836
	chomp($_);
	my @l=split(/\t/,$_);
	if (!defined $rmtr{$l[0]}){
	    print OUT "$_\n";
	}
    }
    close IN;
    close OUT;
    system "mv $input_file-temp $input_file";
    # Sp_protein_ids_exons_pos.txt
    $input_file = "$exons_db_folder/$species/$species"."_protein_ids_exons_pos.txt";
    open (OUT, ">$input_file-temp") || die "It cannot open the temporary output without multi-sk transcripts\n";
    open (IN, $input_file) || die "It cannot open file with annot_exons_prot_ids ($input_file)\n";
    while (<IN>){ #Fomat: FBpp0309061|FBgn0025836  exon_11-113  chrX  129082-129417  +
	chomp($_);
	my @l=split(/\t/,$_);
	if (!defined $rmtr{$l[0]}){
	    print OUT "$_\n";
	}
    }
    close IN;
    close OUT;
    system "mv $input_file-temp $input_file";
    # Sp_protein_ids_intron_pos_CDS.txt
    $input_file = "$exons_db_folder/$species/$species"."_protein_ids_intron_pos_CDS.txt";
    open (OUT, ">$input_file-temp") || die "It cannot open the temporary output without multi-sk transcripts\n";
    open (IN, $input_file) || die "It cannot open file with annot_exons_prot_ids ($input_file)\n";
    while (<IN>){ #Fomat: FBpp0309167|FBgn0004456  chrX  -   Intron_1  13225111-13256617
	chomp($_);
	my @l=split(/\t/,$_);
	if (!defined $rmtr{$l[0]}){
	    print OUT "$_\n";
	}
    }
    close IN;
    close OUT;
    system "mv $input_file-temp $input_file";
    # Sp.exint
    $input_file = "$exons_db_folder/$species/$species.exint";
    open (OUT, ">$input_file-temp") || die "It cannot open the temporary output without multi-sk transcripts\n";
    open (IN, $input_file) || die "It cannot open file with annot_exons_prot_ids ($input_file)\n";
    my $header;
    while (<IN>){ # Fomat: >FBpp0088149|FBgn0039911 X*?\nSEQ
	chomp($_);
	if ($_=~/\>(.+)/){
	    $header = $1;
	}
	else {
	    my $seq = $_;
	    if (defined $header){
		my ($protein) = $header =~/(\S+)/;
		if (defined $protein){
		    print OUT ">$header\n$seq\n" if (!defined $rmtr{$protein});
		}
		else {
		    print "ISSUE: $header\n";
		}
	    }
	    else {
		print "ISSUE: $seq\n";
	    }
	}
    }
    close IN;
    close OUT;
    system "mv $input_file-temp $input_file";
    
    my $number_rm_tr = scalar(keys(%rmtr));
    verbPrint ("Number of transcripts with multiexon skipping removed: $number_rm_tr\n");
}
# Handles intermediate files
my $GTF_to_compress = "$exons_db_folder/$species/$species"."_annot_fake.gtf";
system "gzip $GTF_to_compress";

###########################################################################################
###Getting overlapping exons by species####

my $exposfile=$exons_db_folder."/".$species."/".$species."_protein_ids_exons_pos.txt";
my $outexfile=$exons_db_folder."/".$species."/".$species."_prot_exons.txt";
my $finalout=$exons_db_folder."/".$species."/".$species."_prot_exons.txt";
my %junctions;
my %ov_juncs;
my $sum_junc=0;
my $flag=0;
my ($pos5p,$pos3p)=0;
my $bandera=0;
my $count=0;
my $id;

##sorting exons by gene and position
`cat ${s}  | perl -n -e '$_=~s/\|/\t/; print "$_";' | cut -f2,6 |  sort -k1 | uniq > $outexfile`;
open (EXFILE, ">$finalout"); ##Final output of overlapping exons
open (INFILE, $outexfile) || die "It cannot open $outexfile\n";
while (<INFILE>){ #BL00113	Sc0000002:5027516-5027714:+	Bla	18
    if($_){
	chomp($_);
	my @line=split(/\t/,$_);
	$line[0]=~s/\s+//; ###saving_gene_id
	$line[1]=~s/\s+//;
	my @crs=split(/\:/,$line[1]);
	my @pos=split(/\-/,$crs[1]);
	my $species=$line[2];

	if (!$junctions{$line[0]}){
	    if ($pos5p!=0) {
		if (defined $id){
		    print EXFILE "$ov_juncs{$id}\n" if (defined $ov_juncs{$id});
		}
	    }
	    $count++;
	    $id="OV_EX_".$species."_".$count;
	    $ov_juncs{$id}=$id."\t".$_;
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
		print EXFILE "$ov_juncs{$id}\n" if (defined $ov_juncs{$id});
		$count++;
		$id="OV_EX_".$species."_".$count;
		$ov_juncs{$id}=$id."\t".$_;
	    }
	}
    }
}
close (INFILE);
print EXFILE "$ov_juncs{$id}\n" if (defined $ov_juncs{$id});
close EXFILE;

verbPrint ("Annotations for $species finished!"); 

##END OF SCRIPT##






