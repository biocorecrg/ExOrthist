#!/usr/bin/env perl
use warnings;
use strict;

##SCRIPT FOR SCORING EXONS BTW PAIR OF SPECIES##
my $s1=$ARGV[0]; ##species 1 
my $s2=$ARGV[1]; ##species 2
my $i1=$ARGV[2]; ##Gene clusters
my $i2=$ARGV[3]; ##file protein ids exons species 1
my $i3=$ARGV[4]; ##file protein ids exons species 2
my $i4=$ARGV[5]; ##file with exon positions in aa of species 1
my $i5=$ARGV[6]; ##file with exon positions in aa of species 2
my $i6=$ARGV[7]; ##file with intron positions in aa of species 1
my $i7=$ARGV[8]; ##file with intron positions in aa of species 2
my $i8=$ARGV[9];  ##exint file species 1
my $i9=$ARGV[10];  ##exint file species 
my $part=$ARGV[11];
my $bl=$ARGV[12]; ##blosum62  matrix
my $outf=$ARGV[13]; ## output folder
my $min_sim_prots=$ARGV[14]; ##minimum prot sim in decimal format.
my $cpus=$ARGV[15];
my $prev_folder=$ARGV[16]; # OPTIONAL. Folder with subfolders with pairwise pre-computed info.

$cpus=1 if !$cpus;

my ($dir,$project_dir)=$outf=~/(.+)\/(.+?)\/.+?\//;
my ($sp1,$sp2)=($s1,$s2);
my $idex;
my %onehit;
my %miss;
my (@tn1, @tn2, @ln);
system "mkdir $outf" unless (-e $outf);

my $exsc=$outf."/EX_aln_features_".$s1."_".$s2."_part_".$part.".txt";##outputfile for exon scores
my $insc=$outf."/INT_aln_features_".$s1."_".$s2."_part_".$part.".txt";
my $prsc=$outf."/PROT_aln_features_".$s1."_".$s2."_part_".$part.".txt";##outputfile for scores of whole protein alignment
my $msf=$outf."/EXs_to_split_part_".$part.".txt"; ##generating a file for those exons that need to be realigned
open (EXSC, ">$exsc") || die "It cannot open output file $exsc\n";
open (INSC, ">$insc") || die "It cannot open output file $insc\n";
open (PRSC, ">$prsc") || die "It cannot open output file $prsc\n";
open (MISS, ">$msf") || die "It cannot open output file $msf\n";

##Printing headers 
#Header protein score file
print PRSC "CID\tALN\tQuery\tSubject\t\%Sim_qry_sbj\t\%Sim_sbj_qry\t\%Min_Identity\tGlobal_Score\tSp_query\tSp_subject\n";
#Header exon score file
print EXSC "CID\tProt_query\tExon_number_query\tAA_pos_exon_query\tChr_query\tExon_coords_query\tStrand\tExon_hits_subject\tProt_subject\tAln_AA_pos_subject\t\%Sim_aln_qry_sbj\t\%Id_aln_qry_sbj\tGap_number\t\%Gaps\tProt_subject\tExon_number_subject\tAA_pos_exon_subject\tChr_subject\tExon_coords_subject\tStrand\tSp_query\tSp_subject\n";
#Header exons to realign file
print MISS "CID\tProt_query\tExon_number_query\tAA_pos_exon_query\tChr_query\tExon_coords_query\tStrand\tExon_hits_subject\tProt_subject\tAln_AA_pos_subject\t\%Sim_aln_qry_sbj\t\%Id_aln_qry_sbj\tGap_number\t\%Gaps\tProt_subject\tExon_number_subject\tAA_pos_exon_subject\tChr_subject\tExon_coords_subject\tStrand\tSp_query\tSp_subject\n";
#Header intron score file
print INSC "CID\tProt_query\tIntron_number_query\tAA_pos_intron_query\tIntron_phase_query\tProt_subject\tIntron_number_subject\tAA_pos_intron_subject\tIntron_phase_subject\tAln_dev\tAln_score\tIntron_score\tSp_query\tSp_subject\n";
################
my %sim; ##AA similarity scores
$sim{"FY"}=1;
$sim{"YF"}=1;
$sim{"FW"}=1;
$sim{"WF"}=1;
$sim{"YW"}=1;
$sim{"WY"}=1;
$sim{"VI"}=1;
$sim{"IV"}=1;
$sim{"VL"}=1;
$sim{"LV"}=1;
$sim{"LI"}=1;
$sim{"IL"}=1;
$sim{"RK"}=1;
$sim{"RH"}=1;
$sim{"KR"}=1;
$sim{"HR"}=1;
$sim{"KH"}=1;
$sim{"HK"}=1;
$sim{"DE"}=1;
$sim{"ED"}=1;
$sim{"ST"}=1;
$sim{"TS"}=1;
$sim{"NQ"}=1;
$sim{"QN"}=1;

### Reading Blosum62 matrix
#   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
#A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
my ($pair, $j);
my (@A1,@A2, @line);
my (%sc);
my $c=-1;
open (BL,$bl) || die "Missing Blossum Matrix";
while (<BL>){
    chomp($_); 
    $c++;    
    if ($c==0){
	@A1=split(/\s+/,$_);
	@A2=split(/\s+/,$_);
    }
    else {
	@line=split(/\s+/,$_);
	
	for ($j=1; $j<scalar(@line)-1; $j++){
	    $pair=$A1[$c].$A2[$j];
	    $sc{$pair}=$line[$j];
	}
    }
}
close (BL);
##### Parses Gene cluster (sub)file 
### First to get the valid species
my (%c1,%c2); 
my (@l);
open (IN,"$i1") || die "It cannot open gene cluster file"; ##Gene cluster files
while (<IN>){
    chomp($_);
    @l=split(/\t/,$_);
    if ($l[1] eq $s1) {
	$c1{$l[0]}=1;
    }
    elsif ($l[1] eq $s2) {
	$c2{$l[0]}=1;
    }
}
close IN;

my (%gn, %cl1, %cl2);
open (IN,"$i1") || die "It cannot open gene cluster file (second time)"; ##Gene cluster files
while (<IN>){
    chomp($_);
    @l=split(/\t/,$_);    
    if ($c1{$l[0]} && $c2{$l[0]}){ # only for the valid clusters
	if ($l[1] eq $s1 ) {
	    $gn{$l[2]}=$l[0];
	    if (!$cl1{$l[0]}) {
		$cl1{$l[0]}=$l[2];
	    } 
	    else { 
		$cl1{$l[0]}.=",".$l[2];
	    }
	}
	elsif ($l[1] eq $s2) {
	    $gn{$l[2]}=$l[0];
	    if (!$cl2{$l[0]}) {
		$cl2{$l[0]}=$l[2];
	    } 
	    else { 
		$cl2{$l[0]}.=",".$l[2];  
	    }
	}
    }
}
close IN;

### protein ids exons species 1
### ENSP00000000412|ENSG00000003056
open (IN,"$i2") || die "It cannot open protein ids exon file species 1\n"; 
my (%pid1, %pid2,%spid,%chpr); 
my $cid;
my (%p1, %p2);
while (<IN>){
    chomp($_);
    @l=split(/\|/,$_);   
    if ($gn{$l[1]}){
	$cid=$gn{$l[1]};
	$p1{$_}=1;
	$spid{$_}=$sp1; ##adding protid species info	
	if (!$pid1{$cid}){
	    $pid1{$cid}=$_;
	    $chpr{$_}=1;
	} 
	else { 
	    if (!$chpr{$_}){
		$pid1{$cid}.=",".$_; 
		$chpr{$_}=1;
	    }
	}
    }
}
close IN;

### protein ids exons species 2
### ENSP00000000412|ENSG00000003056
open (IN,"$i3") || die "It cannot open protein ids exon file species 2";
while (<IN>){
    chomp($_);
    @l=split(/\|/,$_);    
    if ($gn{$l[1]}){
	$cid=$gn{$l[1]};
	$p2{$_}=1;
	$spid{$_}=$sp2;##adding protid info
	if (!$pid2{$cid}){
	    $pid2{$cid}=$_;
	    $chpr{$_}=1;
	} 
	else { 
	    if (!$chpr{$_}){
		$pid2{$cid}.=",".$_; 
		$chpr{$_}=1;
	    }
	}
    }
}
close IN;

my %exon;
my $rs1;
my %pos;
my @n;
open (IN,"$i4")|| die "It cannot open exon position files species 1\n";
while (<IN>){
    chomp($_);
    @line=split(/\t/,$_);
    @l=split(/\-/,$line[2]);
    @n=split(/\|/,$line[0]);
    if ($chpr{$line[0]}){    
	if (!$pos{$line[0]}) { $pos{$line[0]}=$_; } else { $pos{$line[0]}.="\n".$_; }    
	for ($rs1=$l[0];  $rs1<=$l[1]; $rs1++){ 
	    $exon{$line[0]."_".$rs1}=$_;
	}
    } 
}
close IN;

open (IN,"$i5") || die "It cannot open exon position files species 2\n"; 
while (<IN>){
    chomp($_);
    @line=split(/\t/,$_);
    @l=split(/\-/,$line[2]);
    @n=split(/\|/,$line[0]);    
    if ($chpr{$line[0]}){
	if (!$pos{$line[0]}) { $pos{$line[0]}=$_; } else { $pos{$line[0]}.="\n".$_; }
	for ($rs1=$l[0];  $rs1<=$l[1]; $rs1++){ 
	    $exon{$line[0]."_".$rs1}=$_;
	}
    }
}
close IN;

### Intron position files for each species
### Format: BL10724_cuf1|BL10724Sc0000009+Intron_13599249-3601287
my (%intron_aa,%intron_gg);
open (INFILE,"$i6") || die "It cannot open file with intron positions for species 1\n"; ### Sp1
while (<INFILE>){
    chomp($_); 
    @line=split(/\s+/,$_);
    $line[3]=lc($line[3]); # turns it into lower case
    
    if ($chpr{$line[0]}){
	if ($line[3] eq "intronless"){
	    $intron_gg{$line[0]."\t".$line[3]}=$line[1]."\t"."NA"."\t".$line[2];
	}
	else {
	    $intron_gg{$line[0]."\t".$line[3]}=$line[1]."\t".$line[4]."\t".$line[2];
	}
    }
}
close INFILE;

open (INFILE,"$i7") || die "It cannot open file with intron positions for species 2\n"; ### Sp2
while (<INFILE>){
    chomp($_); 
    @line=split(/\s+/,$_);
    $line[3]=lc($line[3]);
    
    if ($chpr{$line[0]}){
	if ($line[3] eq "intronless"){
	    $intron_gg{$line[0]."\t".$line[3]}=$line[1]."\t"."NA"."\t".$line[2];
	}
	else {
	    $intron_gg{$line[0]."\t".$line[3]}=$line[1]."\t".$line[4]."\t".$line[2];
	}
    }
}
close INFILE;

##### PROCESSING EXINT FILES
my (@t1, @t2);
my $f1=$i8;
my $f2=$i9;
my $m=0; 
my (%seqs, %tmp_int,%intron);
my ($ic, $x,$id,$z,$dev,$sid); 
my (@t3);
open (EXONE, "$f1") || die "It cannot open exint file species 1 ($f1)";
while (<EXONE>){
    chomp($_);    
    if ($_=~/\>/){
	$m=0;
	@l=split(/\s+/,$_);
	$l[0]=~s/\>//;
	if ($p1{$l[0]}){
	    $tmp_int{$l[0]}=$_;
	    $sid=$_;
	    $m=1;
	    $x=0;    
	    for ($z=1; $z<scalar(@l); $z++){
		$x++;
		@t3=split(/\./,$l[$z]);
		$id=$l[0]."\tintron_".$x;
		$intron_aa{$id}=$t3[0];
	    }
	}
    }
    elsif($m==1) {
	$seqs{$l[0]}=$sid."\n".$_;
    }
}
close (EXONE); 
########
open (EXTWO, "$f2")|| die "It cannot open exint file species 2 ($f2)";
while (<EXTWO>){
    chomp($_);   
    if ($_=~/\>/){
	$m=0;
	@l=split(/\s+/,$_);
	$l[0]=~s/\>//;	
	if ($p2{$l[0]}){
	    $tmp_int{$l[0]}=$_;
	    $sid=$_;
	    $m=1;
	    $x=0;	        
	    for ($z=1; $z<scalar(@l); $z++){
		$x++;
		@t3=split(/\./,$l[$z]);
		$id=$l[0]."\tintron_".$x;
		$intron_aa{$id}=$t3[0];
	    }
	}
    }
    elsif($m==1) {
	$seqs{$l[0]}=$sid."\n".$_;
    }
}
close (EXTWO);

##OPTIONAL: upload pre-computed data to avoid repeating ALNs
my $pre_EX_file="$prev_folder/$s1-$s2/all_EX_aln_features.txt";
my $pre_INT_file="$prev_folder/$s1-$s2/all_INT_aln_features.txt";
my $pre_PROT_file="$prev_folder/$s1-$s2/all_PROT_aln_features.txt";
my $pre_ALN_file="$prev_folder/$s1-$s2/EXINT_aln.gz";
my %pre_EX_data;
my %pre_INT_data;
my %pre_PROT_data;
my %pre_ALN_data;

if ($prev_folder){ # option is provided
    if (-e $pre_EX_file && -e $pre_INT_file && -e $pre_PROT_file && -e $pre_ALN_file){ # all must exist
	### EXON DATA
	open (PRE_EX, $pre_EX_file) || die "It cannot open $pre_EX_file\n";
#	<PRE_EX>;
	while (<PRE_EX>){ 
	    chomp;
	    my @t=split(/\t/);
	    my $prot_sp1 = $t[1];
	    my $prot_sp2 = $t[8];
	    my ($gene_sp1) = $prot_sp1 =~/.+\|(.+)/;
	    my ($gene_sp2) = $prot_sp2 =~/.+\|(.+)/;
	    my $temp_cl_sp1; my $temp_cl_sp2;
	    $temp_cl_sp1 = $gn{$gene_sp1} if defined $gn{$gene_sp1};
	    $temp_cl_sp2 = $gn{$gene_sp2} if defined $gn{$gene_sp2};
	    
	    if ($temp_cl_sp1 && $temp_cl_sp1 eq $temp_cl_sp2){
		my $string = join("\t",$temp_cl_sp1,@t[1..$#t]);
		$pre_EX_data{$prot_sp1}{$prot_sp2}.= "$string\n";
	    } # else: the two genes no longer belong to the same cluster
	}
	close PRE_EX;
	### INTRON DATA
	open (PRE_INT, $pre_INT_file) || die "It cannot open $pre_INT_file\n";
#	<PRE_INT>;
	while (<PRE_INT>){ 
	    chomp;
	    my @t=split(/\t/);
	    my $prot_sp1 = $t[1];
	    my $prot_sp2 = $t[5];
	    my ($gene_sp1) = $prot_sp1 =~/.+\|(.+)/;
	    my ($gene_sp2) = $prot_sp2 =~/.+\|(.+)/;
	    my $temp_cl_sp1; my $temp_cl_sp2;
	    $temp_cl_sp1 = $gn{$gene_sp1} if defined $gn{$gene_sp1};
	    $temp_cl_sp2 = $gn{$gene_sp2} if defined $gn{$gene_sp2};
	    
	    if ($temp_cl_sp1 && $temp_cl_sp1 eq $temp_cl_sp2){
		my $string = join("\t",$temp_cl_sp1,@t[1..$#t]);
		$pre_INT_data{$prot_sp1}{$prot_sp2}.= "$string\n";
	    } # else: the two genes no longer belong to the same cluster
	}
	close PRE_INT;
	### PROT DATA
	open (PRE_PROT, $pre_PROT_file) || die "It cannot open $pre_PROT_file\n";
#	<PRE_PROT>;
	while (<PRE_PROT>){ 
	    chomp;
	    my @t=split(/\t/);
	    my $prot_sp1 = $t[2];
	    my $prot_sp2 = $t[3];
	    my ($gene_sp1) = $prot_sp1 =~/.+\|(.+)/;
	    my ($gene_sp2) = $prot_sp2 =~/.+\|(.+)/;
	    my $temp_cl_sp1; my $temp_cl_sp2;
	    $temp_cl_sp1 = $gn{$gene_sp1} if defined $gn{$gene_sp1};
	    $temp_cl_sp2 = $gn{$gene_sp2} if defined $gn{$gene_sp2};
	    
	    if ($temp_cl_sp1 && $temp_cl_sp1 eq $temp_cl_sp2){
		my $string = join("\t",$temp_cl_sp1,@t[1..$#t]);
		$pre_PROT_data{$prot_sp1}{$prot_sp2}.= "$string\n";
	    } # else: the two genes no longer belong to the same cluster
	}
	close PRE_PROT;
	### ALN DATA
	my ($prot_sp1_A,$prot_sp2_A, $gene_sp1_A, $gene_sp2_A, $valid_cl);
	open (PRE_ALN, "gunzip -c $pre_ALN_file |") || die "It cannot open $pre_ALN_file\n";
	while (<PRE_ALN>){ 
	    chomp;
	    if (/^\>\>\>/){
		my @t = split(/ /,$_);
		$prot_sp1_A = $t[3];
		$prot_sp2_A = $t[4];
		($gene_sp1_A) = $prot_sp1_A =~/.+\|(.+)/;
		($gene_sp2_A) = $prot_sp2_A =~/.+\|(.+)/;
		my $temp_cl_sp1; my $temp_cl_sp2;
		$temp_cl_sp1 = $gn{$gene_sp1} if defined $gn{$gene_sp1};
		$temp_cl_sp2 = $gn{$gene_sp2} if defined $gn{$gene_sp2};
		
		if ($temp_cl_sp1 && $temp_cl_sp1 eq $temp_cl_sp2){
		    $valid_cl = 1;
		    $pre_ALN_data{$prot_sp1_A}{$prot_sp2_A}.="$_\n";
		} 
		else {
		    $valid_cl = 0;
		}
	    }
	    else {
		if ($valid_cl){ # i.e. data from before
		    $pre_ALN_data{$prot_sp1_A}{$prot_sp2_A}.="$_\n";
		}
	    }
	}
	close PRE_ALN;
    }
    else {
	die "*** Some of the precomputed files for $s1-$s2 are missing in $prev_folder\n";
    }
}


##SCORING EXONS OF PROTEIN PAIRS
my @keys=keys(%pid1);
@keys=sort(@keys);
my ($n1,$n2,$e,$i,$g, $zi, $zj);
my $el;
my $tels=scalar(@keys);
my @keys2=keys(%pid2);
@keys2=sort(@keys2);
my ($k, $name);
my (%seq);
my (@s1,@s2);
my (%score);
my ($te,$tg, $tp); ##temporal files
$te=$outf."/tmp_part_".$part.".exint"; ##temporal exint file
$tg=$outf."/tmp_part_".$part.".gde"; ##temporal gde file
$tp=$outf."/pos_part_".$part.".txt"; ##temporal position
my $int_aln=$outf."/tmp_part_".$part.".INT_ALIGN.aln";
my $f_merged_aln =$outf."/$s1-$s2-part_$part.ALL.aln";
open (MERGE_ALN, ">$f_merged_aln");
foreach $el (@keys){
    @t1=(); @t2=();   
    if ($pid1{$el}=~/\,/){ @t1=split(/\,/,$pid1{$el}); }
    else { push(@t1,$pid1{$el});  }    
    if ($pid2{$el}=~/\,/){ @t2=split(/\,/,$pid2{$el}); }
    else { push(@t2,$pid2{$el});  }    
    ##getting pairwise exint file and then make alignment
#    print "$el\n"; ##printing in standard output the gene cluster ID that is being processed
    my $Gclid=$el;
    for ($zj=0; $zj<scalar(@t1); $zj++){ ##opening FOR 1
	for ($zi=0; $zi<scalar(@t2); $zi++) {  ##opening FOR 2
	    my $prot_sp1 = $t1[$zj];
	    my $prot_sp2 = $t2[$zi];
	    if ($pre_EX_data{$prot_sp1}{$prot_sp2} || $pre_EX_data{$prot_sp2}{$prot_sp1}){
		print EXSC $pre_EX_data{$prot_sp1}{$prot_sp2} if $pre_EX_data{$prot_sp1}{$prot_sp2};
		print EXSC $pre_EX_data{$prot_sp2}{$prot_sp1} if $pre_EX_data{$prot_sp2}{$prot_sp1};
		print INSC $pre_INT_data{$prot_sp1}{$prot_sp2} if $pre_INT_data{$prot_sp1}{$prot_sp2};
		print INSC $pre_INT_data{$prot_sp2}{$prot_sp1} if $pre_INT_data{$prot_sp2}{$prot_sp1};
		print PRSC $pre_PROT_data{$prot_sp1}{$prot_sp2} if $pre_PROT_data{$prot_sp1}{$prot_sp2};
		print PRSC $pre_PROT_data{$prot_sp2}{$prot_sp1} if $pre_PROT_data{$prot_sp2}{$prot_sp1};
		print MERGE_ALN $pre_ALN_data{$prot_sp1}{$prot_sp2} if $pre_ALN_data{$prot_sp1}{$prot_sp2}; # should always be
		print MERGE_ALN $pre_ALN_data{$prot_sp2}{$prot_sp1} if $pre_ALN_data{$prot_sp2}{$prot_sp1}; # should not be
	    }
	    else { # business as usual
		## 0) ADDS HEADING TO MERGE FILE
		print MERGE_ALN ">>> $s1 $s2 $t1[$zj] $t2[$zi]\n\n";
		## 1) MAKING TEMPORAL EXINT FILE
		open (TMPALN,">$te");
		print TMPALN "$seqs{$t1[$zj]}\n$seqs{$t2[$zi]}\n"; 
		close (TMPALN);
		# RUNNING ALIGN INTRON POS
		`B0_generate_IPA_prot_aln.pl $te MAFFT $cpus`; 
		# ADDS GDE TO MERGE 
		close MERGE_ALN;
		system "cat $tg >> $f_merged_aln"; # we could already print a post-processed aln instead
		system "cat $int_aln >> $f_merged_aln"; # I have no idea where this is used
		open (MERGE_ALN, ">>$f_merged_aln");
#		print "  > $t1[$zj] $t2[$zi]\n";
		## 2) OPENING OUTPUT GDE FILE
		$n1=""; $n2="";
		%seq=();
		open (ALN, "$tg");
		while (<ALN>){		
		    chomp($_);
		    if ($_=~/\%/){ 
			@l=split(/\s+/,$_);
			$name=$l[0]; 
			$name=~s/\%//;
		    }
		    else {
			$seq{$name}.=$_;
		    }
		}
		close (ALN);
		my @nkeys=keys(%seq);
		my ($seq1,$seq2);
		$n1=$nkeys[0];
		$n2=$nkeys[1];
		$seq1=$seq{$n1};
		$seq2=$seq{$n2};
		$sp1=$spid{$n1};
		$sp2=$spid{$n2};	        
		
		## 3) CALLING SUBROUTINE FOR SCORING PROTEINS
		my ($sim1,$sim2,$idt,$glt)=(0,0,0,0);
		($sim1,$sim2,$idt,$glt)=score_proteins($n1,$n2,$seq1,$seq2);
		print PRSC "$Gclid\tProtein\t$n1\t$n2\t$sim1\t$sim2\t$idt\t$glt\t$sp1\t$sp2\n"; ##printing in the protein scoring file the scores of similarity
		print PRSC "$Gclid\tProtein\t$n2\t$n1\t$sim2\t$sim1\t$idt\t$glt\t$sp2\t$sp1\n"; ##printing in the protein scoring file the scores of similarity
		
		## 4) OPENING OUTPUT INTALN FILE
		open (INTALN, "$int_aln");
		my ($ialn, $inn1, $inn2, $iseq1, $iseq2, $is1, $is2);
		my (@w1);
		my $c=-1;
		while (<INTALN>){
		    chomp($_);
		    if ($_=~/\|/){ 
			$c++; 
			if ($c==0){ $inn1=$_; }
			if ($c==1){ $inn2=$_; }
		    }
		    elsif ($_=~/\>/){ $c=-1; }
		    elsif ($_=~/\w+/) {  
			@w1=split(/\t/,$_); 
			if (scalar(@w1)==3){ 
			    $c++;
			    if ($c==0){  $is1.=$w1[1]; }
			    elsif ($c==1) { $is2.=$w1[1];  }
			} 
			elsif (scalar(@w1)!=3 && $w1[1]){
			    $ialn.=$w1[1];
			}
		    }
		    else { $c=-1; }
		} ## closing INTALN file
		
#	    if (($sim1>=20 && $sim2>=20) && ($sp1 ne $sp2) && !$score{$n1.",".$n2}){ 
		if ((($sim1>=$min_sim_prots*100 && $sim2>=$min_sim_prots*100) || $sim1 >= $min_sim_prots*100*2 || $sim2 >= $min_sim_prots*100*2) && ($sp1 ne $sp2) && !$score{$n1.",".$n2}){ 
		    ##5) CALLING SUBROUTINE FOR SCORING INTRONS
#		print PRSC "$Gclid\tProtein\t$n1\t$n2\t$sim1\t$sim2\t$id1\t$glsc1\t$sp1\t$sp2\n"; # redundant
		    if ($is1 && $is2 && $ialn && $inn1 && $inn2){ # added $inn1 and $inn2
			my $tsp1=$spid{$inn1}; my $tsp2=$spid{$inn2};
			my $tmp1=score_introns($is1,$ialn,$is2,$inn1,$inn2,$tsp1,$tsp2,$Gclid);
		    }
		    ##  6) SCORING EACH EXON PAIR
		    ### 6.1) GETTING RESIDUES ALIGNMENT
		    my $PP1=$pos{$n1}; my $PP2=$pos{$n2};
		    if ($seq1 && $seq2 && $n1 && $n2 && $PP1 && $PP2){
			my $tmp2=score_exons($seq1,$seq2,$n1,$n2,$sp1,$sp2,$glt,$sim1,$sim2,$idt,$Gclid,$PP1,$PP2); # idt not really used
		    }
		    else {} # print "3)$Gclid\t$n1\t$n2\t$seq1\t$seq2\t$PP1\t$PP2\n"; }
		    $score{$n1.",".$n2}=1;
		} 
	    } # CLOSING ELSE
	}##CLOSING FOR 2
    }##CLOSING FOR 1
} ##closing MAIN KEYS
#######SUB-ROUTINES######

# I) SCORING PROTEIN ALIGNMENTS
sub score_proteins {
    my ($n1,$n2,$seq1,$seq2)=@_;
    my (@s1,@s2);
    @s1=split //,$seq1;
    @s2=split //,$seq2;
    
    my $ng=-1;
    my $id_score=0;
    my $f_score=0;
    my $sim_score=0;
    my $e_score=0;
    my $m_score=0;
    my $score=0;
    my $global_score1=0;
    my $global_score2=0;
    my $m1;
    my $m2;
    my $res=0;
    my $g=0;
    my ($l1,$l2);
    my $t1=$seq1;
    $t1=~s/\-//g;
    my $t2=$seq2;
    $t2=~s/\-//g;
    $l1=length($t1); 
    $l2=length($t2); 
    my ($sim1,$sim2,$id1,$id2)=(0,0,0,0);
    my $gap=0;
    my $e=0;
    my $i;
    
    ##scoring proteins
    for ($i=0; $i<scalar(@s1); $i++){
	if ($sc{$s1[$i].$s2[$i]}){
	    $e_score+=$sc{$s1[$i].$s2[$i]};
	}
	if ($s1[$i] eq $s2[$i]){
	    $id_score++;
	    $sim_score++;
	    $ng=0;
	}
	elsif ($sim{$s1[$i].$s2[$i]}){
	    $sim_score++;
	    $ng=0;
	}
	elsif ($s2[$i] eq "-" || $s1[$i] eq "-") { if ($ng==0) { $g++; $ng=1; } else { $e++; $ng=1;} } 
	if ($sc{$s1[$i].$s1[$i]}) { $m1+=$sc{$s1[$i].$s1[$i]}; }
	if ($sc{$s2[$i].$s2[$i]}) { $m2+=$sc{$s2[$i].$s2[$i]}; }
    }
    
    if ($l1>0 && $l2>0){

	$sim1=sprintf("%.2f",(($sim_score/$l1)*100));
	$sim2=sprintf("%.2f",(($sim_score/$l2)*100));
    
	$id1=sprintf("%.2f",(($id_score/$l1)*100));
	$id2=sprintf("%.2f",(($id_score/$l2)*100));
        if ($m1 == 0) {
		$global_score1 = 0;
	} else {
		$global_score1=($e_score+($g*-4)+($e*-1))/$m1;
	}
	if ($m2 == 0) {
                $global_score2 = 0;
        } else {
		$global_score2=($e_score+($g*-4)+($e*-1))/$m2;
    	}
	$global_score1=sprintf("%.2f",$global_score1);
	$global_score2=sprintf("%.2f",$global_score2);
    
	 ##It will report the minimum score
	my ($simt, $idt, $glt);
	if ($sim1<$sim2) { $simt=$sim1; } else { $simt=$sim2; }
	if ($id1<$id2) { $idt=$id1; } else { $idt=$id2; }
	if ($global_score1<$global_score2) { $glt=$global_score1; } else { $glt=$global_score2; }
	
	#RETURNING RESULTS OF SIMILARITY
	return ($sim1,$sim2,$idt,$glt);
    }
    

}

# II) SCORING INTRONS
sub score_introns {
    my ($s1,$aln,$s2,$n1,$n2,$sp1,$sp2,$el)=@_;
    my @nseq1=split //,$s1;
    $aln=~s/\_/ /g;
    my @naln=split //,$aln;
    my @nseq2=split //,$s2;
    my $c=-1;
    my ($nid, $nc, $idin);
    my ($t)="";
    my %name;
    my ($c1,$c2,$c3)=(0,0,0);
    my (%pos1,%pos2,%pos3, %p1, %p2, %in2, %id1, %id2,  %apos1, %apos2, %insc, %indev)=();
    my $p=0;
    my $r=0;
    my $n;
    my ($i1,$i2)=(0,0);
    my ($p1,$p2)=(0,0);
    my (@seq1,@seq2,@aln)=();
    my ($jn, $ks, $paa1, $paa2, $pa1, $pa2, $f, $win);
    my (%in1, %h2);
    my $dev=0;
    my $gp=0;
    my $apos=0;
    my $np=0;
    my $tsimfr=0;
    my $simfr=0;
    my $ngp=0;

    for ($n=0; $n<scalar(@nseq1); $n++){ ##getting aln positions
	if (($nseq1[$n] || $nseq2[$n] || $nseq1[$n] == 0 || $nseq2[$n]==0) && !defined $naln[$n]){
	    $naln[$n]=" ";
	    push (@seq1,$nseq1[$n]);
	    push (@seq2,$nseq2[$n]);
	    push (@aln,$naln[$n]);
	}
	elsif ($nseq1[$n]=~/\s+/ && $nseq2[$n] =~/\s+/ && $naln[$n]=~/\s+/ && defined $naln[$n] && defined $nseq1[$n] && defined $nseq2[$n]){}
	elsif ($nseq1[$n] || $nseq2[$n] || $naln[$n] || $nseq1[$n] == 0 || $nseq2[$n]==0 ){
	    push (@seq1,$nseq1[$n]);
	    push (@seq2,$nseq2[$n]);
	    push (@aln,$naln[$n]);
	}
    }
    my $score=0;
    (%insc, %indev)=();
    ##getting intron scores for Species 1 vs Species 2##
    $gp=0;
    for ($n=0; $n<scalar(@seq1); $n++){ 
	$r=$n+1;
	if ($seq1[$n]=~/\d/){ $i1++; $dev=0; } if ($seq2[$n]=~/\d/){ $i2++; } 
	if ($seq1[$n]=~/\d/){
	    #$i1++;
	    if ($seq2[$n]=~/\d/){
		$dev=1;
		my $m1=$n-1; my  $m2=$n+1;
		$paa1=$intron_aa{$n1."\tintron_".$i1};
		$paa2=$intron_aa{$n2."\tintron_".$i2};
		$score=0; 
		##being in the same phase sums up 25; otherwise it has a value of -25; then if they are not in the same position then is rested -1 for each deviated position.
		##A 0 score, means same phase, but a deviation of 5 aa;
		##A negative score is a non conserved phase 
		##therefore a score >=0 means that the intron is conserved with ot without aa deviation 
		if ($paa1 && $paa2){ ##added this 
		    if ($seq1[$n] == $seq2[$n]){ ##checking if the phase is the same for both introns
			$score=10;
			if ($n1 && $n2){
			    $aln[$m1]="n" if !defined $aln[$m1]; 
			    $aln[$m2]="n" if !defined $aln[$m2];			    
			    print INSC "$el\t$n1\tintron_$i1\t$paa1\t$seq1[$n]\t$n2\tintron_$i2\t$paa2\t$seq2[$n]\t0\t$aln[$m1],$aln[$m2]\t$score\t$sp1\t$sp2\n";
			}
		    }
		    else {
			$score=-10; ##introns in different phase hava a score of -10 (change in the scoring system)
			if ($n1 && $n2){
			    $aln[$m1]="n" if !defined $aln[$m1]; 
			    $aln[$m2]="n" if !defined $aln[$m2];
			    print INSC "$el\t$n1\tintron_$i1\t$paa1\t$seq1[$n]\t$n2\tintron_$i2\t$paa2\t$seq2[$n]\t0\t$aln[$m1],$aln[$m2]\t$score\t$sp1\t$sp2\n";
			}
		    }
		}
	    }
	    else { 
		##First check alignment quality
		$pa1=$n-12; if ($pa1<0){ $pa1=0; }
		$pa2=$n+12; if ($pa2>scalar(@seq1)){ $pa2=scalar(@seq1)-1; }
		$np=0;
		$tsimfr=0;
		$simfr=0;
		$ngp=0;
		$win=10;
		for ($f=$pa1; $f<$pa2; $f++){
		    if ($seq1[$f]=~/[A-Z]/ && $seq2[$f]=~/[A-Z]/){
			$np++;	
			if ($sim{$seq1[$f].$seq2[$f]}){	
			    $simfr+=$sim{$seq1[$f].$seq2[$f]};
			}			
		    }
		    if ($seq1[$f] eq "-" || $seq2[$f] eq "-"){
			$ngp++;	
			$np++;			
		    }
		    
		}		
		if ($np>0){
		    $tsimfr=($simfr/$np)*100;
		    if ($tsimfr<30 || $ngp>=($np*0.3)){ $win=10;  }
		    elsif ($tsimfr>=30 && $tsimfr<50){ $win=8;  }
		    elsif ($tsimfr>=50 && $tsimfr<70){ $win=6;  }
		    elsif ($tsimfr>=70 && $tsimfr<80){ $win=4;  }
		    elsif ($tsimfr>=80 && $tsimfr<90){ $win=3;  }
		    elsif ($tsimfr>=90){ $win=2; }
		}
		##ALN quality
		$j=-$win; ##changing to a window of sliding depending of alignment quality
		$r=$i2;
		$score=0; $apos=0; $gp=0;
		for ($t=$n+$j; $t<=$n+$win; $t++){ ##changing the sliding window according to similarity of upstream and downstream regions
		    if ($t <= $#seq2){
			if ($seq2[$t] eq "-"){ $gp++; $apos++; }
			if ($seq2[$t]=~/[A-Z]/ ) { $apos++; }
			if ($seq2[$t]=~/\d/){		    
			    $dev=1;
			    if ($j>0) { $r++; }
			    my $m1=$t-1; my $m2=$t+1; my $m3=$j-1;
			    $idin=$el."\t".$n1."\tintron_".$i1;
			    $paa1=$intron_aa{$n1."\tintron_".$i1};
			    $paa2=$intron_aa{$n2."\tintron_".$r};
			    $score=0;
			    if ($paa1 && $paa2){ ##added this 
				if ($m3>0){			    
				    if ($seq1[$n] eq $seq2[$t]) { $score=10; } else { $score=-10; } ##change in the scoring system
				    $score=$score-(abs($m3)); ##resting deviating positions to the score
				    ###differentiating from NO_ALN
				    if ($score==0){ $score=-1; }
				    if (!$insc{$idin}){
					if ($n1 && $n2){
					    $aln[$m1]="n" if !defined $aln[$m1]; 
					    $aln[$m2]="n" if !defined $aln[$m2];
					    $insc{$idin}=$el."\t".$n1."\tintron_".$i1."\t".$paa1."\t".$seq1[$n]."\t".$n2."\tintron_".$r."\t".$paa2."\t".$seq2[$t]."\t+".$m3."\t".$aln[$m1].",".$aln[$m2]."\t".$score."\t".$sp1."\t".$sp2;
					    
					    $indev{$idin}=abs($m3);
					}
					#print INSC "$el\t$n1\tintron_$i1\t$paa1\t$seq1[$n]\t$n2\tintron_$r\t$paa2\t$seq2[$t]\t+$m3\t$aln[$m1],$aln[$m2]\t$score\t$sp1\t$sp2\n";
				    }
				    elsif (abs($m3)<$indev{$idin}) {
					if ($n1 && $n2){
					    $aln[$m1]="n" if !defined $aln[$m1]; 
					    $aln[$m2]="n" if !defined $aln[$m2];					    
					    $insc{$idin}=$el."\t".$n1."\tintron_".$i1."\t".$paa1."\t".$seq1[$n]."\t".$n2."\tintron_".$r."\t".$paa2."\t".$seq2[$t]."\t+".$m3."\t".$aln[$m1].",".$aln[$m2]."\t".$score."\t".$sp1."\t".$sp2;
					    
					    $indev{$idin}=abs($m3);
					}
				    }
				}
				else {
				    my $w=$m3+2; 
				    if ($seq1[$n] eq $seq2[$t]) { $score=10; } else { $score=-10; }
				    $score=$score-(abs($w)); ##resting deviating positions to the score
				    ###differentiating from NO_ALN
				    if ($score==0){ $score=-1; }
				    if (!$insc{$idin}){
					
					if ($n1 && $n2){
					    $aln[$m1]="n" if !defined $aln[$m1]; 
					    $aln[$m2]="n" if !defined $aln[$m2];
					    $insc{$idin}=$el."\t".$n1."\tintron_".$i1."\t".$paa1."\t".$seq1[$n]."\t".$n2."\tintron_".$r."\t".$paa2."\t".$seq2[$t]."\t".$w."\t".$aln[$m1].",".$aln[$m2]."\t".$score."\t".$sp1."\t".$sp2;
					    $indev{$idin}=abs($w);
					}
					#print INSC "$el\t$n1\tintron_$i1\t$paa1\t$seq1[$n]\t$n2\tintron_$r\t$paa2\t$seq2[$t]\t+$m3\t$aln[$m1],$aln[$m2]\t$score\t$sp1\t$sp2\n";
				    }
				    elsif (abs($w)<$indev{$idin}) {
					if ($n1 && $n2){
					    $aln[$m1]="n" if !defined $aln[$m1]; 
					    $aln[$m2]="n" if !defined $aln[$m2];
					    $insc{$idin}=$el."\t".$n1."\tintron_".$i1."\t".$paa1."\t".$seq1[$n]."\t".$n2."\tintron_".$r."\t".$paa2."\t".$seq2[$t]."\t".$w."\t".$aln[$m1].",".$aln[$m2]."\t".$score."\t".$sp1."\t".$sp2;
					    $indev{$idin}=abs($w);
					}
				    }
				    #print INSC "$el\t$n1\tintron_$i1\t$paa1\t$seq1[$n]\t$n2\tintron_$r\t$paa2\t$seq2[$t]\t$w\t$aln[$m1],$aln[$m2]\t$score\t$sp1\t$sp2\n";
				}
			    }
			}
		    } #for paa1 & paa2 conditional
		    $j++;
		}
		##Printing only the pair of introns with less deviation and calculating the corresponding score
		if ($dev && $insc{$idin}){ print INSC "$insc{$idin}\n"; }
	    }	    
	    if ($dev==0 || !$paa2){ ##no intron aligned
		$paa1=$intron_aa{$n1."\tintron_".$i1};
		if ($n1 && $n2){
		    if ($gp>=int($apos*0.6)){ ##bad alignment in the protein query
			print INSC "$el\t$n1\tintron_$i1\t$paa1\t$seq1[$n]\t$n2\tNO_ALN\tNA\tNA\tNA\tNA\tNA\t$sp1\t$sp2\n";
			
		    }else {
		   	print INSC "$el\t$n1\tintron_$i1\t$paa1\t$seq1[$n]\t$n2\tNO_INTRON\tNA\tNA\tNA\tNA\t0\t$sp1\t$sp2\n";
		    }
		}
	    }
	    
	}
    }
    ###END OF PROCESSING THE ALIGNMENT SPECIES 1 VS 2     
    ##getting intron scores for Species 2 vs Species 1##           
    ($i1,$i2)=(0,0); 
    $dev=0;
    $gp=0;
    (%insc, %indev)=();
    for ($n=0; $n<scalar(@seq2); $n++){ ##getting aln positions
	$r=$n+1;
	if ($seq2[$n]=~/\d/){ $i1++; $dev=0; } if ($seq1[$n]=~/\d/){ $i2++; }
	if ($seq2[$n]=~/\d/){
	    #$i1++;
	    if ($seq1[$n]=~/\d/){
		$dev=1;
		my $m1=$n-1; my  $m2=$n+1;		
		$paa1=$intron_aa{$n2."\tintron_".$i1};
		$paa2=$intron_aa{$n1."\tintron_".$i2};
		$score=0; 
		##being in the same phase sums up 5; otherwise it rests -1; then if they are not in the same position then is rested -1 for each deviated position.
		##A 0 score, means same phase, but a deviation of 5 aa;
		##A negative score is a non conserved phase 
		##therefore a score >=0 means that the intron is conserved with ot without aa deviation 
		if ($paa1 && $paa2){
		    if ($seq2[$n] == $seq1[$n]){ ##checking if the phase is the same for both introns
			$score=10;
			if ($n1 && $n2){
			    $aln[$m1]="n" if !defined $aln[$m1]; 
			    $aln[$m2]="n" if !defined $aln[$m2];
			    print INSC "$el\t$n2\tintron_$i1\t$paa1\t$seq2[$n]\t$n1\tintron_$i2\t$paa2\t$seq1[$n]\t0\t$aln[$m1],$aln[$m2]\t$score\t$sp2\t$sp1\n";
			}
		    }
		    else {
			$score=-10;
			if ($n1 && $n2){
			    $aln[$m1]="n" if !defined $aln[$m1]; 
			    $aln[$m2]="n" if !defined $aln[$m2];
			    print INSC "$el\t$n2\tintron_$i1\t$paa1\t$seq2[$n]\t$n1\tintron_$i2\t$paa2\t$seq1[$n]\t0\t$aln[$m1],$aln[$m2]\t$score\t$sp2\t$sp1\n";
			}
		    }
		}
	    }
	    else {
		##First check alignment quality
		$pa1=$n-12; if ($pa1<0){ $pa1=0; }
		$pa2=$n+12; if ($pa2>scalar(@seq1)){ $pa2=scalar(@seq1)-1; }
		$np=0;
		$tsimfr=0;
		$simfr=0;
		$ngp=0;
		$win=10;
		for ($f=$pa1; $f<$pa2; $f++){
		    if ($seq1[$f]=~/[A-Z]/ && $seq2[$f]=~/[A-Z]/){
			$np++;	
			if ($sim{$seq1[$f].$seq2[$f]}){	
			    $simfr+=$sim{$seq1[$f].$seq2[$f]};
			}			
		    }
		    if ($seq1[$f] eq "-" || $seq2[$f] eq "-"){
			$ngp++;	
			$np++;			
		    }
		    
		}
		if ($np>0){
		    $tsimfr=($simfr/$np)*100;
		    if ($tsimfr<30 || $ngp>=($np*0.3)){ $win=10;  }
		    elsif ($tsimfr>=30 && $tsimfr<50){ $win=8;  }
		    elsif ($tsimfr>=50 && $tsimfr<70){ $win=6;  }
		    elsif ($tsimfr>=70 && $tsimfr<80){ $win=4;  }
		    elsif ($tsimfr>=80 && $tsimfr<90){ $win=3;  }
		    elsif ($tsimfr>=90){ $win=2; }
		}
		##ALN quality
		$j=-$win; ##changing deviation accordint to alignment quality
		my $r=$i2; $apos=0;
		$score=0; $gp=0;
		for ($t=$n+$j; $t<=$n+$win; $t++){
		    if ($t <= $#seq1){ # in some cases, $t is bigger than $#seq1
			if ($seq1[$t] eq "-"){ $gp++; $apos++; } # @seq1 is the positions of the aln
			if ($seq1[$t] =~/[A-Z]/) { $apos++; }
			if ($seq1[$t]=~/\d/){ 		    
			    $dev=1;
			    if ($j>0) { $r++; }
			    my $m1=$t-1; my $m2=$t+1; my $m3=$j-1;
			    $idin=$el."\t".$n2."\tintron_".$i1;
			    $paa1=$intron_aa{$n2."\tintron_".$i1};
			    $paa2=$intron_aa{$n1."\tintron_".$r};
			    $score=0;
			    if ($paa1 && $paa2){
				if ($m3>0){
				    if ($seq2[$n] eq $seq1[$t]) { $score=10; } else { $score=-10; }
				    $score=$score-(abs($m3)); ##resting deviating positions to the score
				    if (!$insc{$idin}){
					if ($n1 && $n2){
					    $aln[$m1]="n" if !defined $aln[$m1]; 
					    $aln[$m2]="n" if !defined $aln[$m2];
					    $insc{$idin}=$el."\t".$n2."\tintron_".$i1."\t".$paa1."\t".$seq2[$n]."\t".$n1."\tintron_".$r."\t".$paa2."\t".$seq1[$t]."\t+".$m3."\t".$aln[$m1].",".$aln[$m2]."\t".$score."\t".$sp2."\t".$sp1;
					    $indev{$idin}=abs($m3);
					}
				    }
				    elsif (abs($m3)<$indev{$idin}){
					if ($n1 && $n2){
					    $aln[$m1]="n" if !defined $aln[$m1]; 
					    $aln[$m2]="n" if !defined $aln[$m2];
					    $insc{$idin}=$el."\t".$n2."\tintron_".$i1."\t".$paa1."\t".$seq2[$n]."\t".$n1."\tintron_".$r."\t".$paa2."\t".$seq1[$t]."\t+".$m3."\t".$aln[$m1].",".$aln[$m2]."\t".$score."\t".$sp2."\t".$sp1;
					    $indev{$idin}=abs($m3);
					}
				    }
				    #print INSC "$el\t$n2\tintron_$i1\t$paa1\t$seq2[$n]\t$n1\tintron_$r\t$paa2\t$seq1[$t]\t+$m3\t$aln[$m1],$aln[$m2]\t$score\t$sp2\t$sp1\n";
				}
				else {
				    my $w=$m3+2; 
				    if ($seq2[$n] eq $seq1[$t]) { $score=10; } else { $score=-10; }
				    $score=$score-(abs($w)); ##resting deviating positions to the score
				    
				    if (!$insc{$idin}){
					if ($n1 && $n2){
					    $aln[$m1]="n" if !defined $aln[$m1]; 
					    $aln[$m2]="n" if !defined $aln[$m2];
					    $insc{$idin}=$el."\t".$n2."\tintron_".$i1."\t".$paa1."\t".$seq2[$n]."\t".$n1."\tintron_".$r."\t".$paa2."\t".$seq1[$t]."\t".$w."\t".$aln[$m1].",".$aln[$m2]."\t".$score."\t".$sp2."\t".$sp1;
					    $indev{$idin}=abs($w);
					}
				    }
				    elsif (abs($w)<$indev{$idin}){
					if ($n1 && $n2){
					    $aln[$m1]="n" if !defined $aln[$m1]; 
					    $aln[$m2]="n" if !defined $aln[$m2];
					    $insc{$idin}=$el."\t".$n2."\tintron_".$i1."\t".$paa1."\t".$seq2[$n]."\t".$n1."\tintron_".$r."\t".$paa2."\t".$seq1[$t]."\t".$w."\t".$aln[$m1].",".$aln[$m2]."\t".$score."\t".$sp2."\t".$sp1;
					    $indev{$idin}=abs($w);
					}
				    }
				}
			    }
			}
		    }##for paa1&&paa2 cond
		    $j++;
		}
		##Printing only the pair of introns with less deviation and calculating the corresponding score
		if ($dev && $insc{$idin}){ print INSC "$insc{$idin}\n"; }
	    }
	    if ($dev==0 || !$paa2) {  ##no intron aligned
		$paa1=$intron_aa{$n2."\tintron_".$i1};
		if ($n1 && $n2){
		    if ($gp>=int($apos*0.6)){ ##bad alignment in the protein query, giving a score of 0
			print INSC "$el\t$n2\tintron_$i1\t$paa1\t$seq2[$n]\t$n1\tNO_ALN\tNA\tNA\tNA\tNA\tNA\t$sp2\t$sp1\n"; # bug corrected
			
		    }else {   		   
			print INSC "$el\t$n2\tintron_$i1\t$paa1\t$seq2[$n]\t$n1\tNO_INTRON\tNA\tNA\tNA\tNA\t0\t$sp2\t$sp1\n";  
		    } 
		}
	    }
	}
    }
###END OF PROCESSING THE ALIGNMENT SPECIES 2 VS 1
} ##End of subroutine score_introns

###III) SCORING EXONS 

sub score_exons {
    my ($seq1,$seq2,$n1,$n2,$sp1,$sp2,$glsc1,$sim1,$sim2,$id1,$el,$PP1,$PP2)=@_;    
    @s1=split(//,$seq1);
    @s2=split(//,$seq2);
    open (POS, ">$tp");
    print POS "$PP1\n$PP2\n";
    close (POS);
    $sp1=$spid{$n1}; $sp2=$spid{$n2};    
    ####SCORING EXONS OF EACH PROTEIN PAIR
    my ($res1,$res2,$res)=(0,0);
    my ($k1,$k2, $s2)=(0,0,0);
    my $r2;
    my %Iscore_exons=();
    my %Sscore_exons=();
    my %Fscore_exons=();
    my %s1_scores_sim=();
    my %s1_scores_id=();
    my %s2_scores_sim=();
    my %s2_scores_id=();
    my %s1_s2_res=();
    my %s2_s1_res=();
    my $rs2;
    my $r;
    my $png;
    my @l1;
    
    for ($i=0; $i<scalar(@s1); $i++){
	if ($s2[$i] ne "-") { $res2++; $k2=$res2; }
	if ($s1[$i] ne "-"){ $res++; $res1++; }
	$s1_s2_res{$res1}=$res2; $s2_s1_res{$res2}=$res1;	    
	##similarity score
	if ($sim{$s1[$i].$s2[$i]}){  $s1_scores_sim{$res1}=$sim{$s1[$i].$s2[$i]}; $s2_scores_sim{$res2}=$sim{$s1[$i].$s2[$i]}; }	    
	##identity score
	if ($s1[$i] eq $s2[$i]) {  $s1_scores_sim{$res1}=1; $s2_scores_sim{$res2}=1; $s1_scores_id{$res1}=1; $s2_scores_id{$res2}=1; }
	if ($s1[$i]=~/[A-Z]/ && $s2[$i] eq "-"){ $s1_scores_sim{$res1}="gap"; $s1_scores_id{$res1}="gap"; }
	if ($s1[$i]eq "-" && $s2[$i]=~/[A-Z]/ ){ $s2_scores_sim{$res2}="gap"; $s2_scores_id{$res2}="gap"; }
    }    
    my $fres1=$res1;
    my $fres2=$res2;
    
    ##Opening temporal position file
    my (%check,%print, %cin);
    my ($tex, $ex, $in1, $in2, $pr, $size, $ne, $x, $tstr);
    my (@tmp, @nex);
    my $id_score=0;
    my $f_score=0;
    my $sim_score=0;
    my $e_score=0;
    my $m_score=0;
    my $score=0;
    my $ng=0; 
    my $g=0;
    my $e=0;
    my $nr1=-1;    
    if ($fres1!=0 && $fres2!=0){ ##first if
	open (POS, "$tp");
	while (<POS>){
	    chomp($_);
	    @line=split(/\s+/,$_);
	    ##checking N1
	    if ($line[0] eq $n1){
		$id_score=0;
		$f_score=0;
		$sim_score=0;
		$e_score=0;
		$m_score=0;
		$score=0;
		$ng=0; 
		$g=0;
		$e=0;
		$pr="$n1";		    
		if (!$print{$pr}){
#		    print PRSC "$el\tProtein\t$n1\t$n2\t$sim1\t$sim2\t$id1\t$glsc1\t$sp1\t$sp2\n"; ##Adding species info in the last 2 columns => it was redundant
		    %cin=();
		    $print{$pr}=1;
		}
		@l=split(/\-/,$line[2]);		    
		if ($l[1] && $l[0]) { $size=($l[1]-$l[0])+1; } else { $size=1000000; }		    
		$nr1=-1;
		for ($r=$l[0]; $r<=$l[1]; $r++){
		    if ($s1_scores_sim{$r}){ 
			if ($s1_scores_sim{$r} ne "gap"){
			    $sim_score++;
			    $nr1=$r;
			}
			elsif ($s1_scores_sim{$r} eq "gap") { 
			    $ng++;  $nr1=$r; 
			}
		    }
		    if ($s1_scores_id{$r}){
			if ($s1_scores_id{$r} ne "gap" ){
			    $id_score++; 
			}
		    }
		    if ($s1_s2_res{$r}){ $rs2=$s1_s2_res{$r}; }
		}
		$rs1=$s1_s2_res{$l[0]}; 
		if ($size!=0){
		    $id_score=sprintf("%.2f",(($id_score/$size)*100));
		    $sim_score=sprintf("%.2f",(($sim_score/$size)*100));
		    $png=sprintf("%.2f",(($ng/$size)*100));
		} 
		else { 
		    $id_score=-1; $sim_score=-1; $png=-1;  ###when there was nothing aligned to the given exon
		}
		$tex="";
		%check=();
		if ($sim_score>=0){ ##modified to print all exons that aligned even when sim_score is low, they will be filtered in other programs
		    if (defined $rs1 && defined $rs2){
			for ($k=$rs1; $k<=$rs2; $k++){
			    if ($exon{$n2."_".$k}){
				$ex=$exon{$n2."_".$k};
				if (!$check{$ex}){
				    $check{$ex}=1;				    
				    if (!$tex){
					$tex=$ex;
				    } 
				    else { 
					$tex.=",".$ex; 
				    }
				}
			    }
			}
			if ($tex){
			    if ($rs1==0){ $rs1=1; }
			    if ($tex=~/\,/){ ##if the exon need to be realign is written in a special file for further processing
				@nex=split(/\,/,$tex); $ne=scalar(@nex); ##changing output format
				for ($x=0; $x<scalar(@nex); $x++){
				    $tstr=$el."\t".$_."\t".$ne."\t".$n2."\t".$rs1."-".$rs2."\t".$id_score."\t".$sim_score."\t".$ng."\t".$png."\t".$nex[$x]."\t".$sp1."\t".$sp2;
				    $miss{$tstr}=1;
				    #print MISS "$el\t$_\t$ne\t$n2\t$rs1-$rs2\t$id_score\t$sim_score\t$ng\t$png\t$nex[$x]\t$sp1\t$sp2\n";
				    print EXSC "$el\t$_\t$ne\t$n2\t$rs1-$rs2\t$id_score\t$sim_score\t$ng\t$png\t$nex[$x]\t$sp1\t$sp2\n";
				}
			    }
			    else {
				print EXSC "$el\t$_\t1\t$n2\t$rs1-$rs2\t$id_score\t$sim_score\t$ng\t$png\t$tex\t$sp1\t$sp2\n";
				@ln=split(/\t/,$_);
				@tn1=split(/\|/,$ln[0]);
				@tn2=split(/\|/,$n2);
				$idex=$el."\t".$tn1[1]."\t".$ln[4]."\t".$tn2[1];
				$onehit{$idex}=1;
			    }
			}
		    }
		    else { 
			print EXSC "$el\t$_\t0\t$n2\tNO_EXON_ALN\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$sp1\t$sp2\n"; 
		    }
		} 
		else { 
		    print EXSC "$el\t$_\t0\t$n2\tNO_EXON_ALN\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$sp1\t$sp2\n"; 
		}
	    }	    
	    ##checking N2
	    elsif ($line[0] eq $n2){		    
		$id_score=0;
		$f_score=0;
		$sim_score=0;
		$e_score=0;
		$m_score=0;
		$score=0;
		$ng=0; 
		$g=0;
		$e=0;		    
		$pr="$n2";		    
		if (!$print{$pr}){
		    %cin=();
#		    print PRSC  "$el\tProtein\t$n2\t$n1\t$sim2\t$sim1\t$id1\t$glsc1\t$sp2\t$sp1\n"; ##Adding species info in the last 2 columns => corrected sim1 and sim2, but also eliminated as redundant
		    $print{$pr}=1;
		}		    
		@l=split(/\-/,$line[2]);		    
		if ($l[1] && $l[0]) { $size=($l[1]-$l[0])+1; } 
		else { $size=1000000; }		    
		$nr1=-1;
		for ($r=$l[0]; $r<=$l[1]; $r++){
		    if ($s2_scores_sim{$r}){
			if ($s2_scores_sim{$r} ne "gap"){
			    $sim_score++; 
			    $nr1=$r;
			}
			elsif ($s2_scores_sim{$r} eq "gap"){ 
			    $ng++; $nr1=$r; 
			}
		    }
		    if ($s2_scores_id{$r}){
			if ($s2_scores_id{$r} ne "gap"){
			    $id_score++;
			}
		    }
		    if ($s2_s1_res{$r}){ 
			$rs2=$s2_s1_res{$r}; 
		    }
		}
		$rs1=$s2_s1_res{$l[0]};		    
		if ($size!=0){
		    $id_score=sprintf("%.2f",(($id_score/$size)*100));
		    $sim_score=sprintf("%.2f",(($sim_score/$size)*100));
		    $png=sprintf("%.2f",(($ng/$size)*100));
		} 
		else { 
		    $id_score=-1; $sim_score=-1; $png=-1; ###when there was nothing aligned to the given exon
		}		    
		$tex="";
		%check=();
		if ($sim_score>=0){ ###modified to print all exons that aligned event when sim_score is low, they will be filtered in other programs
		    if (defined $rs1 && defined $rs2){
			for ($k=$rs1; $k<=$rs2; $k++){
			    if ($exon{$n1."_".$k}){ # $n1 => prot|gene
				$ex=$exon{$n1."_".$k};
				if (!$check{$ex}){
				    $check{$ex}=1;
				    if (!$tex){
					$tex=$ex;
				    } 
				    else { 
					$tex.=",".$ex; 
				    }
				}
			    }
			}
			if ($rs1==0){ $rs1=1; }

			if ($tex=~/\,/){ ##If the exon needs to be realigned is written in an special file for its further processing
			    @nex=split(/\,/,$tex); $ne=scalar(@nex); ##changing output format
			    for ($x=0; $x<scalar(@nex); $x++){
				$tstr=$el."\t".$_."\t".$ne."\t".$n1."\t".$rs1."-".$rs2."\t".$id_score."\t".$sim_score."\t".$ng."\t".$png."\t".$nex[$x]."\t".$sp2."\t".$sp1;
				$miss{$tstr}=1;
				print EXSC "$el\t$_\t$ne\t$n1\t$rs1-$rs2\t$id_score\t$sim_score\t$ng\t$png\t$nex[$x]\t$sp2\t$sp1\n";
			    }
			} 
			else {
			    print EXSC "$el\t$_\t1\t$n1\t$rs1-$rs2\t$id_score\t$sim_score\t$ng\t$png\t$tex\t$sp2\t$sp1\n";
			    @ln=split(/\t/,$_);
			    @tn1=split(/\|/,$ln[0]);
			    @tn2=split(/\|/,$n2);
			    $idex=$el."\t".$tn1[1]."\t".$ln[4]."\t".$tn2[1];
			    $onehit{$idex}=1;
			}
		    }
		    else {
			print EXSC "$el\t$_\t0\t$n1\tNO_EXON_ALN\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$sp2\t$sp1\n"; # 1 nt exons at the Cterm
		    }
		} 
		else { 
		    print EXSC "$el\t$_\t0\t$n1\tNO_EXON_ALN\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$sp2\t$sp1\n"; 
		}
	    } 
	}##WHILE POS
    }#CLOSING IF
} ##End of subroutine score_exons
##printing file of exons to realign##
my @ks=sort(keys(%miss));
my $mex;
foreach $mex(@ks){
    @ln=split(/\t/,$mex);
    @tn1=split(/\|/,$ln[1]);
    @tn2=split(/\|/,$ln[8]);
    $idex=$ln[0]."\t".$tn1[1]."\t".$ln[5]."\t".$tn2[1];
    if (!$onehit{$idex}){
	print MISS "$mex\n";
    }
}
close (PRSC);
close (EXSC);
close (INSC);
close (MISS);
