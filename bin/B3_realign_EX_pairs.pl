#!/usr/bin/env perl
use warnings;
use strict;
my $s1=$ARGV[0]; # species 1
my $s2=$ARGV[1]; # species 2
my $i1=$ARGV[2]; # Exons to realign file ( infile = exons_to_realign.getFileName())
my $e1=$ARGV[3]; # exint sp1
my $e2=$ARGV[4]; # exint sp2
my $part=$ARGV[5]; # ${cls_parts[1]}
my $outf=$ARGV[6]; # realigned_file = "realigned_" + exons_to_realign.getBaseName() + "_${cls_parts[1]}.txt"
my $odir=$ARGV[7]; # ${sp1}_${sp2}
my $bl=$ARGV[8]; #path to blossum 62 matrix  
my $cpus=$ARGV[9]; #cpus
$cpus=1 if !$cpus;

my  %sim; ##AA similarity
mkdir $odir; #Create output directory 
$odir.="/"; #Make sure output directory ends with a slash
my $texf=$odir."tex_part_".$part.".exint"; #define temporary output exint file (for exon sequences)
my $tmpgde=$odir."tex_part_".$part.".gde"; #define temporary alignment files
print "$texf\n$tmpgde\n"; #print the filenames in the log

#set similarity scores
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

#Reading Blosum62 matrix
#   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
#A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
my ($pair, $j);
my (@A1,@A2, @line, @l);
my (%sc);
my $c=-1;
my $cdir=`pwd`; chomp($cdir);
my ($dir,$project_dir)=$outf=~/(.+)\/(.+?)\/.+?\//;
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
	    $sc{$pair}=$line[$j]; #score for each pair of AA
	}
    }
}
close (BL);

##PROCESSING EXINT FILES
my (@t1, @t2);
my $sid; 
my $m=0; 
my (%seqs); #create an hash for the sequences
my $z;
my %intron; #create an hash for the introns
my $ic; 
my $dev;
open (EXONE, "$e1") || die "Missing exint file species 1"; #Open exint file for species1
while (<EXONE>){
    chomp($_); #remove the newline character at the end of the line
    if ($_=~/\>/){ #if the line contains the fasta header
	@l=split(/\s+/,$_); #split the header on the white spaces: first element = protein/geneID; other elements = exon positions (if present) 
	$l[0]=~s/\>//; #remove the > from the protein/geneID
	$sid=$_; #set the current ID to the whole first line
    }
    else { #if the line is not header, but sequence
	$seqs{$l[0]}.=$_; #associate the sequence to the header in the hash
    }    
}
close (EXONE); #cloes the first exint

open (EXTWO, "$e2");#|| die "Missing exint file species 2"; #Open exint file for species2
while (<EXTWO>){ 
    chomp($_);
    if ($_=~/\>/){
	@l=split(/\s+/,$_);
	$l[0]=~s/\>//;
	$sid=$_;
    }
    else {	
	$seqs{$l[0]}.=$_;	
    }    
}
close (EXTWO);

##PROCESSING EXON SCORING FILES
#Exon	ENSP00000349874|ENSG00000166321	exon_6	199-235	chr10	74885156-74885267	+	ENSMUSP00000074263|ENSMUSG00000021809	199-233	59.46%	59.46%	25.41%	ENSMUSP00000074263|ENSMUSG0000
my ($n1,$n2,$sq1, $sq2, $size1, $size2, $st,$idex, $tmpex, $exonid, $sp1, $sp2, $png);
my (@cr1,@cr2);
my (%score,%ex, %seq, %rex);
my ($name);
my %pair_already_aln=();
open (IN, "$i1")|| die "Missing exons to realign file"; #opens the exons to realign file
<IN>; # removes header
while (<IN>){ 
    chomp($_); #remove newline character at the end of the line
    @l=split(/\t/,$_); #split line into array using tab as separator
    @cr1=split(/\-/,$l[3]); #split the fourth element (i.e. the AA positions of the query exon)
    $size1=$cr1[1]-$cr1[0]+1; #compute the size of the exon in AA
    $sq1=substr($seqs{$l[1]},($cr1[0]-1),$size1); #get the seq of the query exon from the hash previously created using the protein/geneID as key
    $m=16; #set m to 16 to generate indexes in a few lines
    $sp1=$l[20]; $sp2=$l[21]; #select the query and target species from the header
    @cr2=split(/\-/,$l[16]); #split the AA positions of the target exon
    $size2=$cr2[1]-$cr2[0]+1; #compute the lenght of the target exon
    $sq2=substr($seqs{$l[8]},($cr2[0]-1),$size2); ## get the seq of the target exon from the hash previously created using the protein/geneID as key
    $st=substr($l[19],0,1); #get the strand
    $idex=$l[8]."\t".$l[$m-1]."\t".$l[$m]."\t".$l[$m+1]."\t".$l[$m+2]."\t".$st; #print a tab separated list of: target_protein, target_exon_number, target_exon_position, chr_target, target_exon_coords 

    ### Adds a check to save realn the same exons multiple times
    my $pair_of_seqs = "$sq1=$sq2";

    if (!$pair_already_aln{$pair_of_seqs}){ #if the two exons have not been realigned already:
	##aligning exons locally and getting their sim scores
	open (TMPALN,">$texf"); #open the temporary exint file
	print TMPALN ">$l[1]\n$sq1\n>$l[8]\n$sq2\n"; #print the two exon sequences to compare, each with its fasta header
	close (TMPALN); #close the temporary exint file
	`B0_generate_IPA_prot_aln.pl $texf MAFFT $cpus`; #align the two exons; the output filename is built from the input filename
	## 2) OPENING OUTPUT GDE FILE
	$name="";
	%seq=();
	my (@ls);
	open (ALN, "$tmpgde"); #open the temporary alignment file
	while (<ALN>){
	    chomp($_); #remove newline character at the end of the line
	    if ($_=~/\%/){ #if header (the header of the alignment output starts with %)
		@ls=split(/\s+/,$_); #split the header on the whitespaces
		$name=$ls[0]; #get the first element of the array (proteinID)
		$name=~s/\%//; #remove the %from the proteinID
	    }
	    else { #if sequence
		$seq{$name}.=uc($_); #add the sequence to the hash with the corresponding ID, after converting it to uppercase
	    }
	}
	close (ALN); #close the alignment file
	my @nkeys=keys(%seq); #get an array with all the keys in the seq hash (meaning the two protein/geneIDs)
	my ($seq1,$seq2)=""; #initiliaze the sequences	
	($n1,$n2)=""; 
	$n1=$nkeys[0]; #query protein/geneID
	$n2=$nkeys[1]; #target protein/geneID
	$seq1=$seq{$n1}; #query exon aligned sequence (with gaps and everything)
	$seq2=$seq{$n2}; #target exon aligned sequence (with gaps and everything)
	
	## 3) CALLING SUBROUTINE FOR SCORING PROTEINS	
	if (($n1 && $n2) && ($seq1 && $seq2)){ #if we have proteinID and sequences for both the exons
	    my ($sim,$id,$gp,$glsc)=score_proteins($n1,$n2,$seq1,$seq2); #use the score proteins subroutine to get the sequence similarity, the identity, number of gaps, global score.
	    @{$pair_already_aln{$pair_of_seqs}}=($sim,$id,$gp,$glsc); ## info store for the future
	    $png=sprintf("%.2f",(($gp/$size1)*100)); #compute the percentage of gaps
	    if (!$score{$l[1]."_".$l[2]."_".$l[8]}){ ## if there is no score for that exon yet 
		$score{$l[1]."_".$l[2]."_".$l[8]}=$sim;
		$ex{$l[1]."_".$l[2]."_".$l[8]}=$l[8]."\t".$l[$m]."\t".$sim."\t".$id."\t".$gp."\t".$png."\t".$idex;
	    }
	    elsif($sim>$score{$l[1]."_".$l[2]."_".$l[8]}){ # if there is score, but the new one is higher
		$score{$l[1]."_".$l[2]."_".$l[8]}=$sim;
		$ex{$l[1]."_".$l[2]."_".$l[8]}=$l[8]."\t".$l[$m]."\t".$sim."\t".$id."\t".$gp."\t".$png."\t".$idex;
	    }	    
	}
    }
    else { # i.e. the pair of sequences have already been aligned
	my ($sim,$id,$gp,$glsc)=@{$pair_already_aln{$pair_of_seqs}};
	$png=sprintf("%.2f",(($gp/$size1)*100));
	if (!$score{$l[1]."_".$l[2]."_".$l[8]}){ # if there is no score for that exon yet
	    $score{$l[1]."_".$l[2]."_".$l[8]}=$sim;
	    $ex{$l[1]."_".$l[2]."_".$l[8]}=$l[8]."\t".$l[$m]."\t".$sim."\t".$id."\t".$gp."\t".$png."\t".$idex;
	}
	elsif($sim>$score{$l[1]."_".$l[2]."_".$l[8]}){ # if there is score, but the new one is higher
	    $score{$l[1]."_".$l[2]."_".$l[8]}=$sim;
	    $ex{$l[1]."_".$l[2]."_".$l[8]}=$l[8]."\t".$l[$m]."\t".$sim."\t".$id."\t".$gp."\t".$png."\t".$idex;
	}	    	
    }

    ### Final storage of data. If an exon with a higher match, it will re-write the entry. Otherwise, it will write the previous entry again
    $tmpex=$l[1]."_".$l[2]."_".$l[8]; # ProtID1|GeneID1   exon_N   ProtID2|GeneID2
    if ($score{$tmpex}){
	if ($score{$tmpex}>=5){ # originally it was 20
	    $rex{$tmpex}=$l[0]."\t".$l[1]."\t".$l[2]."\t".$l[3]."\t".$l[4]."\t".$l[5]."\t".$l[6]."\t"."1\t".$ex{$tmpex}."\t".$sp1."\t".$sp2;
	} else { 
	    $rex{$tmpex}=$l[0]."\t".$l[1]."\t".$l[2]."\t".$l[3]."\t".$l[4]."\t".$l[5]."\t".$l[6]."\t0\t$l[8]\tNO_EXON_ALN\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$sp1\t$sp2";
	}
    }
}   

open (OUT, ">$outf") || die "It cannot open output file ($outf)\n";
#printing header
print OUT "CID\tProt_query\tExon_number_query\tAA_pos_exon_query\tChr_query\tExon_coords_query\tStrand\tExon_hits_subject\tProt_subject\tAln_AA_pos_subject\t%Sim_aln_qry_sbj\t\%Id_aln_qry_sbj\tGap_number\t\%Gaps\tProt_subject\tExon_number_subject\tAA_pos_exon_subject\tChr_subject\tExon_coords_subject\tStrand\tSp_query\tSp_subject\n";
my @ks=sort(keys(%rex));
my $eid;
foreach $eid(@ks){
    print OUT "$rex{$eid}\n";
}
close (OUT);

###SUB-ROUTINES###
##I) SCORING EXON REGIONS
sub score_proteins {
    my ($n1,$n2,$seq1,$seq2)=@_; #as input we have the proteinIDs of the two sequences and the sequences of the two exons to be aligned.
    my (@s1,@s2);
    @s1=split //,$seq1; #get each character of exon1 sequence
    @s2=split //,$seq2; #get each character of exon2 sequence
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
    for ($i=0; $i<scalar(@s1); $i++){ #cycle on each position
	if ($sc{$s1[$i].$s2[$i]}){ #if there is a BLOSUM score for the pair
	    $e_score+=$sc{$s1[$i].$s2[$i]}; #add that to the e_score
	}
	if ($s1[$i] eq $s2[$i]){ #if the two positions are identical
	    $id_score++; #add 1 to the identity score
	    $sim_score++; #add 1 to the similarity score
	    $ng=0; #set number of gaps=0
	}
	elsif ($sim{$s1[$i].$s2[$i]}){ #if the two positions are not identical, but they are among the ones considered similar
	    $sim_score++; #add one to the similarity score
	    $ng=0; #set the number of gaps to zero
	}
	elsif ($s2[$i] eq "-" || $s1[$i] eq "-") { if ($ng==0) { $g++; $ng=1; } else { $e++; $ng=1;} } 
	if ($sc{$s1[$i].$s1[$i]}) { $m1+=$sc{$s1[$i].$s1[$i]}; } #Add the blosum score to m1 
	if ($sc{$s2[$i].$s2[$i]}) { $m2+=$sc{$s2[$i].$s2[$i]}; } #Add the blosum score to m2
    }
    if ($l1>0 && $l2>0){
   	$sim1=sprintf("%.2f",(($sim_score/$l1)*100));
    	$sim2=sprintf("%.2f",(($sim_score/$l2)*100));

    	$id1=sprintf("%.2f",(($id_score/$l1)*100));
        $id2=sprintf("%.2f",(($id_score/$l2)*100));
	#correct for weird cases when computing the global score
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
    }
    ##It will report the minimum score
    my ($simt, $idt, $glt);
    if ($sim1<$sim2) { $simt=$sim1; } else { $simt=$sim2; }
    if ($id1<$id2) { $idt=$id1; } else { $idt=$id2; }
    if ($global_score1<$global_score2) { $glt=$global_score1; } else { $glt=$global_score2; }
    #RETURNING RESULTS OF SIMILARITY
    return ($simt,$idt,$g,$glt);
}
##END OF SUBROUTINE I)
