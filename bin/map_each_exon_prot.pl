#!/usr/bin/perl
use warnings;
use strict;

##SCRIPT FOR MAP EACH EXON OF A PROTEIN IN SPECIES 1 TO ALL THE EXONS OF AN ORTHOLOGOUS PROTEINS IN SPECIES 2 AND VICEVERSA##

##
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
my $bin=$ARGV[12]; ##bin folder with alignintron pos script
my $bl=$ARGV[13]; ##blosum62  matrix
my $outf=$ARGV[14];
my ($dir,$project_dir)=$outf=~/(.+)\/(.+?)\/.+?\//;
my ($sp1,$sp2)=($s1,$s2);

my $exsc=$outf."/map_exons_".$s1."_".$s2."_part_".$part.".txt";##outputfile for exon scores
open (EXSC, ">$exsc");

#PRINTING HEADER
print EXSC "PID|GID_Sp1\tExon_#Sp1\tExon_coords\tExon_size\tExon_pos_aa\tExon_pos_aln\tPh_up_int\tPh_dw_int\t";
print EXSC "PID|GID_Sp2\tExon_#Sp2\tExon_coords\tExon_size\tExon_pos_aa\tExon_pos_aln\tPh_up_int\tPh_dw_int\t";
print EXSC "Dev_aln_pos_sp2-sp1\t%Sim\tScore_aln\tGap_num\tSp1\tSp2\n";
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
open (IN,"$i1") || die "Missing gene cluster file"; ##Gene cluster files
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

my (%gn, %cl1, %cl2);
open (IN,"$i1") || die "Missing gene cluster file"; ##Gene cluster files
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
### protein ids exons species 1
### ENSP00000000412|ENSG00000003056
open (IN,"$i2") || die "Missing protein ids exon file species 1"; 
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
### protein ids exons species 2
### ENSP00000000412|ENSG00000003056
open (IN,"$i3") || die "Missing protein ids exon file species 2";
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
my %exon;
my $rs1;
my %pos;
my @n;
open (IN,"$i4")|| die "Missing exon position files species 1";
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
open (IN,"$i5") || die "Missing exon position files species 2"; 
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
### Intron position files for each species
### Format: BL10724_cuf1|BL10724Sc0000009+Intron_13599249-3601287
my (%intron_aa,%intron_gg);
open (INFILE,"$i6"); ### Sp1
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
open (INFILE,"$i7"); ### Sp2
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

##### PROCESSING EXINT FILES
my (@t1, @t2);
my $f1=$i8;
my $f2=$i9;
my $m=0; 
my (%seqs, %tmp_int,%intron, %phin);
my ($ic, $x,$id,$z,$dev,$sid); 
my (@t3);
open (EXONE, "$f1") || die "Missing exint file species 1 ($f1)";
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
		$x=$z;
		@t3=split(/\./,$l[$z]);
		$id=$l[0]."\tintron_".$x;
		$intron_aa{$id}=$t3[0];
		$phin{$id}=$t3[1];
	    }
	}

    }
    elsif($m==1) {
	$seqs{$l[0]}=$sid."\n".$_;
    }
}
close (EXONE); 
########
open (EXTWO, "$f2")|| die "Missing exint file species 2 ($f2)";
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
		$x=$z;
		@t3=split(/\./,$l[$z]);
		$id=$l[0]."\tintron_".$x;
		$intron_aa{$id}=$t3[0];
		$phin{$id}=$t3[1];
	    }
	}
    }
    elsif($m==1) {
	$seqs{$l[0]}=$sid."\n".$_;
    }
}
close (EXTWO);

######EVALUATING THE PROTEIN ALIGNMENTS######
my @keys=keys(%pid1);
@keys=sort{$a<=>$b}(@keys);
my ($n1,$n2,$e,$i,$g, $zi, $zj);
my $el;
my $tels=scalar(@keys);
my @keys2=keys(%pid2);
@keys2=sort(@keys2);
my ($k, $name);
my (%seq);
my (@s1,@s2);
my (%score);
my $te=$outf."/tmp_part_".$part.".exint"; ##temporal exint file
my $tg=$outf."/tmp_part_".$part.".gde"; ##temporal gde file
my $tp=$outf."/pos_part_".$part.".txt"; ##temporal position
foreach $el (@keys){
    @t1=(); @t2=();    
    if ($pid1{$el}=~/\,/){ @t1=split(/\,/,$pid1{$el}); }
    else { push(@t1,$pid1{$el});  }    
    if ($pid2{$el}=~/\,/){ @t2=split(/\,/,$pid2{$el}); }
    else { push(@t2,$pid2{$el});  }    
    ##getting pairwise exint file and then make alignment
    print "$el\n"; ##printing in standard output the gene cluster ID that is being processed
    my $Gclid=$el;    
    for ($zj=0; $zj<scalar(@t1); $zj++){ ##opening FOR 1
	for ($zi=0; $zi<scalar(@t2); $zi++) {  ##opening FOR 2	        
	    ## 1) MAKING TEMPORAL EXINT FILE
	    open (TMPALN,">$te");
	    print TMPALN "$seqs{$t1[$zj]}\n$seqs{$t2[$zi]}\n"; 
	    close (TMPALN);
	    # RUNNING ALIGN INTRON POS
	    `perl $bin/AlignIntronPos.pl $te`; 	       
	    ## 2) OPENING OUTPUT GDE FILE
	    $n1=""; $n2=""; #print ">$t1[$zj]\t$t2[$zi]\n";
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
	    my ($sim1,$sim2,$id1,$glsc1)=(0,0,0,0);
	    ($sim1,$sim2,$id1,$glsc1)=score_proteins($n1,$n2,$seq1,$seq2);
	    if (($sim1>=20 && $sim2>=20) && ($sp1 ne $sp2) && !$score{$n1.",".$n2}){ 		
		##5) CALLING SUBROUTINE FOR SCORING INTRONS
		##  6) SCORING EACH EXON PAIR
		### 6.1) GETTING RESIDUES ALIGNMENT
		my $PP1=$pos{$n1}; my $PP2=$pos{$n2};
		if ($seq1 && $seq2 && $n1 && $n2 && $PP1 && $PP2){
		    my $tmp2=check_exons($seq1,$seq2,$n1,$n2,$sp1,$sp2,$glsc1,$sim1,$sim2,$id1,$Gclid,$PP1,$PP2);
		}
		else {} # print "3)$Gclid\t$n1\t$n2\t$seq1\t$seq2\t$PP1\t$PP2\n"; }
		$score{$n1.",".$n2}=1;
	    } 
	}
    }
}

########SUB-ROUTINES######

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
    my $tg=0;
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
	elsif ($s2[$i] eq "-" || $s1[$i] eq "-") { $tg++; if ($ng==0) { $g++; $ng=1; } else { $e++; $ng=1;} } 
	if ($sc{$s1[$i].$s1[$i]}) { $m1+=$sc{$s1[$i].$s1[$i]}; }
	if ($sc{$s2[$i].$s2[$i]}) { $m2+=$sc{$s2[$i].$s2[$i]}; }
    }    
    if ($l1>0 && $l2>0){
	$sim1=sprintf("%.2f",(($sim_score/$l1)*100));
	$sim2=sprintf("%.2f",(($sim_score/$l2)*100));    
	$id1=sprintf("%.2f",(($id_score/$l1)*100));
	$id2=sprintf("%.2f",(($id_score/$l2)*100));   
	$global_score1=($e_score+($g*-4)+($e*-1))/$m1;
	$global_score2=($e_score+($g*-4)+($e*-1))/$m2;   
	$global_score1=sprintf("%.2f",$global_score1);
	$global_score2=sprintf("%.2f",$global_score2);  
	##It will report the minimum score
	my ($simt, $idt, $glt);
	if ($sim1<$sim2) { $simt=$sim1; } else { $simt=$sim2; }
	if ($id1<$id2) { $idt=$id1; } else { $idt=$id2; }
	if ($global_score1<$global_score2) { $glt=$global_score1; } else { $glt=$global_score2; }    
	#RETURNING RESULTS OF SIMILARITY
	return ($sim1,$sim2,$global_score1,$tg);
    }
}
##II) CHECK ALL PROTEIN EXONS SP1 VS ALL PROTEIN EXONS SP2
sub check_exons {
    my ($seq1,$seq2,$n1,$n2,$sp1,$sp2,$glsc1,$sim1,$sim2,$id1,$el,$PP1,$PP2)=@_;   
    @s1=split(//,$seq1);
    @s2=split(//,$seq2);    
    my ($res1,$res2,$res)=(0,0);
    my %s1_s2_res=();
    my %s2_s1_res=();
    my %s1_s2_aln=();
    my %s2_s1_aln=();
    my %coords=();
    my ($se1, $se2)="";
    my ($k2, $pa);
    for ($i=0; $i<scalar(@s1); $i++){
	if ($s2[$i] ne "-") { $res2++; $k2=$res2; $se2.=$s2[$i]; }
	if ($s1[$i] ne "-"){ $res++; $res1++; $se1.=$s1[$i]; }
	$s1_s2_res{$res1}=$res2; #aa info
	$s2_s1_res{$res2}=$res1;
	$pa=$i+1;
	$s1_s2_aln{$res1}=$pa; ##aln info
	$s2_s1_aln{$res2}=$pa;
    }
    $res1++; $res2++;
    $pa=$i+1;
    $s1_s2_aln{$res1}=$pa; ##aln info
    $s2_s1_aln{$res2}=$pa;
    #ENSP00000416492|ENSG00000096093exon_4193-279chr652317486-52317746+
    my @p1=split(/\n/,$PP1); ##getting exon positions Prot1
    my @p2=split(/\n/,$PP2); ##getting exon positions Prot2
    my ($k, $l, $m, $sz1, $sz2, $rs1, $rs2, $as1, $as2, $bs1, $bs2, $e1, $e2, $ui1, $ud1, $ui2, $ud2, $dev, $inid);
    my (@exs1, @exs2, @l, @t1, @t2);
    for ($k=0; $k<scalar(@p1); $k++){
	@l=split(/\t/,$p1[$k]);
	$exs1[$k]=$l[2];
	$coords{$l[0]."\t".$l[1]}=$l[3].":".$l[4].":".$l[5];
    }
    for ($k=0; $k<scalar(@p2); $k++){
	@l=split(/\t/,$p2[$k]);
	$exs2[$k]=$l[2];
	$coords{$l[0]."\t".$l[1]}=$l[3].":".$l[4].":".$l[5];
    }
    for ($l=0; $l<scalar(@exs1); $l++) { ## Exons prot sp1 vs exons prot sp2
	for ($m=0; $m<scalar(@exs2); $m++) {
	    @t1=split(/\-/,$exs1[$l]);
	    @t2=split(/\-/,$exs2[$m]);
	    $sz1=$t1[1]-$t1[0]+1;
	    $sz2=$t2[1]-$t2[0]+1;	    
	    my $te1=substr($se1,($t1[0]-1),$sz1);
	    my $te2=substr($se2,($t2[0]-1),$sz2);	    
	    my $efile=$outf."/tmp_exons_part_".$part.".fa";
	    my $ofile=$outf."/tmp_exons_part_".$part.".aln";
	    #making a global alignment of only the given exons, one of each species
	    open (EX, ">$efile");
	    print EX ">$n1\n$te1\n>$n2\n$te2\n";
	    `mafft --quiet $efile > $ofile`;
	    open (OF, "$ofile");
	    my $cs=0;
	    my ($es1,$es2);
	    while (<OF>){
		chomp($_);
		if($_=~/>/){
		    $cs++;		    
		}
		elsif ($cs==1){
		    $es1.=$_;
		}
		elsif ($cs==2){
		    $es2.=$_;
		}
	    }
	    #print "$n1\t$n2#$es1\n*$es2\n";
	    my ($sime1,$sime2,$gle1,$ge1)=score_proteins($n1,$n2,$es1,$es2);	    
	    $rs1=$s1_s2_res{$t1[0]};
	    $rs2=$s1_s2_res{$t1[1]};
	    $as1=$s1_s2_aln{$t1[0]};
	    $as2=$s1_s2_aln{$t1[1]};
	    $bs1=$s2_s1_aln{$t2[0]};
	    $bs2=$s2_s1_aln{$t2[1]};	    
	    $e1=$l+1; $e2=$m+1;
	    if ($e1==1){ $ui1="NA"; } else { $inid=$n1."\tintron_".($e1-1); $ui1=$phin{$inid}; }
	    if ($e2==1){ $ui2="NA"; } else { $inid=$n2."\tintron_".($e2-1); $ui2=$phin{$inid};  }
	    if ($e1==(scalar(@exs1))){ $ud1="NA"; }  else { $inid=$n1."\tintron_".($e1); $ud1=$phin{$inid};  }
	    if ($e2==(scalar(@exs2))){ $ud2="NA"; } else { $inid=$n2."\tintron_".($e2); $ud2=$phin{$inid};  }
	    $dev=$bs1-$as1; ##deviation of aln position of exon sp1 vs exon sp2 
	    my $ie1=$n1."\texon_".$e1;
	    my $ie2=$n2."\texon_".$e2;
	    if ($dev>0){ $dev="+".$dev; }	    
	    print EXSC "$n1\texon_$e1\t$coords{$ie1}\t$sz1\t$exs1[$l]\t$as1-$as2\t$ui1\t$ud1\t";
	    print EXSC "$n2\texon_$e2\t$coords{$ie2}\t$sz2\t$exs2[$m]\t$bs1-$bs2\t$ui2\t$ud2\t";
	    print EXSC "$dev\t$sime1\t$gle1\t$ge1\t$sp1\t$sp2\n";
	}
    }
}

