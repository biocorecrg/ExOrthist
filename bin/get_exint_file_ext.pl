#!/usr/bin/env perl -w
use strict;
use Getopt::Long;

#Declaration of variables
my $gtf_folder;
my $genome_folder;
my $exons_db_folder="./";
my $species_string;
my $verboseFlag=1;
my $help;

#Arguments

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(  "GTF=s" => \$gtf_folder,
	     "G=s" => \$genome_folder,
	     "EX_DB=s" => \$exons_db_folder,
	     "sp=s" => \$species_string,
	     "verbose=s" => \$verboseFlag,
	     "h" => \$help,
	     "help" => \$help
    );


sub verbPrint {
    my $verbMsg = shift;
    unless ($verboseFlag == 0 || $verboseFlag eq "F" || $verboseFlag eq "FALSE") {
	chomp($verbMsg);
	print STDERR "[exorter mod I]: $verbMsg\n";
    }
}

if (!defined $genome_folder || !defined $gtf_folder || !defined $species_string || defined $help){
    die "\nUsage: Run_Module_I.pl -GTF path_to_gtfs/ -G path_to_genomes/ -sp Sp1,Sp2 [-EX_DB path_to_EXONS_DB/]

COMPULSORY
     -GTF            Path where GTFs are stored (they should be named Sp1_annot.gtf)
     -G              Path where gDNAs are stores (they should be named Sp1_gDNA.fasta)
     -sp Sp1,Sp2     String of species.

OPTIONAL
     -EX_DB          Path to EXONS_DB/ folder (default ./)
     -verbose T/F    Verbose (default TRUE) 
     -h/--help       This help message.

";
    
}

### just giving the info
verbPrint("EXONS_DB path set to $exons_db_folder\n");

### loops for each species
my @SPECIES = split(/\,/, $species_string);
foreach my $species (@SPECIES){
##Opening genome file
    my (%seq);
    my $chr;
    my @g;
    my $genome_file = "$genome_folder/$species"."_gDNA.fasta";
    open (GENOME, $genome_file) || die "Cannot open $genome_file for $species\n";
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
#Scaffold1348	protein_coding	exon	162899	162997	.	+	.	gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "5";
#Scaffold1348	protein_coding	CDS	162899	162997	.	+	.	gene_id "WHL22.100348"; transcript_id "WHL22.100348.0"; exon_number "5"; protein_id "WHL22.100348.0";

    my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
    my (%pid, %header, %ntseq);
    my $b=0;
    my (@l, @l1, @line, @l2, @l3, @s);
    my %rseq;
    my ($res, $size, $tmpseq);
    my %fex; ##saving the phase of the first exon
    my $gtf_file = "$gtf_folder/$species"."_annot.gtf";

    open (GTF,"$gtf_file") || die "Cannot open $gtf_file for $species\n";
    while (<GTF>){
	chomp($_); 
	@line=split(/\t/,$_);
	if ($seq{$line[0]}){
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
    open (EXINT_OUT, ">$exint_output") || die "Cannot open $exint_output\n";

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
    open (SIZE_OUT, ">$size_output") || die "Cannot open $size_output\n";
    print SIZE_OUT "ProteinID\tGeneID\tSize\n";
    
    open (EXINT_IN, $exint_output) || die "Cannot open $exint_output\n";
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
    open (L_P_OUTPUT, ">$longest_prot_output") || die "Cannot open $longest_prot_output\n";
    open (L_E_OUTPUT, ">$longest_exint_output") || die "Cannot open $longest_exint_output\n";
    
    print L_P_OUTPUT "GeneID\tProteinID\n";
    foreach my $temp_exint (sort keys %longest_prot){
	print L_P_OUTPUT "$temp_exint\t$longest_prot{$temp_exint}\n";
	print L_E_OUTPUT "$exint_data{$longest_prot{$temp_exint}}\n";
    }
    close L_P_OUTPUT;
    close L_E_OUTPUT;
}


### get_trs_gtf.pl could go straight here, although there could be some clashes with variables => OK if redo species loop
