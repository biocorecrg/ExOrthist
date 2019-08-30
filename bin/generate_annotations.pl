#!/usr/bin/env perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path cwd);

#Declaration of variables
my $gtf_folder;
my $genome_folder;
my $exons_db_folder="./";
my $species_string;
my $verboseFlag=1;
my $help;
my $vastdb_refs;

#Arguments

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(  "GTF=s" => \$gtf_folder,
	     "G=s" => \$genome_folder,
	     "EX_DB=s" => \$exons_db_folder,
	     "sp=s" => \$species_string,
	     "vastdb=s" => \$vastdb_refs,
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

if (!defined $genome_folder || !defined $gtf_folder || !defined $species_string || defined $help){
    die "\nUsage: generate_annotations.pl -GTF path_to_gtfs/ -G path_to_genomes/ -sp Sp1,Sp2 [-EX_DB path_to_EXONS_DB/ -vastdb REF1,REF2]

Script that creates all annotation files needed for the second module of the pipeline

COMPULSORY
     -GTF              Path where GTFs are stored (they should be named Sp1_annot.gtf)
     -G                Path where gDNAs are stores (they should be named Sp1_gDNA.fasta)
     -sp Sp1,Sp2       String of species.

OPTIONAL
     -vastdb r1,r2     Comma-separated list of VASTDB reference files (must match species list in -sp)
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

### Gets @SPECIES array
my @SPECIES = split(/\,/, $species_string);

### prepares VASTDB references files
my %VASTDB_files;
if ($vastdb_refs){
    my @temp_vts=split(/\,/,$vastdb_refs);
    die "You need to provide the same number of comma-separated elements for -vastdb and -sp\n" if $#SPECIES != $#temp_vts;
    foreach my $i (0..$#SPECIES){
	$VASTDB_files{$SPECIES[$i]}=$temp_vts[$i];
    }
    verbPrint("VASTDB Reference files provided\n");
}

### loops for each species (Loop 1: previous get_ref_proteins.pl)
foreach my $species (@SPECIES){
    # In this first loop, it checks if the folders for each species exist or create them
    system "mkdir $exons_db_folder/$species" unless (-e "$exons_db_folder/$species");

    # Opening genome file
    my (%seq);
    my $chr;
    my @g;
    my $genome_file = "$genome_folder/$species"."_gDNA.fasta";
    open (GENOME, $genome_file) || die "Cannot open $genome_file for $species (loop 1)\n";
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

    open (GTF,"$gtf_file") || die "Cannot open $gtf_file for $species (loop 1)\n";
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
#    print SIZE_OUT "ProteinID\tGeneID\tSize\n"; # maintain without header, as the original
    
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
}


### loops for each species (loop 2: get_trs_gtf.pl)
foreach my $species (@SPECIES){
    my (@line, @coords1, @coords2, @gene, @l1, @l2, @l3);
    my ($l, $grep, $gid, $int, $r, $prot, $tmp, $trid, $n);
    my (%tr1, %tr2, %tpid);
    my %prot;
    my (%coords1, %coords2);
    my (%strand1, %strand2);
    my %intron;
    my (%chr1, %chr2);

    my $infile1 = "$gtf_folder/$species"."_annot.gtf";
    open (INFILEONE, $infile1) || die "Cannot open $infile1 (loop 2)\n"; 
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


### Only if VASTDB files provided
if ($vastdb_refs){
    ### loops for each species (Loop 3: get_prot_isof_exon_2.pl)
    #   -i1 Hsa_prot_sizes.txt -i2 Hsa_tr_coords_CDS.txt -i3 REFERENCE-ALL_ANNOT-Hs2136.tab 
    #   -o1 Ref_protein_exons_Hsa_1.txt -o2 missing_exons.txt
    foreach my $species (@SPECIES){  
        # Declaration of variables
	my (@line);
	my %pid;
	my %lprot;
	my ($ex, $l, $intron, $ex2);
	my (@t1,@t2, @c);
	my %psize;
	my $pr;
	
	# Skips the step for the species if NA was provided
	next if $VASTDB_files{$species} eq "NA";
	
	# Declaration of input files
	my $i1 = "$exons_db_folder/$species/$species"."_prot_sizes.txt";
	my $i2 = "$exons_db_folder/$species/$species"."_tr_coords_CDS.txt";
	my $i3 = $VASTDB_files{$species};
	
	# Declaration of output files
	my $out = "$exons_db_folder/$species/Ref_protein_exons_".$species."_1.txt";
	my $out2 = "$exons_db_folder/$species/missing_exons.txt";
	open (OUT,">$out") || die "Cannot open $out (loop 3)\n";
	open (MISS,">$out2") || die "Cannot open $out2 (loop 3)\n";
	
	verbPrint("Generating Ref_protein_exons_".$species."_1.txt and missing_exons.txt from VASTDB\n");

	open (INONE, $i1) || die "Cannot open $i1 (loop 3)\n";
        # Format: BL21561_evm0	BL21561	476
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
	#Format: BL20899_cuf1|BL20899	Sc0000317	-	255783,255866-258240,258374-259202,259294
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
	
	my ($id,$id4, $b, $trcoords);
	my @ex;
	open (INTHREE, $i3) || die "Cannot open $i3 (loop 3)\n"; 
	# Format: FAM13A	BlaEX0015343	Sc0000095:756974-757072	99	Sc0000095:757428,756974-757072,756243	S	
	#         Sc0000095:757428,756974-757072,756243:-	99	
	#         Sc0000095:757428-757560	Sc0000095:756974-757072	Sc0000095:756125-756243
	while (<INTHREE>){ 
	    chomp($_);		
	    @line=split(/\t/,$_);
	    if ($line[1]=~/EX/){ ##processing only exons
		if ($pid{$line[6]}){
		    $trcoords=$coords{$pid{$line[6]}};
		    print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$pid{$line[6]}\tAnnotated_exon_C1,A,C2\t$trcoords\n";
		}
		elsif ($pid{$line[2]}) {
		    $trcoords=$coords{$pid{$line[2]}};
		    print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$pid{$line[2]}\tAnnotated_exon_A\t$trcoords\n";
		}
		else {
		    @t1=split(/\,/,$line[6]);	
		    $id=$t1[0].",".$t1[2]; ##searching then only the intron, the exon might be not annotated
		    if ($pid{$id}){
			$trcoords=$coords{$pid{$id}};
			print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$pid{$id}\tNot_annotated_exon|intron:$id\t$trcoords\n"; 
		    }	
		    else { print MISS "$_\n"; }
		}
	    }
	}
	close (INTHREE);
	close OUT;
	close MISS;
    }
    
    
    ### loops for each species (Loop 4: get_prot_isof_exon_3.pl) only if VASTDB refs
    #   -sp Hsa 
    #   -i1 Hsa_prot_sizes.txt -i2 Hsa_trid_protid.txt -i3 Hsa_tr_coords.txt -i4 Ref_protein_exons_Hsa_1.txt -i5 missing_exons.txt 
    #   -o1 Ref_protein_exons_Hsa_2.txt -o2 missing_exons_2.txt -o3 Final_exons_ref_proteins_Hsa.txt -o4 Hsa_vast_exons_prot_ids.txt
    foreach my $species (@SPECIES){
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
	next if $VASTDB_files{$species} eq "NA";

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
	my $out4 = "$exons_db_folder/$species/$species"."_vast_exons_prot_ids.txt";
	my $mex = "$exons_db_folder/$species/$species"."_not_annotated_exons.txt";
	
	open (OUT,">$out") || die "Cannot open $out (loop4)\n";
	open (MISS,">$out2") || die "Cannot open $out2 (loop4)\n";
	open (OUTF,">$out4") || die "Cannot open $out4 (loop 4)\n"; # out3 is done by a cat
	open (MEX, ">$mex") || die "Cannot open $mex (loop 4)\n";
	
	### Parsing first input
	open (INONE, $i1) || die "Cannot open $i1 (loop 4)\n";
	# Format: BL21561_evm0	BL21561	476
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
	# Format: BL20899_cuf1|BL20899	Sc0000317	-	255783,255866-258240,258374-259202,259294
	while (<INTHR>){ 
	    chomp($_);		
	    @line=split(/\t/,$_);
	    $s=$line[3]; $s=~s/\,/\|/g; $s=~s/\-/\,/g; $s=~s/\|/\-/g;
	    @n=split(/\|/,$line[0]);
	    ### This gave WARNINGS: basically, nc transcripts do not have prot?
#	    $trid{$n[0]}="NA" if (!$trid{$n[0]}); # this solution adds lines in the output for nc exons
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
			if ($psize{$trid{$n[0]}} > $lprot{$ex}){
			    $pid{$ex}=$trid{$n[0]};
			    $lprot{$ex}=$psize{$trid{$n[0]}};
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
			if ($psize{$trid{$n[0]}} > $lprot{$intron}){
			    $pid{$intron}=$trid{$n[0]};
			    # This gave warnings (likely for nc transcripts) 
			    $lprot{$intron}=$psize{$trid{$n[0]}} if (defined $trid{$n[0]});
			}
		    }
		}
	    }	
	}
	close (INTHR);
	
	### Parsing fifth input
	verbPrint("Generating Ref_protein_exons_$species"."_2.txt\n");
	verbPrint("Generating missing_exons_2.txt from VASTDB");
	
	my ($id,$g, $b, $trcoords);
	my @ex;
	open (INFIVE, $i5) || die "Cannot open $i5 (loop 4)\n";
	# Format: FAM13A	BlaEX0015343	Sc0000095:756974-757072	99	Sc0000095:757428,756974-757072,756243	S
	#         Sc0000095:757428,756974-757072,756243:-	99	
	#         Sc0000095:757428-757560	Sc0000095:756974-757072	Sc0000095:756125-756243
	while (<INFIVE>){ 
	    chomp($_);		
	    @line=split(/\t/,$_);
	    if ($pid{$line[6]}){
		$g=$gene{$pid{$line[6]}};
		$trcoords=$coords{$pid{$line[6]}};
		print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$pid{$line[6]}|$g\tAnnotated_exon_C1,A,C2*\t$trcoords\n";
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
	verbPrint("Generating $species"."_vast_exons_prot_ids.txt");
	verbPrint("Generating $species"."_not_annotated_exons.txt");
	
	my %print;
	open (IN, $out3) || die "Cannot open $out3 (loop 4)\n";
	while (<IN>){ 
	    chomp($_);		
	    @line=split(/\t/,$_);
	    if (!$print{$line[6]}){
		print OUTF "$line[6]\n";
		$print{$line[6]}=1;
	    }
	    if ($line[7]=~/Not_annotated_exon/){
		print MEX "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
	    }	
	}
	close (MEX);
	close (OUTF);
	close IN;
    }
}

### loops for each species (Loop 5: get_aa_pos_exon.pl)
#   -i Hsa_annot_exons_prot_ids.txt –v Hsa_vast_exons_prot_ids.txt -GTF Hsa_annot.gtf 
#   -out Hsa_protein_ids_exons_pos.txt
foreach my $species (@SPECIES){  
    #Declaration of variables
    my (%pid, %protein);
    my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m, $ex);
    my (%header, %ntseq);
    my $b=0;
    my (@l, @l1, @line, @l2, @l3, @s);
    
    # Declaration of input files
    my $p = "$exons_db_folder/$species/$species"."_annot_exons_prot_ids.txt";
    my $v = "$exons_db_folder/$species/$species"."_vast_exons_prot_ids.txt" if (-e "$exons_db_folder/$species/$species"."_vast_exons_prot_ids.txt");
    my $gtf_file = "$gtf_folder/$species"."_annot.gtf";
    
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
    open (GTF,"$gtf_file") || die "Cannot open $gtf_file (loop 5)\n";
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


### loops for each species (Loop 6: get_genomic_coords_intron.pl
#   -p Hsa_annot_exons_prot_ids.txt -v Hsa_vast_exons_prot_ids.txt -t Hsa_trid_protid.txt 
#   -c1 Hsa_tr_coords.txt -c2 Hsa_tr_coords_CDS.txt 
#   -o1 Hsa_protein_ids_intron_pos.txt -o2 Hsa_protein_ids_intron_pos_CDS.txt
foreach my $species (@SPECIES){  
    # Declaration of variables
    my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
    my (%pid, %tid, %ntseq);
    my $b=0;
    my (@l, @l1, @line, @l2, @l3, @coords);
    my ($l, $size, $tmpseq, $c);

    # Declaration of input files
    my $p = "$exons_db_folder/$species/$species"."_annot_exons_prot_ids.txt";
    my $v = "$exons_db_folder/$species/$species"."_vast_exons_prot_ids.txt" if (-e "$exons_db_folder/$species/$species"."_vast_exons_prot_ids.txt");
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
    # Format: ENSBTAT00000000005	ENSBTAP00000000005
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


### loops for each species (Loop 7: get_ex_size_phase.pl) => try to intergrate in exint?
# -GTF Hsa_annot.gtf -out Hsa_annot_exon_sizes.txt
foreach my $species (@SPECIES){  
    # Declaration of variables (same as loop 1)
    my ($j,$tmp, $gid, $prot, $phase, $nuc, $aa, $pid, $tfile, $m);
    my (%pid, %header, %ntseq);
    my $b=0;
    my (@l, @l1, @line, @l2, @l3, @s);
    my %rseq;
    my ($res, $size, $tmpseq);
    my %fex; ##saving the phase of the first exon

    my $outfile = "$exons_db_folder/$species/$species"."_annot_exon_sizes.txt";
    open (OUT, ">$outfile") || die "Cannot create $outfile (loop 7)\n";
    
    verbPrint("Generating $species"."_annot_exon_sizes.txt\n");

    ### Parsing GTF file
    my $gtf_file = "$gtf_folder/$species"."_annot.gtf";
    open (GTF,"$gtf_file") || die "Cannot open $gtf_file for $species (loop 1)\n";
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
}


