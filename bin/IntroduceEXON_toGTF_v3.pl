#!/usr/bin/perl

# started on 19/11/16
# unlike IntroduceEXON.pl, this scripts introduces non-annotated exons into the GTFs
# then these can be used to get the non-annotated proteins (fB).

die "\nUsage: IntroduceEXON_toGTF_v3.pl REFERENCE GTF Sp Dir_vastb_files genome_file\n\n" if !$ARGV[2];

$file_exons=$ARGV[0]; # REFERENCE table
$file_GTF=$ARGV[1];
$sp=$ARGV[2];
$dir=$ARGV[3]; ##path to vastdb files
$file_gDNA=$ARGV[4];


@f=glob("$dir/REFERENCE_SUPER_C1C2-$sp*.tab");
$file_SUPER = $f[0];
$file_NewIDs = "$dir/New_ID-$sp.txt";
$file_Genes = "$dir/$sp.Event-Gene.IDs.txt";
$file_ref_tr = "$dir/$sp"."_ref.txt";

### parse keys
open (IDs, $file_NewIDs) || die "*** DIE: Cannot find NewIDs\n\n";
while (<IDs>){
    chomp;
    @t=split(/\t/);
    next if $t[0]!~/EX/;
    $old_new{$t[1]}=$t[0];
}
close IDs;

open (GENES, $file_Genes) || die "*** DIE: Cannot find Events to Genes\n\n";
while (<GENES>){
    chomp;
    @t=split(/\t/);
    next if $t[0]!~/EX/;
    $Ev_Gene{$t[0]}=$t[1];
}
close GENES;

open (RefTR, $file_ref_tr) || die "*** DIE: Cannot find Ref tr\n";
while (<RefTR>){
    chomp;
    @t=split(/\t/);
    $ref_tr{$t[0]}=$t[1];
}
close RefTR;


### Parse gDNA
open (gDNA, $file_gDNA) || die "Can't get the genomic sequence: $file_gDNA";
while (<gDNA>){
    chomp;
    if (/\>/){
	($chr)=/\>(.+)/;
    }
    else {
	$seq{$chr}.=$_;
    }
}
close gDNA;

### parse superREF
open (SUPER, $file_SUPER) || die "*** DIE: Cannot find SuperREF file\n\n";
<SUPER>;
while (<SUPER>){
    chomp;
    @t=split(/\t/);
    $newID=$old_new{$t[0]};
    $newID=$t[0] if !$newID && $t[0]=~/EX/;
    next if !$newID;
    # head: EVENT C1_far C2_far C1_inc C2_inc C1_exc C2_exc
    # only INC are taken from here (otherwise, I'll take the REF one)
    $C1_inc{$newID}=$t[3];
    $C2_inc{$newID}=$t[4];
}
close SUPER;

#### Get the exons and their info
open (EXONS, $file_exons) || die "*** DIE: Cannot open $file_exons\n";
<EXONS>;
while (<EXONS>){
    chomp;
    @t=split(/\t/);
    $g=$Ev_Gene{$t[1]};
    $ev=$t[1];
    next if !$g;

    $C1_ref{$ev}=$t[8];
    $A_ref{$ev}=$t[9];
    $C2_ref{$ev}=$t[10];

    $co_ev{$t[9]}=$ev;
    ($chr,$i,$f)=$t[9]=~/(.+?)\:(.+?)\-(.+)/;
    $coA="$chr:$i";
    $coB="$chr:$f";
    $coA_ev{$coA}=$ev;
    $coB_ev{$coB}=$ev;
}
close EXONS;

#### Starts checking for annotation
open (GTF, $file_GTF) || die "*** DIE: Cannot open GTF\n";
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
	
#	push(@{$TR{$g}},$tr); # for each g, gets all tr
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

### If more GTFs are provided:
open (GTF2, $ARGV[3]); ######## GTF2
while (<GTF2>){
    chomp;
    @t=split(/\t/);

    ($trans_check)=$t[8]=~/transcript_id \"(.+?)\"/;
    next if $done{$trans_check} || !$trans_check;
    
    if ($t[2] eq "exon"){
	($g)=$t[8]=~/gene_id \"(.+?)\"/;
	
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
	
#	push(@{$TR{$g}},$tr); # for each g, gets all tr
	$TR{$g}{$tr}=1;
	$tr_co{$tr}{$co}=1; # co exists in that tr
	$tr_coA{$tr}{$coA}=1; # coA exists in that tr
	$tr_coB{$tr}{$coB}=1; # coB exists in that tr
	
	push(@{$tr_lines{$tr}},$_); # all lines of a transcript. Index 0 should be exon 1
	
	$str{$g}=$t[6];
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
close GTF2;
#################### GTF2

### v3 => vB to avoid confusions with Ensembl versions (10/02/17)
open (O, ">FakeTranscripts-$sp-vB.gtf");
open (LOG, ">LOG_FakeTranscripts-$sp-vB.tab");
print LOG "EVENT\tGENE_ID\tEXON\tFAKE_TR\tTYPE\tSOLUTION\n";
open (LOG2, ">STATS_FakeTranscripts-$sp-vB.tab");

#### Now loops through every exon, and focuses on the non-annotated
# Rules:
# 1) if C1inc and C2inc are present OK.
# ...
# In some cases, the exact exon is wrongly annotated, but there is no skipping form

foreach $ev (sort keys %A_ref){
    $g=$Ev_Gene{$ev};
    next if !$g;

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
#		print "$ev\t$A_ref{$ev}\t$t_C2\tOK\n";
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
			
#			if ($length_A%3==0){
#			    $comment="Alt_prot"; # all set
#			}
#			else {
#			    $out_frame=1; # and moves the problem to the "next" CDS
#			}
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
#	    print "$g\t$ev\t$A_ref{$ev}\n";
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
