#!/usr/bin/perl

$i1=$ARGV[0];
$i2=$ARGV[1];
##change EXONS_DB accordingly!!!
$e1="/users/mirimia/ymarquez/ANALYSIS/ORTOLOGUES/EXORTH_PAPER/NOVA/EXONS_DB/Hs2/Hs2_prot_exons_overlap.txt";
$e2="/users/mirimia/ymarquez/ANALYSIS/ORTOLOGUES/EXORTH_PAPER/NOVA/EXONS_DB/Mm2/Mm2_prot_exons_overlap.txt";
$e3="/users/mirimia/ymarquez/ANALYSIS/ORTOLOGUES/EXORTH_PAPER/NOVA/EXONS_DB/Dme/Dme_prot_exons_overlap.txt";

#GF00001BtaENSBTAG00000045550TSPAN6
open (LIST, $i1);
while (<LIST>){ ##checking list of consistent events
    chomp($_);
    @l=split(/\t/,$_); 
    if ($l[4] ne "Intronless"){
	$cl{$l[2]}=$l[0];
	$sp{$l[2]}=$l[1];
    }
}

#OV_EX_Mm2_186054ENSMUSG00000067377133892701-133892766
open (LIST, $e1);
while (<LIST>){ ##checking list of consistent events
    chomp($_);
    @l=split(/\t/,$_);
    if ($cl{$l[1]}){
	$s=$sp{$l[1]};
	if (!$ex{$l[1]}{$l[0]}){
	    $ex{$l[1]}{$l[0]}=1;
	    $tex{$s}++;
	}
	$id=$l[1]."\t".$l[2];
	$ov{$id}=$l[0];
    } 
}
open (LIST, $e2);
while (<LIST>){ ##checking list of consistent events
    chomp($_);
    @l=split(/\t/,$_);
    if ($cl{$l[1]}){
	$s=$sp{$l[1]};
	if (!$ex{$l[1]}{$l[0]}){
	    $ex{$l[1]}{$l[0]}=1;
	    $tex{$s}++;
	}
	$id=$l[1]."\t".$l[2];
	$ov{$id}=$l[0];
    } 
}
open (LIST, $e3);
while (<LIST>){ ##checking list of consistent events
    chomp($_);
    @l=split(/\t/,$_);
    if ($cl{$l[1]}){
	$s=$sp{$l[1]};
	if (!$ex{$l[1]}{$l[0]}){
	    $ex{$l[1]}{$l[0]}=1;
	    $tex{$s}++;
	}
	$id=$l[1]."\t".$l[2];
	$ov{$id}=$l[0];
	#print "$id\t$l[0]\n";
    } 
}

open (OUT, ">Exon_cluster_NR.tab");
open (STATS, ">Stats_exon_cluster.txt");
#GF00002.003ENSMUSG00000031250chrX:133859723-133859863:+Mm2
open (LIST, $i2);
while (<LIST>){ ##checking list of consistent events
    chomp($_);
    @l=split(/\t/,$_);
    if ($cl{$l[1]}){ 
	@t1=split(/\:/,$l[2]);
	$id2=$l[1]."\t".$t1[1];
	$oid=$ov{$id2};
	#print "#$oid#\t$id2\n";
	if (!$nex{$oid}){
	    $clex{$l[3]}++;
	    $nex{$oid}=1;
	    print OUT "$_\n";
	}
	
    }
    
}


@keys=sort(keys(%tex));
foreach $el (@keys){
    $av=($clex{$el}/$tex{$el})*100;
    $av=sprintf("%.2f",$av);
    print STATS "Annotated exons $el\t$tex{$el}\tExons in clusters\t$clex{$el}\t$av%\n";
}




