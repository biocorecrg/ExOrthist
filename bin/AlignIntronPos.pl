#!/usr/bin/env perl

# From Scott's original script and adaptations in 2014
# Adapted and modified AGAIN on 08/03/16 --MI

use English;

### variables
$exint=$ARGV[0];
$aligner=$ARGV[1];
$cpus=$ARGV[2];
$aligner="MAFFT" if !$aligner;

#($root) = $exint =~ /([^\/]+)\./; # Exint file
($root) = $exint =~ /(.+)\./; # Exint file
print "\nRunning intron position alignment for $root...\n\n";

### with clustal
if ($aligner=~/clustal/i){
    system "clustalw2 $exint -output=gde";
}
elsif ($aligner=~/MAFFT/i){
    system "mafft --auto --thread $cpus --quiet $exint | sed 's/>/%/g' > $root.gde";  
}

open (OUT, ">$root.INT_ALIGN.aln") || die "Can't open output file\n";

$window = 10;
@code = ("u", "u", "u", "v", "w", "x", "y", "z", "z", "z", "z");
# u = 0-30%, v = 30-40%, w = 40-50%, x =50-60%, y = 60-70%, z = 70-100%.

$alength = 140;
%c = (0,0,1,1,2,2,"","*");

#  GET INTRONS FROM .EXINT FILE
$/ = ">";
open (EXINT, $exint);
while (<EXINT>) {
    ($gene, $introns) = /^\s*(\S+)\s+\[?([ \t\d\.]+)\]?/;
    $introns{$gene} = $introns;
    print OUT "$gene\n";
}
print OUT "\n";
close EXINT;

# OPEN THE ALN GDE
open (ALN, "$root.gde");
$_ = <ALN>;
%ind = ();
%seq = ();
%S = ();
%intpos = ();
%m = ();
%pairwise = ();

$ch = 0;
while (s/\%\s*(\S+).*\n([A-Za-z\-\n\s]+)(\n\%|\s*$)/$3/) {
    ($gene, $seq) = ($1,$2);
    $gene =~ /.{0,6}/;
    $list[$ch++] = $gene;
    $seq{$gene} = $seq;
    $seq{$gene} =~ s/\s//g;
    $seq{$gene}=~tr/a-z/A-Z/; # makes it all big caps
    @{$S{$gene}} = split (//, $seq{$gene});
    foreach $intron ($introns{$gene} =~ /\d+\.\d/g) {
	($pos, $phase) = $intron =~ /(\d+)\.(\d)/;
	$intpos{$gene}[$pos] = $phase;
    }
}

$length = length $seq{$gene};

#FILTERING HAPPENS HERE
for $i (0..$#list) {
    for $j ($i+1..$#list) {
	($gene, $gene2) = @list[$i,$j];
	foreach $pos (0..$length-1) {
	    if (($S{$gene}[$pos] =~ /[a-zA-Z]/) && ($S{$gene2}[$pos] =~ /[a-zA-Z]/)) { # only allowed small caps?? => it should be the opposite?
		$m{"$gene:$gene2"} .= int ($S{$gene}[$pos] eq $S{$gene2}[$pos]);
	    }
	    else {$m{"$gene:$gene2"} .= " ";}
	}
	$m = $m{"$gene2:$gene"} = $m{"$gene:$gene2"};
    }	      	       
    for $j (0..$#list) {
	($i == $j) && next;
	($gene, $gene2) = @list[$i,$j];
	for $pos ($window..$length-1-$window) {
	    $m{"$gene:$gene2"} =~ /.{$pos}/;
	    ($up, $down) = ($&, $POSTMATCH);
	    $up =~ s/\s+//g; $down =~ s/\s+//g;
	    ($up =~ /(.{$window})$/) || next;
	    $upregion = $1;
	    $up = $upregion =~ tr/1/1/;
	    ($down =~ /(.{$window})/) || next;
	    $downregion = $1;
	    $down = $downregion =~ tr/1/1/;
	    $pairwise{"$gene:$gene2"}[$pos] = 
		$pairwise{"$gene2:$gene"}[$pos] = min ($up, $down);
	}
    }		   
}

@cons = ();
for $i (0..$window, $length-$window..$length) {$cons[$i] = " ";} # it gives " " to the first and last 10 positions of the Aln
for $pos ($window..$length-1-$window) { # now it works on the remaining aln
    $min = 50;
    for $i (0..$#list) { # $#list only has two elements here (geneA and geneB)
	for $j ($i+1..$#list) {
	    $min = min ($min, $pairwise{"$list[$i]:$list[$j]"}[$pos]);
	}
    }
    $cons[$pos] = $code[$min];
    $tally[$min]++;
}

print OUT "INTRON_CODES: ";
for $i (0..$length-1) {
    foreach $gene (keys %seq) {
	$ind{$gene} += ($S{$gene}[$i] =~ /[a-zA-Z]/);
    }
    $is_there_an_int = 0;
    foreach $gene (keys %seq) {
	$is_there_an_int += ($intpos{$gene}[$ind{$gene}] ne "");
    }
    if ($is_there_an_int) {
	$shared = "";
	$intcode = 0;;
	foreach $g (0..$#list) {
	    $gene = $list[$g];
	    if ($intpos{$gene}[$ind{$gene}] =~ /\d/) {
		    $intcode += 2**$g;}
	    $S{$gene}[$i] .= " $c{$intpos{$gene}[$ind{$gene}]} ";
	    $shared .= $c{$intpos{$gene}[$ind{$gene}]};		       
	    $intpos{$gene}[$ind{$gene}] = "";
	}
	
	while ($shared =~ s/(.)\1/$1/) {}
	if ((length $shared) > 1) {$shared = "*";}
	if ($cons[$i] =~ /[wxyz]/) {
	    $cons[$i] .= "_" . "$shared ";
	    print OUT "$intcode ";
	}
	else {$cons[$i] .= " $shared ";}
    }
}
foreach $gene (keys %seq) {
    $seq{$gene} = join ("", @{$S{$gene}});
}

$cons = join ("", @cons);

print OUT "@tally\n";
print OUT "> $root\n\n";
while (length ($seq{$list[0]})) {
    foreach $gene (@list) {
	($spec) = $gene =~ /(.{0,6})/; # number of letters in name
	$seq{$gene} =~ s/.{1,$alength}//;
	print OUT "$spec\t$&\t$spec\n";
    }
    print OUT "\n";
    $cons =~ s/.{1,$alength}//;
    print OUT "\t$&\n\n";
}


#### Subrutines
sub max {
    if ($_[0] < $_[1]) {$_[1]}
    else {$_[0]}
}

sub min {
    if ($_[0] < $_[1]) {$_[0]}
    else {$_[1]}
}
