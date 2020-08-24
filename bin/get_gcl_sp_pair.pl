#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Cwd qw(abs_path cwd);
use Data::Dumper;

my $cluster_file;

die "\nUsage: perl get_gcl_sp_pair.pl -f Gene_cluster_file\n\n" if !$ARGV[0];

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(  "f=s" => \$cluster_file);

# Format of gene cluster file:
#GF000001	Mmu	ENSMUSG00000090353	Gm17555	OG0000001	O=0,R=0
my %sp;
my %csp;

open (CL, $cluster_file) || die "Cannot open gene cluster file ($cluster_file)\n";
while (<CL>){
    chomp($_);
    my @l=split(/\t/,$_);
    $sp{$l[1]}=1;
    $csp{$l[0]}{$l[1]}=1;	
}
close (CL);

my %pair;
my @k=sort(keys(%sp));
for(my $m=0; $m<scalar(@k); $m++){
    for(my $n=0; $n<scalar(@k); $n++){
	if ($m!=$n){
	    if(!$pair{$k[$m]}{$k[$n]} && !$pair{$k[$n]}{$k[$m]}){
		my $s1=$k[$m]; 
		my $s2=$k[$n];
		if ($s1 && $s2){
		    open (SCL, ">$s1-$s2.cls.tab");
		    open (CL, $cluster_file) || die "Cannot open gene cluster file ($cluster_file)\n";
		    while (<CL>){
			chomp($_);
			my @l=split(/\t/,$_);
			if ($csp{$l[0]}{$s1} && $csp{$l[0]}{$s2}){
			    if ($l[1] eq $s1 || $l[1] eq $s2){
				print SCL "$_\n";
			    }
			}
		    }
		    #print Dumper \%csp;
		    close (SCL);
		    close (CL);
		    $pair{$k[$m]}{$k[$n]}=1;
		}
	    }
	}
    }	
}





