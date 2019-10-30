#!/usr/bin/perl
#use warnings;
#use strict;


die "\nUsage: perl get_gcl_sp_pair.pl Gene_cluster_file\n\n" if !$ARGV[0];

$cl=$ARGV[0]; ##Gene cluster files

#GF000001	Mmu	ENSMUSG00000090353	Gm17555	OG0000001	O=0,R=0
open (CL, "$cl");
while (<CL>){
	chomp($_);
    	@l=split(/\t/,$_);
	$sp{$l[1]}=1;
	$csp{$l[0]}{$l[1]}=1;	
}


@k=sort(keys(%sp));
for($m=0; $m<scalar(@k); $m++){
	for($n=0; $n<scalar(@k); $n++){
		if ($m!=$n){
			if(!$pair{$k[$m]}{$k[$n]} && !$pair{$k[$n]}{$k[$m]}){
				$s1=$k[$m]; 
				$s2=$k[$n];
				if ($s1 && $s2){
					open (SCL, ">GCL_$s1-$s2.tab");
					open (CL, "$cl");
					while (<CL>){
						chomp($_);
    						@l=split(/\t/,$_);
						if ($csp{$l[0]}{$s1} && $csp{$l[0]}{$s2}){
							if ($l[1] eq $s1 || $l[1] eq $s2){
								print SCL "$_\n";
							}
						}
					}
					close (SCL);
					$pair{$k[$m]}{$k[$n]}=1;
				}
			}
		}
	}	
}





