#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use File::Basename;

if (@ARGV != 3) {
    die "usage: $0 <real network out> <simulation network out> <output>";
}

sub compHGScoreCutoff {
    my ($realFile, $simFile, $outFile, $fdr) = @_;

    $fdr //= 0.05; # by default look for cutoff for 0.05 FDR

    $fdr = 100*$fdr - 1; #conver the FDR to column index from the file read in 


    my @real_edges = ();
    open(IN, "<$realFile");
    while (my $line = <IN>) {
	if ($line =~ "^FBgn"){
	    chomp($line);
	    my @field = split(/\t/, $line);
	    push @real_edges, $field[2];
	}
    }
    close(IN);

    my @simulated_edges = ();
    open(IN, "<$simFile");
    while (my $line = <IN>) {
	if ($line =~ "^FBgn"){
	    chomp($line);
	    my @field = split(/\t/, $line);
	    push @simulated_edges, $field[2];
	}
    }
    close(IN);

    my $FDR_cutoff = 0;
    my $FDR = 0;
    my $s = 0;
    my $r;
    while ($FDR < 0.05 && $s < $#simulated_edges) {
	$r = 0;
	while ($real_edges[$r] >= $simulated_edges[$s]) {
	    $r ++;
	}
	$FDR = ($s+1)/($r+1);
	$s ++;
    }
    $FDR_cutoff = $real_edges[$r+1];
    #$r = $r+1;
    my $cmd = "~/DPiM/scripts/scoreCutoff.pl -in $realFile -out $outFile -cutoff $FDR_cutoff -noheader ";
    system($cmd);
    #my $range = "8,".$r."p";
    #	`echo "qsub sed -n $range $realFile >> $outFile" | qsub -l h_rt=00:01:00`; # wrtie out the interactions in one file
    #`sed -n $range $realFile > $outFile`; # write out the interactions in one file
    return ($FDR_cutoff);
}

my $all_res = "All_Permutation_results.txt";

my $realFile = $ARGV[0];
my $simFile = $ARGV[1];
my $outFile = $ARGV[2];

my ($hg_cutoff, $num_sims_new) = compHGScoreCutoff($realFile, $simFile, $outFile);
#print "HGScore cutoff is $hg_cutoff\n";

my $sim = basename($realFile);
$sim =~ s/\.o\d+//;
print "$sim\t$hg_cutoff\n";

