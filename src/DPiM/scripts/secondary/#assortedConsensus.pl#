#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Path qw(make_path);
use File::Basename;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineAPMS);

{
    my $baseDir = '/home/glocke/DPiM/dpim4/witInstr';
    my $in = '/home/glocke/DPiM/dpim4/witInstr/cons1_01-30-2016.tsv';
    my $baseOut = basename($in);
    $baseOut =~ s/\.tsv//;
    
    my $cutoff = "$ENV{DPSCR}/supportCutoff.pl";
    
    my @perc = qw( Min 30 34 40 50 60 );
    my @networks;
    for my $p (@perc) {
	my $outDir = "$baseDir/percentile$p";
	make_path($outDir) unless -d $outDir;

	my $outFile = "$outDir/$baseOut"."_$p"."Percent.tsv";
	if ( -e $outFile) {
	    warn "found $outFile. skipping";
	}
	my $cmd = "$cutoff -in $in -out $outFile -percentile $p";
	if ($p eq 'Min') {
	    $cmd ="$cutoff -in $in -out $outFile -percentile 1";
	} elsif ($p == 34) {
	    $cmd ="$cutoff -in $in -out $outFile -minsupport 950";
	}
	say $cmd;
	system($cmd);
    }

    # map corum
}

exit;
