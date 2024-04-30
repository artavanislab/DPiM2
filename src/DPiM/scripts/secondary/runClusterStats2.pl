#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);

# compare putative clusterings to each other

{
    my @mcode = qw'
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.05.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.1.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.25.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.5.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.75.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.9.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.95.list
';
    my @mcl = readList('/home/glocke/DPiM/MCL_network-enrichment-PPI/ppis-mcl-dpim3.09-25-2015.nrBait.77.44.network/mcl.clusters/short.list');
    checkExist('f', $_) for (@mcode, @mcl);

    my $ref = $mcode[0];
    #my $ref = '/home/glocke/DPiM/oldDpim/dpim3.1/corumClusters.maxJaccard0.9.tab';
    #my $ref = '/home/glocke/DPiM/oldDpim/dpim3.1/corum.plus.go.complexes.tab';

    #my $test = "/home/glocke/DPiM/scripts/secondary/pairwiseClusterStats.pl -ref $ref";    
    #my $test = "/home/glocke/DPiM/scripts/secondary/jaccardClusterStats.pl -ref $ref";
    #my $test = "/home/glocke/DPiM/scripts/secondary/overlappingNMIClusterStats.pl -ref $ref";
    my $test = "/home/glocke/DPiM/scripts/secondary/overlappingNMIClusterStats.pl -ref $ref -norm avg -clusterref";

    my @stats; # that which was to be determined
    
    for my $putative (@mcode) {
	my $cmd = "$test -in $putative";
	say $cmd;
	my @read = `$cmd`;
	chomp @read;
	push @stats, \@read;
    }

    for my $putative (@mcl) {
	my $cmd = "$test -in $putative -mode mcl";
	say $cmd;
	my @read = `$cmd`;
	chomp @read;
	push @stats, \@read;
    }

    my $maxRow = $#{ $stats[0] };
    for my $i (0..$maxRow) {
	say join " ", map { $_->[$i] } @stats;
    }
}
