#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);

# compare putative clusterings to each other
# version 02-01-2016 comparing different percentage support

{
    my $baseDir = '/home/glocke/DPiM/dpim4/witInstr';
    my @perc = qw( Min 30 34 40 50 60 );
    my %clusterLists = qw(
Min cluster_pMin.list
30 cluster_p30.list 
34 cluster_p34.list   
40 cluster_p40.list  
50 cluster_p50.list
60 cluster_p60.list 
);
    my %refs;
    {
	my @refs = readList('/home/glocke/DPiM/dpim4/witInstr/corum.list');
	%refs = map { $perc[$_] => $refs[$_] } 0..$#refs;
    }
    
    my @tests = (
	"/home/glocke/DPiM/scripts/secondary/pairwiseClusterStats.pl",
	"/home/glocke/DPiM/scripts/secondary/jaccardClusterStats.pl -verbose",
	"/home/glocke/DPiM/cpp/oNMI/wrapOnmi.pl",
	);
    
    #shift(@tests);
    #shift(@tests);
    
    my @stats; # that which was to be determined
    # one row for each cluster
    for my $perc (@perc) {
	my @covers = readList($baseDir."/".$clusterLists{$perc});
	@covers = sortByMCLStat(@covers);
	my $ref = $refs{$perc};
	say $perc;
	for my $putative (@covers) {
	    say "\t$putative";
	    my @results = ($perc, '');
	    for my $test (@tests) {
		my $cmd = "$test -in $putative -mode mcl -ref $ref";
		say $cmd;
		my @read = `$cmd`;
		chomp @read;
		push @results, @read, '';
	    }
	    push @stats, \@results;
	}
    }

    # print the transpose of the matrix
    my $maxRow = $#{ $stats[0] };
    for my $i (0..$maxRow) {
	say join " ", map { $_->[$i] } @stats;
    }
}

sub sortByMCLStat {
    my @covers = @_;

    my @sortable;
    for my $c (@covers) {
	$c =~ m/i([\d\.]+)\.txt/ or die "Can't parse $c";
	push @sortable, [$1, $c];
    }
    return map {$_->[1]} sort {$a->[0] <=> $b->[0]} @sortable
}

