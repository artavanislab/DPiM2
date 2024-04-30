#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);

# simple script that produces cluster statistics such that they can be easily 
# copied into excel

my %opts = getCommandLineOptions();

{
    my @mcode = qw(
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.05.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.1.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.25.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.5.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.75.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.9.list
/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/d0.95.list
);
    my @mcl = readList('/home/glocke/DPiM/MCL_network-enrichment-PPI/ppis-mcl-dpim3.09-25-2015.nrBait.77.44.network/mcl.clusters/short.list');
    my @louvain = readList('/home/glocke/DPiM/oldDpim/dpim3.1/louvain/list.list');
    checkExist('f', $_) for (@mcode, @mcl, @louvain);

    my %refFiles = (
	corum => '/home/glocke/DPiM/oldDpim/dpim3.1/corumClusters.tab',
	corumGO => '/home/glocke/DPiM/oldDpim/dpim3.1/corum.plus.go.complexes.tab',
	'corum0.9' => '/home/glocke/DPiM/oldDpim/dpim3.1/corumClusters.maxJaccard0.9.tab',
	GO => '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOComplexes.min4',
    );
    my $ref = $refFiles{$opts{ref}};

    my %tests = (
	pair => "/home/glocke/DPiM/scripts/secondary/pairwiseClusterStats.pl -ref $ref",
	jaccard => "/home/glocke/DPiM/scripts/secondary/jaccardClusterStats.pl -ref $ref",
	onmi => "/home/glocke/DPiM/cpp/oNMI/wrapOnmi.pl -ref $ref",
	);
    my $test = $tests{$opts{test}};
    #my $test = "/home/glocke/DPiM/scripts/secondary/overlappingNMIClusterStats.pl -ref $ref -norm avg"; # must be run 3 times for different normalizations

    my @stats; # that which was to be determined
    if (0) {
    {
	# truth vs truth
	my $cmd = "$test -in $mcode[0] -mode truthonly";
	say $cmd;
	my @read = `$cmd`;
	chomp @read;
	push @stats, \@read;
    }
    
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
    }

    for my $putative (@louvain) {
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

sub getCommandLineOptions {
    # test:
    my @tests = qw(pair jaccard onmi);
    my %tests = map {$_ => 1} @tests;
    my @arr = @tests;
    $arr[0] = "*".$arr[0]."*";
    my $testString = "-test ".join("/", @arr);

    # reference data:
    my @refs = qw(corum corum0.9 corumGO GO);
    my %refs = map {$_ => 1} @refs;
    @arr = @refs;
    $arr[0] = "*".$arr[0]."*";
    my $refString = "-ref ".join("/", @arr);

    my %defaults = (
	);
    my $defaultString = join " ", map { "-$_ $defaults{$_}" } keys %defaults;
    my $usage = "usage: $0 $testString $refString \n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "test=s", "ref=s");
    die $usage unless exists $opts{ref} && exists $opts{test};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    die "you must select one of the following tests: ", 
        join(", ",keys(%tests)), "\n" 
	if (exists $opts{test} && ! exists $tests{$opts{test}});
    $opts{test} //= $tests[0];

    die "you must select one of the following reference datasets: ", 
        join(", ",keys(%refs)), "\n" 
	if (exists $opts{ref} && ! exists $refs{$opts{ref}});
    $opts{ref} //= $refs[0];

    return %opts;
}
