#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);

# compare putative clusterings to each other
# version 04-08-2016 comparing different support cutoffs and DPIM1

{
    my $baseDir = '/home/glocke/DPiM/dpim4/withInstr/consensus4_nets';
    my $out = "$baseDir/clustStats.tab";
    my @supp = qw( support0001
support0500
support1000
support1500
support1800
support1900
support1980
support2000 
DPIM1);
    my @i = qw(1.2 1.4 1.6 1.8 2 3 4 6 10 15 20 30);
    
    #@supp = qw(support2000);
    #@i = qw(1.2 1.4 1.6);

    
    my %refs = map { $_ => "/home/glocke/DPiM/dpim4/withInstr/consensus4_nets/$_/$_.net.corum.allMin3.clust" } @supp;
    $refs{DPIM1} = "/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.corum.allMin3.clust";

    checkExist('f', $_) for values %refs;
    my @cols = ('stat');
    for my $supp (@supp) {
	for my $i (@i) {
	    my $putative = clustFile($supp, $i);
	    checkExist('f', $putative);
	    push @cols, "$supp-$i";
	}
    }

    my %tests = (
	pair => "/home/glocke/DPiM/scripts/secondary/pairwiseClusterStats.pl",
	jacc => "/home/glocke/DPiM/scripts/secondary/jaccardClusterStats.pl",
	onmi => "/home/glocke/DPiM/cpp/oNMI/wrapOnmi.pl",
	basic => "/home/glocke/DPiM/scripts/secondary/jaccardClusterStats.pl -basic",
	);
    
    #shift(@tests);
    #shift(@tests);
    
    my %results;
    my @stats; # that which was to be determined
    # one row for each cluster
    for my $supp (@supp) {
	warn "$supp\n";
	my $ref = $refs{$supp};
	for my $i (@i) {
	    warn "\t$i\n";
	    my $putative = clustFile($supp, $i);
	    checkExist('f', $putative);
	    for my $t (keys %tests) {
		warn "\t\t$t\n";
		my $test = $tests{$t};
		my $cmd = "$test -in $putative -mode mcl -ref $ref";
		#say $cmd;
		my @read = `$cmd`;
		chomp @read;
		$results{$t}{$supp}{$i} = \@read;
	    }
	}
    }

    
    
    
    #open my $OUT, ">", $out or die "can't write to $out. $!";
    say join "\t", @cols;
    writeBasic($results{basic}, \@supp, \@i);
    writePair($results{pair}, \@supp, \@i);
    writeJacc($results{jacc}, \@supp, \@i);
    writeOnmi($results{onmi}, \@supp, \@i);
    
    
    # print the transpose of the matrix
    my $maxRow = $#{ $stats[0] };
    for my $i (0..$maxRow) {
	say join " ", map { $_->[$i] } @stats;
    }
}

sub clustFile {
    my ($supp, $i) = @_;
    if ($supp eq 'DPIM1') {
	return "/home/glocke/DPiM/prevDPIM/DPIM1_MCL/mcl.clusters/mcl.clusters.DPIM1_scores.updateFBgn.ppi.abc-format.i$i.txt";	    
    }
    return "/home/glocke/DPiM/dpim4/withInstr/consensus4_nets/$supp/mcl2/mcl.clusters/mcl.clusters.$supp.ppi.abc-format.i$i.txt";
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

sub writeBasic {
    my ($stats, $supp, $ii) = @_;

    my @rows = qw(nCluster nProtein min max median);
    writer($stats, $supp, $ii, \@rows)

}

sub writePair {
    my ($stats, $supp, $ii) = @_;

    my @rows = qw(precision recall F MCC);
    writer($stats, $supp, $ii, \@rows)

}

sub writeJacc {
    my ($stats, $supp, $ii) = @_;

    my @rows = qw(precision recall F matches);
    writer($stats, $supp, $ii, \@rows)

}

sub writeOnmi {
    my ($stats, $supp, $ii) = @_;

    my @rows = qw(NMI_max NMI_lfk NMI_avg);
    writer($stats, $supp, $ii, \@rows)

}
sub writer {
    my ($stats, $supp, $ii, $rows) = @_;
    for my $r (0..$#$rows) {
	my @row = ($rows->[$r]);
	for my $s (@$supp) {
	    for my $i (@$ii) {
		push @row, sprintf("%.5f", $stats->{$s}{$i}[$r]);
	    }
	}
	say join "\t", @row;
    }
}
