#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readHeader readColsRef);
use DpimLib qw(networkHashFromEdgeList readHS);

# annotate a network so that it gives the recall of the discoverable edges*
#   in all droid sources
# * a discoverable edge is an edge in the droid network that connects two nodes
#   present in our network
#
# algorithm:
#   ingest our network
#   ingest all droid edges
#   take the subset of the droid edges connecting nodes in our network
#   find the column sums for each droid source
#   finally, go through our edges in order of high to low score and annotate
#     the recall of each source at each step of the way
#
# if -raw flag is set, merely report whether each edge is seen in the given
#   network

my %opts = getCommandLineOptions();

{
    my $netFile = $opts{in};
    my $droidFile = $opts{droid};
    my $out = $opts{out};
    

    my %droidNet;
    my %droidSum;
    readDroid(\%droidNet, $droidFile);
    pruneDroid(\%droidNet, \%droidSum, $netFile);

    my %recallSum;
    
    open my $IN, "<", $netFile or die "can't read $netFile. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 summing the recall of $droidFile among edges in $netFile";
    my @sources = qw(multidb flybase curagen_yth     finley_yth      hybrigenics_yth perrimon_apms   sc_interologs   hs_interologs   ce_interologs   all);
    say $OUT join "\t", qw(protein1 protein2 score), @sources;
    while (my $line = readHS($IN)) {
	chomp $line;
	my ($p1, $p2, $score) = (split /\s+/, $line);
	($p1, $p2) = sort ($p1, $p2);
	my @row;
	for my $src (@sources) {
	    $recallSum{$src}++ if exists $droidNet{$p1}{$p2}{$src};
	    push @row, $recallSum{$src} // 0;
	    $row[-1] /= $droidSum{$src} unless exists $opts{raw};
	}
	say $OUT join "\t", $p1, $p2, $score, @row;
    }
    close $OUT; 
    close $IN; 
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	droid => '/home/glocke/DPiM/droid/comprehensive.net',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString -raw >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "droid=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{droid});

    return %opts;
}

sub readDroid {
    my ($ret, $inFile) = @_;

    say "get header '$inFile'";
    my @cols = readHeader($inFile);
    my @read;
    say "readCols";
    readColsRef(\@read, $inFile, \@cols);
    shift @cols;
    shift @cols;
    for my $row (@read) {
	my $p1 = $row->{protein1};
	my $p2 = $row->{protein2};
	for my $c (@cols) {
	    $ret->{$p1}{$p2}{$c} = 1 if $row->{$c} > 0;
	}
    }

    return;
}

sub pruneDroid {
    my ($droid, $sum, $netFile) = @_;

    my %net;
    networkHashFromEdgeList(\%net, $netFile, undef, 'symmetric');

    for my $p1 (keys %$droid) {
	if (! exists $net{$p1}) {
	    delete $droid->{$p1};
	    next;
	}
	for my $p2 (keys %{ $droid->{$p1} }) {
	    if (! exists $net{$p2}) {
		delete $droid->{$p1}{$p2};
		next;
	    }
	    for my $src (keys %{ $droid->{$p1}{$p2} }) {
		$sum->{$src}++;
	    }
	}
    }
}
