#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsRef);
use DpimLib qw(readHS);

# given one network, find the scores given by another network for the same edge
# not found is defined as -1

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $cmp = $opts{cmp};
    my $out = $opts{out};

    my %cmpScore;
    my @cols = qw(node1 node2 score);
    my @input;
    readColsRef(\@input, $in, \@cols, 'line');
    for my $row (@input) {
	my ($n1, $n2) = sort ($row->{node1}, $row->{node2});
	$cmpScore{$n1}{$n2} = -1;
    }
    
    
    open my $IN, "<", $cmp or die "Can't read from $cmp. $!";
    my %missing;
    while (my $line = readHS($IN)) {
        $_ = $line;
        my @spl = split;
        my ($n1, $n2) = sort @spl[0..1];
        next unless exists $cmpScore{$n1}{$n2};
	$cmpScore{$n1}{$n2} = $spl[2];
    }

    my $header = join "\t", qw(node1 node2 score1 score2);
    my $preComments = "# started with $in\n# added scores from $cmp";

    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT $preComments;
    say $OUT $header;
    for my $row (@input) {
	chomp $row->{line};
	say $OUT $row->{line}, "\t", $cmpScore{$row->{node1}}{$row->{node2}};
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	cmp => '/home/glocke/DPiM/data/dpim3.12232014.nrBait.noCut.network',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in small.network -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "cmp=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{cmp});

    return %opts;
}
