#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);

# given a network and a list of proteins, pull out any edge that contains
# only proteins in the list

my %opts = getCommandLineOptions();

{
    my $netFile = $opts{net};
    my $nodeString = $opts{nodes};
    my $out = $opts{out};
    my $testString = $opts{test}; # if defined, only report edges between
    my @grepCols = split ',', $opts{grepcols};
    # these nodes and $opts{nodes}, along with edges among these

    my %nodes;
    {
	my @nodes = split ",", $nodeString;
	if (@nodes == 1) {
	    $_ = $nodeString;
	    @nodes = split;
	}
	%nodes = map {$_ => 1} @nodes;
    }
    my %testNodes;
    if (defined $testString) {
	my @testNodes = split ",", $testString;
	if (@testNodes == 1) {
	    $_ = $testString;
	    @testNodes = split;
	}
	%testNodes = map {$_ => 1} @testNodes;
    }

    my $testRow;
    if (defined $testString) {
	$testRow = sub {
	    my ($n1, $n2) = @_;
	    return (exists $nodes{$n1} && exists $testNodes{$n2}) ||
		(exists $nodes{$n2} && exists $testNodes{$n1});
	}
    } elsif (exists $opts{inclusive}) {
	$testRow = sub {
	    my ($n1, $n2) = @_;
	    return exists $nodes{$n1} || exists $nodes{$n2};
	}
    } else {
	$testRow = sub {
	    my ($n1, $n2) = @_;
	    return exists $nodes{$n1} && exists $nodes{$n2};
	}
    }

    open my $IN, "<", $netFile or die "can't read from $netFile. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";

    #if (@testNodes > 0) 
    while (my $line = <$IN>) {
	$_ = $line;
	my ($n1, $n2) = (split)[@grepCols];
	if ($testRow->($n1, $n2) || $line =~ /^protein1/) {
	    print $OUT $line;
	}
    }
    close $IN;
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	grepcols => '0,1',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -net net.in -nodes node1,node2,... -out output < ".
	" $defaultString -test test1,test2 -inclusive >\n";

    my %opts = ();
    GetOptions(\%opts, "net=s", "nodes=s", "out=s", "grepcols=s", "test=s", "inclusive");
    die $usage unless exists $opts{net} && exists $opts{nodes} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});

    return %opts;
}
