#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(networkHashFromEdgeList);

# given a network and a list of proteins, pull out any edge that contains
# only proteins in the list

my %opts = getCommandLineOptions();

{
    my $netFile = $opts{net};
    my $nodeString = $opts{nodes};
    my $out = $opts{out};
    my $testString = $opts{test}; # if defined, only report edges between
    # these nodes and $opts{nodes}, along with edges among these

    my @nodes;
    {
	my %x = map { $_ => 1 } split ",", $nodeString;
	@nodes = sort keys %x;
    }
    my @testNodes;
    if (defined $testString) {
	my %x = map { $_ => 1 } split ",", $testString;
	@testNodes = sort keys %x;
    }
    my %net;
    networkHashFromEdgeList(\%net, $netFile, 0, 1, 1);
    open my $OUT, ">", $out or die "can't write to $out. $!";

    #if (@testNodes > 0) 
    for my $i (0..($#nodes-1)) {
	my $n1 = $nodes[$i];
	for my $j (($i+1)..$#nodes) {
	    my $n2 = $nodes[$j];
	    if (exists $net{$n1}{$n2}) {
		say $OUT join "\t", $n1, $n2, $net{$n1}{$n2}
	    }
	}
    }

}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -net net.in -nodes node1,node2,... -out output < ".
	" -test test1,test2 >\n";

    my %opts = ();
    GetOptions(\%opts, "net=s", "nodes=s", "out=s", "test=s");
    die $usage unless exists $opts{net} && exists $opts{nodes} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});

    return %opts;
}
