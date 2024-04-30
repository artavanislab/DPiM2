#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);
use DpimLib qw(networkHashFromEdgeList);

# for each network find
#   number of edges
#   number of nodes

my %opts = getCommandLineOptions();

{
    my $netListFile = $opts{in};
    my $out = $opts{out};
    my $nameString = $opts{names};
    
    my @netFiles = readList($netListFile);
    my @cols;
    if (defined $nameString) {
	@cols = split /,/, $nameString;
	die "different number of names (".(0+ @cols).") vs. networks (".
	    (0+ @cols).")" if @cols != @netFiles;
    } else {
	@cols = map { "net$_" } 1..(0+ @netFiles);
    }


    my @stats; # stats[i] = {edges, nodes};
    for my $f (@netFiles) {
	my %network;
	networkHashFromEdgeList(\%network, $f);
	push @stats, netStats(\%network);
    }

    my @rows = qw(nodes edges);

    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# found edge statistics for ".(join (",", @netFiles));
    say $OUT join "\t", '', @cols;
    for my $r (@rows) {
	say $OUT join "\t", $r, map { $_->{$r} } @stats;
    }
    close $OUT;
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

    my $usage = "usage: $0 -in net.list -out output < -names net1,net2,...\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "names=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

sub netStats {
    my ($net) = @_;

    my %nodes;
    my $edgeCnt = 0;
    for my $n1 (keys %$net) {
	$nodes{$n1} = 1;
	for my $n2 (keys %{ $net->{$n1} }) {
	    $edgeCnt++;
	    $nodes{$n2} = 1;
	}
    }

    return { nodes=> (0+ keys %nodes), edges=> $edgeCnt };
}
