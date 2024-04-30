#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(networkHashFromEdgeList);

# given a support network, find the support such that there are N edges

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $n = $opts{n};
    #my $out = $opts{out};

    my %net;
    #say "read";
    networkHashFromEdgeList(\%net, $in, undef, 'keepScore');
    my @sort = sortByScore(\%net);

    my $edgeN = $sort[$n-1] // die "edge $n undefined (".(0+ @sort)." edges total)";
    my $support = $net{$edgeN->[0]}{$edgeN->[1]};
    my ($prevN, $prevSupp) = findBoundary($n, $support, \%net, \@sort, -1);
    my ($nextN, $nextSupp) = findBoundary($n, $support, \%net, \@sort, 1);

    say join "\t", $prevN, $prevSupp;
    say join "\t", $n, $support;
    say join "\t", $nextN, $nextSupp;
    #open my $OUT, ">", $out or die "can't write to $out. $!";
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

    my $usage = "usage: $0 -in input -n nEdges\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "n=i");
    die $usage unless exists $opts{in} && exists $opts{n};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

# ret = [[n1, n2], [n1, n2],...] 
#   where $net{$ret[i][0]}{$ret[i][1]} >= $net{$ret[i+1][0]}{$ret[i+1][1]}
sub sortByScore {
    my $net = shift;

    my %sortable;
    for my $n1 (keys %$net) {
	for my $n2 (keys %{ $net->{$n1} }) {
	    my $key = "$n1-$n2";
	    $sortable{$key} = $net->{$n1}{$n2};
	}
    }
    
    my @doubleKeys = sort {$sortable{$b} <=> $sortable{$a}} sort keys %sortable;

    return map { [split /-/, $_] } @doubleKeys;
}

# find the next edge with support ne to $support
# if direction is < 0, search downards
sub findBoundary {
    my ($n, $support, $net, $edges, $direction) = @_;

    $direction //= 1;
    my $pip = ($direction<1)?(-1):1;
    
    my ($thisSupport, $edge);

    do {
	$n+=$pip;
	$edge = $edges->[$n-1];
	$thisSupport = $net->{$edge->[0]}{$edge->[1]};
    } while ($n > 1 && $n < @$edges && $thisSupport == $support);
    return($n, $thisSupport);
}
