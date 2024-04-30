#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(networkHashFromEdgeList);

# find the fraction of edges within a network that connect clusters

my %opts = getCommandLineOptions();

# cluster format: 
#   clust{$fbgn} = { $cluster1 => 1, $cluster2, => 1... }
#     where $clusterX is a numeric identifier for the cluster
# the values are in "member format"

{
    my $netFile = $opts{net};
    my $clustFile = $opts{clust};
    my $minCluster = 2;
    
    my %clusters;
    my $nClusters = clusterFromEachLine(\%clusters, $clustFile, $minCluster);

    my %net;
    networkHashFromEdgeList(\%net, $netFile);
    my $total=0;
    my $intra=0;
    for my $n1 (keys %net) {
	for my $n2 (keys %{ $net{$n1} }) {
	    $total++;
	    next unless exists $clusters{$n1} && exists $clusters{$n2};
	    for my $clusterID (keys %{ $clusters{$n1} }) {
		if (exists $clusters{$n2}{$clusterID}) {
		    $intra++;
		    last;
		}
	    }
	}
    }
    say $intra/$total;
    say $intra;
    say $total;
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

    my $usage = "usage: $0 -net in.net -clust mcl.txt \n";

    my %opts = ();
    GetOptions(\%opts, "net=s", "clust=s");
    die $usage unless exists $opts{net} && exists $opts{clust};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});
    checkExist('f', $opts{clust});

    return %opts;
}

# $ret is a ref to hash in cluster format
# assume each line is a cluster, with members in a tab-delimited list
sub clusterFromEachLine {
    my ($ret, $in, $minCluster) = @_;
    $minCluster //= 3;

    open my $IN, "<", $in or die "can't read from $in. $!";
    my $i=0;
    while (<$IN>) {
	next if /^#/;
	chomp;
	my @cluster = split /\t/;
	next if @cluster < $minCluster;
	$ret->{$_}{$i} = 1 for @cluster;
	$i++;
    }

    return $i;
}
