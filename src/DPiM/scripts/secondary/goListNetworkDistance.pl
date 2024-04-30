#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColRef writeCols);
use DpimLib qw(networkHashFromEdgeList);

# find the HG-score of all members of a complex

my %opts = getCommandLineOptions();

{
    #my $in = $opts{in};
    my $distFile = $opts{dist};
    my $listFile = $opts{list};
    my $out = $opts{out};
    my $mode = $opts{mode};

    my @clusters;
    listsFromCommaCol(\@clusters, $listFile, 'proteins');
    my %network;
    networkHashFromEdgeList(\%network, $distFile, 'symmetric', 'score');

    if ($mode eq 'hgscore') {
	my @distances;
	my @connected;
	my $cutoff = 0;
	
	for my $c (@clusters) {
	    push @connected, findDistances($c, \%network, \@distances, $cutoff);
	}
	my $whichConnect = sum(@connected);

	open my $OUT, ">", $out or die "can't read from $out. $!";
	say $OUT "# $whichConnect out of ", (0+ @connected), " clusters are connected at cutoff = $cutoff";
	say $OUT $_ for @distances;
	close $OUT;
    } elsif ($mode eq 'inout') {
	my @memberEdges;
	my @extraEdges;
	my @fraction;
	for my $c (@clusters) {
	    next if @$c < 4;
	    my ($m, $e) = countEdges($c, \%network);
	    push @memberEdges, $m;
	    push @extraEdges, $e;
	    push @fraction, $m/($m+$e);
	}

	my $inTot = sum(@memberEdges);
	my $outTot = sum(@extraEdges);
	
	my $header = join "\t", qw(inward outward fraction);
	my $format = join "\t", qw(%d %d %.3f);
	my $preComments= "# $inTot inward edges vs $outTot outward edges ( ".
	    sprintf("%.1f", 100*$inTot/($inTot+$outTot)). "% )\n";
	writeCols($out, [\@memberEdges, \@extraEdges, \@fraction], 
		  $header, $format, $preComments);
    }
    
}
exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    my @modes = qw(hgscore inout);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	dist => '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.noCutoff.network',
	list => '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOComplexes.min4',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out raw.list < $modeString $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "dist=s", "list=s", "mode=s");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    return %opts;
}

# $ret is a ref to hash in cluster format
sub listsFromCommaCol {
    my ($ret, $in, $col) = @_;
   
    my @read;
    readColRef(\@read, $in, 'proteins', "\t");
    for my $i (0..$#read) {
	my @cluster = split /,/, $read[$i];
	push @$ret, \@cluster;
    }
}

# find the HGScore among all pairs within a set of nodes 
sub findDistances {
    my ($nodes, $network, $distances, $cutoff) = @_;
    $cutoff //= 0;

    my $allConnected = 1;
    for my $i (0..$#$nodes) {
	my $connected = undef;
	my $n1 = $nodes->[$i];
	for my $j (($i+1)..$#$nodes) {
	    my $n2 = $nodes->[$j];
	    my $dist = $network->{$n1}{$n2} // -1;
	    push @$distances, $dist;
	    $connected = 1 if $dist >= $cutoff;
	}
	$allConnected = $allConnected && $connected;
    }

    return $allConnected // 0;
}

# count how many connections club members have to other members and to 
#   non-members
sub countEdges {
    my ($club, $network) = @_;
    
    #                 ####
    my ($in, $out) = (0, 0);
    #                \\()//
    #                 |  |
    #                 -Edvard Munch

    my %clubHash = map { $_ => 1 } @$club;
    for my $mem (@$club) {
	for my $k (keys %{ $network->{$mem} }) {
	    if (exists $clubHash{$k}) {
		$in++;
	    } else {
		$out++;
	    }
	}
    }

    return($in, $out);
}
