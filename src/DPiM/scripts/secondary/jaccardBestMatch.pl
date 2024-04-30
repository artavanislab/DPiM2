#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min max sum);
use Statistics::Basic qw(median);
use HomeBrew::IO qw(checkExist readColRef readList readHeader readColsRef writeCols);
#use DpimLib qw(getLineAPMS);

# Calculate the jaccard distance between every putative cluster and every true 
#  cluster.  For each putative cluster, count a match when jaccard >= 0.5
#  Recall is the fraction of true clusters with a match
#  Precision is the fraction of putative clusters with a match

my %opts = getCommandLineOptions();

# cluster format: 
#   clust{$fbgn} = { $cluster1 => 1, $cluster2, => 1... }
#     where $clusterX is a numeric identifier for the cluster
# the values are in "member format"
{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};
    my $refFile = $opts{ref};

    my %truth;
    my (%metaData, @cols);
    my $trueClusters = clusterFromCommaCol(\%truth, $refFile, 'proteins', 
					   \%metaData, \@cols);
    my ($nClusters, %prediction);
    if ($mode eq 'mcode') {
	$nClusters = clusterFromNetworkList(\%prediction, $in);
    } elsif ($mode eq 'mcl') {
	$nClusters = clusterFromEachLine(\%prediction, $in);
    } elsif ($mode eq 'truthonly') {
	$nClusters = $trueClusters;
	%prediction = %truth;
    }
    
    my %truthInv = %{ invertCluster(\%truth)};
    my %predInv = %{ invertCluster(\%prediction)};
    {
	my @k = keys %predInv;
	@k = grep { @{ $predInv{$_} } < 3 } @k;
	delete $predInv{$_} for @k;
    }
    
    for my $i (keys %truthInv) {
	#die Dumper($truthInv[$i], $predInv[0]);
	my ($j, $jacc) = bestMatch($truthInv{$i}, \%predInv);
	$metaData{$i}{match} = join ",", @{$predInv{$j}};
	$metaData{$i}{jaccard} = sprintf("%.3f", $jacc);
    }

    my @k = sort {$metaData{$b}{jaccard} <=> $metaData{$a}{jaccard}} 
        keys %metaData;
    push @cols, 'jaccard', 'match';

    my @outData;
    for my $c (@cols) {
	push @outData, [map { $metaData{$_}{$c} } @k];
    }
    my $preComments = "# found best matches for primary clusters defined in $refFile among comparison clusters in $in";
    my $header = join "\t", @cols;
    my $format = join "\t", ('%s') x (0+ @cols);
    writeCols($out, \@outData, $header, $format, $preComments);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my @modes = qw(mcl mcode truthonly);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	ref => '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOComplexes.min4',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;
    
    my $usage = "usage: $0 -in putative.list -out output < $modeString ".
	"$defaultString -verbose >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "ref=s", "minjaccard=f", 
	"verbose");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{ref});

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    return %opts;
}

# $ret is a ref to hash in cluster format
sub clusterFromCommaCol {
    my ($ret, $in, $col, $metaData, $cols) = @_;
    
    my @read;
    readColRef(\@read, $in, 'proteins', "\t");
    for my $i (0..$#read) {
	my @cluster = split /,/, $read[$i];
	$ret->{$_}{$i} = 1 for @cluster;
    }

    if (defined $metaData) {
	push @$cols, readHeader($in);
	my @meta;
	readColsRef(\@meta, $in, $cols, undef, "\t");
	$metaData->{$_} = $meta[$_] for 0..$#meta;
    }
    
    return 0+ @read;
}

# assume first two columns of files are node1, node2
sub clusterFromNetworkList {
    my ($ret, $listFile) = @_;

    my @files = readList($listFile);
    checkExist('f', $_) for @files;
    
    for my $i (0..$#files) {
	my @header = readHeader($files[$i]);
	pop @header while @header > 2;
	my @read;
	readColsRef(\@read, $files[$i], \@header);
	for my $row (@read) {
	    $ret->{$row->{$header[0]}}{$i}=1;
	    $ret->{$row->{$header[1]}}{$i}=1;
	}
    }

    return 0+ @files;
}


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

sub clusterStats {
    my ($clus) = @_;
    
    my %clusCount; # clusCount{$id} = # of members
    for my $c (values %$clus) {
	for my $clusID (keys %$c) {
	    $clusCount{$clusID}++;
	}
    }
    my @cnts = values %clusCount;
    return (min(@cnts), max(@cnts), median(@cnts))
}

  ##
 # Jaccard-similarity based statistics
## 

# takes a cluster format as input
# return a hashref of lists: $ret->{$clusterID} = [member1, member2, ...]
sub invertCluster {
    my ($hash) = @_;

    my %inverted; # return value
    for my $protein (keys %$hash) {
	for my $cl (keys %{$hash->{$protein}}) {
	    $inverted{$cl} //= [];
	    push @{$inverted{$cl}}, $protein;
	}
    }

    for my $cl (values %inverted) {
	$cl = [ sort @$cl ];
    }
    
    return \%inverted;
}

# find the best match for needle in haystack
sub bestMatch {
    my ($needle, $haystack) = @_;

    my ($bestK, $bestJac) = (-1, -1);
    for my $k (keys %$haystack) {
	my $jac = jaccard($needle, $haystack->{$k});
	if ($jac > $bestJac) {
	    $bestK = $k;
	    $bestJac = $jac;
	}
    }
    #die Dumper($needle, $haystack->{$bestI]) if $bestJac > 0.8;
    return ($bestK, $bestJac);
}

# intersection of elements divided by union of elements
sub jaccard {
    my ($set1, $set2) = @_;
    
    my %set2 = map { $_ => 1 } @$set2;
    my %union = map {$_ => 1} @$set1, @$set2;
    my %intersect;
    for my $k (@$set1) {
	$intersect{$k} = 1 if exists $set2{$k}
    }

    return (0+ keys %intersect) / (0+ keys %union);
}

## f statistic := harmonic mean of $p and $r
sub fStat {
    my ($r, $p) = @_;
    if (($r+$p) == 0) {
	return "nan";
    }
    return 2*$r*$p/($r+$p); 
}
