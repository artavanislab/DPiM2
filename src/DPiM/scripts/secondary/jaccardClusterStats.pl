#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min max sum);
use Statistics::Basic qw(median);
use HomeBrew::IO qw(checkExist readColRef readList readHeader readColsRef);
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
    my $seed = $opts{seed};
    my $refFile = $opts{ref};
    my $minJaccard = $opts{minjaccard};

    my %truth;
    my $trueClusters = clusterFromCommaCol(\%truth, $refFile, 'proteins');
    #die Dumper(\%truth);
    
    my ($nClusters, %prediction);
    if ($mode eq 'mcode') {
	$nClusters = clusterFromNetworkList(\%prediction, $in);
    } elsif ($mode eq 'mcl') {
	$nClusters = clusterFromEachLine(\%prediction, $in);
    } elsif ($mode eq 'truthonly') {
	$nClusters = $trueClusters;
	%prediction = %truth;
    }
    my ($min, $max, $median) = clusterStats(\%prediction);
    
    if ($opts{basic}) {
	say $nClusters;
	say 0+ keys %prediction;
	say $min;
	say $max;
	say $median;
	exit;
    }
    
    my @truth = invertCluster(\%truth); 
    my @pred = invertCluster(\%prediction);
    @pred = grep { @$_ >= 3 } @pred;

    #warn "pred\ttruth\n";
    my @rawPrecision = jaccardTest(\@pred, \@truth, $minJaccard);
    my $p = sum(@rawPrecision)/$nClusters;

    #warn "truth\tpred\n";
    my @rawRecall = jaccardTest(\@truth, \@pred, $minJaccard);
    my $r = sum(@rawRecall)/$trueClusters;

    say $p; 
    say $r;
    say fStat($r, $p); ## f statistic == harmonic mean of $p and $r
    say sum(@rawPrecision);
    if ($opts{verbose}) {
	say $nClusters;
	say 0+ keys %prediction;
	say $min;
	say $max;
	say $median;
	say sum(@rawPrecision);
	say sum(@rawRecall);
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my @modes = qw(mcode mcl truthonly);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	ref => '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOComplexes.min4',
	minjaccard => 0.5,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;
    
    my $usage = "usage: $0 -in putative.list < $modeString ".
	"$defaultString -verbose -basic >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "ref=s", "minjaccard=f", 
	"verbose", "basic");
    die $usage unless exists $opts{in};

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
    my ($ret, $in, $col) = @_;
    
    my @read;
    readColRef(\@read, $in, 'proteins', "\t");
    for my $i (0..$#read) {
	my @cluster = split /,/, $read[$i];
	$ret->{$_}{$i} = 1 for @cluster;
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
# return a list of lists: @ret = ([cluster1_1, cluster1_2,...], [cluster2_1...])
sub invertCluster {
    my ($hash) = @_;

    my @array; # return value
    my %clusterMap;
    my $clusterCount = 0;
    for my $protein (keys %$hash) {
	for my $cl (keys %{$hash->{$protein}}) {
	    if (! exists $clusterMap{$cl}) {
		$clusterMap{$cl} = $clusterCount++;
		$array[$clusterMap{$cl}] = [];
	    }
	    push @{$array[$clusterMap{$cl}]}, $protein;
	}
    }

    for my $cl (@array) {
	$cl = [ sort @$cl ];
    }
    
    return @array;
}

sub jaccardTest {
    my ($test, $against, $min) = @_;

    my @ret;
    for my $cl1 (@$test) {
	my $flag = 0;
	for my $cl2 (@$against) {
	    my $j = jaccard($cl1, $cl2);
	    if ($j >= $min) {
		$flag = 1;
	    }
	    #warn "\t".(join(",", @$cl1))."\t".(join(",", @$cl2))."\t$j\n" 
	#	if $j != 0;
	}
	push @ret, $flag;
    }

    return @ret;
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
