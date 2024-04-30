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

# Use Try2Step to find optimal parameters for clustering algorithms

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

    my %truth;
    my $trueClusters = clusterFromCommaCol(\%truth, $refFile, 'proteins');

    
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
    
    if ($opts{jaccard}) {
	my @truth = invertCluster(\%truth); 
	my @pred = invertCluster(\%prediction);
	@pred = grep { @$_ >= 3 } @pred;

	my @rawPrecision = jaccardTest(\@pred, \@truth);
	my $p = sum(@rawPrecision)/$nClusters;

	my @rawRecall = jaccardTest(\@truth, \@pred);
	my $r = sum(@rawRecall)/$trueClusters;

	say $p; 
	say $r;
	say fStat($r, $p); ## f statistic == harmonic mean of $p and $r
	say $nClusters;
	say 0+ keys %prediction;
	say $min;
	say $max;
	say $median;
	say sum(@rawPrecision);
	say sum(@rawRecall);
	exit;
    }
    

    my %conf = confusionMatrix(\%prediction, \%truth);
    
    my $p = precision(\%conf);
    my $r = recall(\%conf);

    say $p; 
    say $r;
    say fStat($r, $p); ## f statistic == harmonic mean of $p and $r
    say $nClusters;
    say 0+ keys %prediction; # number of proteins among all clusters
    say $min; # min size of cluster
    say $max; # max size of cluster
    say $median; # median ''     '''
    say MCC_formula(\%conf);
    say $conf{true}{positive};
    say $conf{true}{negative};
    say $conf{false}{positive};
    say $conf{false}{negative};
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
	seed => 'system_time',
	ntry => 20,
	nstep => 100,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;
    
    my $usage = "usage: $0 -in input -out output < $modeString ".
	"$defaultString -jaccard >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "ref=s", "seed=i", "ntry=i", 
	       "nstep=i", "jaccard");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

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

# MCC := (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
# this is calculated over all pairs of proteins, where "positive" means these 
#   proteins cluster together and "negative" means they don't
# $pred and $truth arguments are in cluster format
sub MCC {
    my ($pred, $truth) = @_;

    my %conf = confusionMatrix($pred, $truth);
    return MCC_formula(\%conf);
}

# MCC := (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
sub MCC_formula {
    my ($conf) = @_;
    
    my $num = ($conf->{true}{positive}*$conf->{true}{negative} - 
	       $conf->{false}{positive}*$conf->{false}{negative});
    my $denom = sqrt(
	($conf->{true}{positive} + $conf->{false}{positive}) *
	($conf->{true}{positive} + $conf->{false}{negative}) *
	($conf->{true}{negative} + $conf->{false}{positive}) *
	($conf->{true}{negative} + $conf->{false}{negative}) );
    return $num / $denom ;
}

sub confusionMatrix {
    my ($pred, $truth) = @_;

    my %ret = ( true => {positive=>0, negative=>0}, 
		false => {positive=>0, negative=>0} );
    my @k = keys %$pred;
    for my $i (0..($#k-1)) {
	my $p1 = $k[$i];
	for my $j (($i+1)..$#k) {
	    my $p2 = $k[$j];
	    say "$p1,$p2";

	    my $posNeg = coCluster($pred->{$p1}, $pred->{$p2});
	    my $posNegKey = ($posNeg)?'positive':'negative';

	    # trueF is true if p1 and p2 coCluster in both pred and truth
	    # trueF is true if p1 and p2 do not coCluster in both
	    # false otherwise
	    my $trueF = coCluster($truth->{$p1}, $truth->{$p2});
	    # right now, trueF is merely telling us whether they coCluster in 
	    #   truth.  This is the correct value if they coCluster in pred.
	    #   If they do not coCluster in pred, then coCluster(truth) is 
	    #   the opposite of the correct result, so just invert it.
	    $trueF = !$trueF if ! $posNeg;
	    my $trueFKey = ($trueF)?'true':'false';
	    $ret{$trueFKey}{$posNegKey}++;
	}
    }
    exit;
    return %ret;
}

# return 1 iff the arguments are in the same cluster
# return undef iff they are not
# note that when there are overlapping clusters, the arguments could be in 
#   different clusters and the same cluster, in which case coCluster returns 1
#
# arguments are in member format
sub coCluster {
    my ($protein1, $protein2) = @_;

    for my $k (keys %$protein1) {
	return 1 if exists $protein2->{$k};
    }

    return undef;
}

  ##
 # the following classifier statistics take a confusion matrix as input
##
sub precision {
    my ($conf) = @_;
    
    return $conf->{true}{positive} /
	($conf->{true}{positive} + $conf->{false}{positive});
}

sub recall {
    my ($conf) = @_;
    
    return $conf->{true}{positive} /
	($conf->{true}{positive} + $conf->{false}{negative});
}

##
 #  end confusion statistics
  ##

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
    my ($test, $against) = @_;

    my @ret;
    for my $cl1 (@$test) {
	my $flag = 0;
	for my $cl2 (@$against) {
	    if (jaccard($cl1, $cl2) > 0.5) {
		$flag = 1;
	    }
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
    return 2*$r*$p/($r+$p); 
}
