#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min max sum);
use Statistics::Basic qw(median);
use HomeBrew::IO qw(checkExist readColRef readList readHeader readColsRef);
use DpimLib qw(networkHashFromEdgeList);

# calculate the following confusion matrix (and related statistics)
#            True  False
# Positive  | +,+ | -,+ |
# Negative  | -,- | -,+ |
# 
# where + means "cocluster" and - means "don't cocluster"
# the first refers to truth, the second to prediction
# this confusion matrix is counted over every pair of proteins in the network

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
    my $netFile = $opts{net};

    my %universe;
    networkHashFromEdgeList(\%universe, $netFile, 'symmetric');
    
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
    
    my %conf = confusionMatrix(\%prediction, \%truth, \%universe);
    
    my $p = precision(\%conf);
    my $r = recall(\%conf);

    say $p; 
    say $r;
    say fStat($r, $p); ## f statistic == harmonic mean of $p and $r
    say MCC_formula(\%conf);
    if ($opts{verbose}) {
	my ($min, $max, $median) = clusterStats(\%prediction);
	say $nClusters;
	say 0+ keys %prediction; # number of proteins among all clusters
	say $min; # min size of cluster
	say $max; # max size of cluster
	say $median; # median ''     '''
	say $conf{true}{positive};
	say $conf{true}{negative};
	say $conf{false}{positive};
	say $conf{false}{negative};
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
	net => '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;
    
    my $usage = "usage: $0 -in input < $modeString $defaultString -verbose >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "ref=s", "net=i", "verbose");
    die $usage unless exists $opts{in};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{ref});
    checkExist('f', $opts{net});

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

# see head of document for definition
# implementation:
# for proteins in either $pred and $truth, do a complete count
# for proteins not in pred or truth, all pairs (internal and between outside 
#   and inside) are true negatives
sub confusionMatrix {
    my ($pred, $truth, $universe) = @_;

    my %conf = ( true => {positive=>0, negative=>0}, 
		false => {positive=>0, negative=>0} );
    
    my %union = map { $_ => 1 } keys %$pred, keys %$truth;
    my $outside = 0+ grep { ! exists $union{$_} } keys %$universe;
    $conf{true}{negative} = (($outside)*($outside-1)/2) + 
	$outside*(keys %union);

    #my @proteins = sort keys %union;
    my @proteins = keys %union;
    for my $i (0..($#proteins-1)) {
	my $p1 = $proteins[$i];
	for my $j (($i+1)..$#proteins) {
	    my $p2 = $proteins[$j];

	    my $posNegPred = coCluster($pred->{$p1}, $pred->{$p2});
	    my $posNegKey = ($posNegPred)?'positive':'negative';

	    my $posNegTruth = coCluster($truth->{$p1}, $truth->{$p2});
	    my $trueF = ($posNegPred == $posNegTruth);

	    #say "$p1, $p2" if ($trueF) && (!$posNegPred);

	    my $trueFKey = ($trueF)?'true':'false';
	    $conf{$trueFKey}{$posNegKey}++;
	}
    }

    return %conf;
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

    return 0;
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


# is the confusion matrix calculated correctly?
sub unitTest0 {    
    my @trueClusters = (
	[qw(0 1 2 3)],
	[qw(4 5 6 7)] );
    my @predClusters = (
	[qw(0 1 2 4)],
	[qw(3 5 6)],
	);
    my %universe = map {$_ => 1} 0..9;

    my %truthConf = (
	true => {
	    positive => (4*3/2)*2,
	    negative => 4*4 + 2*8 + 1,
	    # inter-cluster negatives + inside-outside + outside outside
	},
	false => { positive => 0, negative => 0 }
	);
    my %predConf = (
	true => {
	    positive => 3 + 1, # 0,1,2 and 5,6
	    negative => 2*8+1 + 3*3 + 2, # outside + (0,1,2) x (5,6,7) + 3x4,7,
	},
	false => {
	    positive => 3 + 2, # 4x0,1,2 and 3x5,6
	    negative => 3 + 3 + 2, # 3x0,1,2, 4x5,6,7, 7x5,6
	}
	);
    

    my %truth = listsToClusterFormat(\@trueClusters);
    my %pred = listsToClusterFormat(\@predClusters);

    #use Test::More tests=>2;
    my %conf = confusionMatrix(\%truth, \%truth, \%universe);
    #is_deeply(\%conf, \%truthConf);
    %conf = confusionMatrix(\%pred, \%truth, \%universe);
    #is_deeply(\%conf, \%predConf);


    exit;
}

sub unitTest1 {    

    my $nClust = 5;
    my $minClust = 4;
    my $maxClust = 8;
    my $diff = $maxClust - $minClust;
    
    my $cnt = 0;
    my @trueClusters = ();
    my $truePos = 0;
    my @sizes;
    while ($nClust > @trueClusters) {
	my $size = $minClust + int(rand($diff+1));
	push @sizes, $size;
	my $upBound = $cnt+$size;
	my @clust = ($cnt..($upBound-1));
	push @trueClusters, \@clust;
	$truePos += $size * ($size-1)/2;
	$cnt = $upBound;
    }
    my $outside = 10 * $nClust;
    my $maxID = $cnt + $outside-1;
    my %universe = map { $_ => 1 } 0..$maxID;
    my $nProt = $maxID+1;
    my $totalPairs = $nProt*($nProt-1)/2;
    my $trueNeg = $totalPairs - $truePos;
    my %trueConf = (
	true => { positive => $truePos, negative => $trueNeg },
	false => { positive => 0, negative => 0 }
	);
    
    my %true = listsToClusterFormat(\@trueClusters);
    
    #use Test::More tests=>3;
    {
	my %conf = confusionMatrix(\%true, \%true, \%universe);
	#is_deeply(\%conf, \%trueConf, "can I get the true x true confusion?" );
    }

    {
	# make every cluster correct except with one false positive and 
	# one false negative
	my @predClusters;
	
	my $top = $maxID;
	for my $c (@trueClusters) {
	    my @x = @$c;
	    pop @x;
	    push @x, $top;
	    push @predClusters, \@x;
	}
	my $tP = sum ( map { ($_-1)*($_-2)/2 } @sizes );
	my $fP = sum ( map { $_-1 } @sizes );
	my $fN = $fP;
	my $tN = $totalPairs - $tP - $fP - $fN;
	
	my %predConf = (
	    true =>  { positive => $tP, negative => $tN },
	    false => { positive => $fP, negative => $fN },
	);
	my %pred = listsToClusterFormat(\@predClusters);
	
	my %conf = confusionMatrix(\%pred, \%true, \%universe);
	#is_deeply(\%conf, \%predConf, "can I get -1+1 confusion?");

	{
	    my $p = precision(\%conf);
	    my $r = recall(\%conf);
	    my $f = fStat($r, $p);
	    say "MCC = ", MCC_formula(\%conf);
	    say "F = $f";
	}
	
	{
	    my $trueClust = $nClust;
	    my @rawPrecision = jaccardTest(\@predClusters, \@trueClusters);
	    my $p = sum(@rawPrecision)/$nClust;
	    
	    my @rawRecall = jaccardTest(\@trueClusters, \@predClusters);
	    my $r = sum(@rawRecall)/$trueClust;

	    my $f = fStat($r, $p);
	    say "jaccard f = $f";
	}
    }

    {
	# make every cluster correct except with one false positive and 
	# one false negative
	# omit one cluster entirely
	my @predClusters;
	
	my $top = $maxID;
	for my $c (@trueClusters) {
	    my @x = @$c;
	    pop @x;
	    push @x, $top;
	    push @predClusters, \@x;
	}
	pop @predClusters;
	my @allButLast = @sizes;
	my $last = pop @allButLast;
	my $tP = sum ( map { ($_-1)*($_-2)/2 } @allButLast );
	my $fP = sum ( map { $_-1 } @allButLast );
	my $fN = $fP + $last*($last-1)/2;
	my $tN = $totalPairs - $tP - $fP - $fN;
	
	my %predConf = (
	    true =>  { positive => $tP, negative => $tN },
	    false => { positive => $fP, negative => $fN },
	);
	my %pred = listsToClusterFormat(\@predClusters);
	
	my %conf = confusionMatrix(\%pred, \%true, \%universe);
	#is_deeply(\%conf, \%predConf, "can I get -1+1 missing one confusion?");

	{
	    my $p = precision(\%conf);
	    my $r = recall(\%conf);
	    my $f = fStat($r, $p);
	    say "MCC = ", MCC_formula(\%conf);
	    say "F = $f";
	}
	
	{
	    my $trueClust = $nClust;
	    my @rawPrecision = jaccardTest(\@predClusters, \@trueClusters);
	    my $p = sum(@rawPrecision)/$nClust;
	    
	    my @rawRecall = jaccardTest(\@trueClusters, \@predClusters);
	    my $r = sum(@rawRecall)/$trueClust;

	    my $f = fStat($r, $p);
	    say "jaccard f = $f";
	}
    }

    exit;
}

sub listsToClusterFormat {
    my ($lists) = @_;
    
    my %ret;
    for my $i (0..$#$lists) {
	my $list = $lists->[$i];
	for my $j (@$list) {
	    $ret{$j}{$i} = 1;
	}
    }
    return %ret;
}
