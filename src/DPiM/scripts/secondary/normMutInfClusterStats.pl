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

# calculate the normalized mutual information between known complexes and
#  putative clusters
# Define T = {t_1, t_2, ...} the set of true complexes
# Define P = {p_1, p_2, ...} the set of putative clusters
# x_i (x=t,p) is a set of genes
# I is the mutual information, H is the information entropy, 
#   N is the number of proteins
# NMI := 2* I/(H(T) + H(P))
# I := Sum{T, P}( P(t_i ^ p_j) log2(P(t_i ^ p_j) / (P(t_i) P(p_i))
#   where probabilities are calcluated empirically as |set|/N
# Note that when the intersection is zero, the sense of the formula is 
#   preserved by setting the contribution from i, j to zero
# H is calculated in the usual way
# H := - Sum(P(x_k) * log2 P(x_k)

test2();

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

    my $totalSize;
    {
	my %universe;
	networkHashFromEdgeList(\%universe, $netFile, 'symmetric');
	$totalSize = 0+ keys %universe;
    }

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
    #my ($min, $max, $median) = clusterStats(\%prediction);
    
    my @prediction = invertCluster(\%prediction);
    my @truth = invertCluster(\%truth);
    #@truth = ([1,2], [3,4,5,6]);
    #@prediction = ([3,4,5,6], [1,2]);
    say NMI(\@prediction, \@truth, $totalSize);
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
    
    my $usage = "usage: $0 -in input < $modeString $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "ref=s", "net=i");
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

  ##
 # Normalized Mutual Information 
##

# find the normalized mutual information between two sets of clusters C1 and C2
sub NMI {
    my ($C1, $C2, $N) = @_;

    # find one-point probabilities and entropies
    my ($H1, %P1);
    my $log2 = 1/log(2); # strictly speaking, this isn't necessary because
    #                  # eventually it cancels out.  Ok, fine, but do it anyway.
    for my $k (0..$#$C1) {
	my $nk = 0+ @{$C1->[$k]};
	my $p = $nk/$N;
	$P1{$k} = $p;
	$H1 -= $p * $log2 * log($p);
    }
    my ($H2, %P2);
    for my $k (0..$#$C2) {
	my $nk = 0+ @{$C2->[$k]};
	my $p = $nk/$N;
	$P2{$k} = $p;
	$H2 -= $p * $log2 * log($p);
    }
    
    # find mutual information (to be normalized)
    my $MI = 0;
    for my $k1 (0..$#$C1) {
	my $p1 = $P1{$k1};
	for my $k2 (0..$#$C2) {
	    my $p12 = intersection($C1->[$k1], $C2->[$k2]) / $N;
	    next if $p12 == 0;
	    my $p2 = $P2{$k2};
	    $MI += $p12 * $log2 * log($p12/($p1*$p2));
	}
    }
    say "2 * $MI / ($H1 + $H2);";
    return 2 * $MI / ($H1 + $H2);
}

# set intersection of two arrays
sub intersection {
    my ($set1, $set2) = @_;
    
    my %set2 = map { $_ => 1 } @$set2;
    my %intersect;
    for my $k (@$set1) {
	$intersect{$k} = 1 if exists $set2{$k}
    }

    return 0+ keys %intersect;
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


# like the test in the paper, but the test case has an extra (unmatched) cluster
sub test2 {
    my (@true, @tests);

    my @badClust = ('a'..'j');
    
    my $nClust = 20;
    for my $i (1..$nClust) {
	my $start = 10*$i;
	push @true, [$start..($start+9)];
	my @test = (\@badClust);
	push @test, $_ for @true;
	push @tests, \@test;
    }
    my $totalSize = 10*($nClust+1);
    
    my @NMI;
    for my $test (@tests) {
	push @NMI, NMI($test, \@true, $totalSize);
    }

    die Dumper(\@NMI);

}
