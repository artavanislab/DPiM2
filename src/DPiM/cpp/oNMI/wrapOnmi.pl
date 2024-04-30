#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(shuffle);
use File::Temp qw(tempfile);
use HomeBrew::IO qw(checkExist readColRef readList readHeader readColsRef);

# send correctly formatted files to onmi

#test4();
#exit;

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};
    my $minCluster = $opts{mincluster};
    my $refFile = $opts{ref};

    die "truthonly and clusterref cannot be used simultaneously"
	if exists $opts{clusterref} && $mode eq 'truthonly';
    
    my (%truth, $trueClusters);
    unless (exists $opts{clusterref}) {
	$trueClusters = clusterFromCommaCol(\%truth, $refFile, 'proteins');
    }

    
    my ($nClusters, %prediction);
    if ($mode eq 'mcode') {
	$nClusters = clusterFromNetworkList(\%prediction, $in, $minCluster);
	if (exists $opts{clusterref}) {
	    $trueClusters = clusterFromNetworkList(\%truth, $refFile, $minCluster);
	}
    } elsif ($mode eq 'mcl') {
	$nClusters = clusterFromEachLine(\%prediction, $in, $minCluster);
	if (exists $opts{clusterref}) {
	    $trueClusters = clusterFromEachLine(\%truth, $refFile, $minCluster);
	}
    } elsif ($mode eq 'truthonly') {
	$nClusters = $trueClusters;
	%prediction = %truth;
    }
    #my ($min, $max, $median) = clusterStats(\%prediction);
    
    say for NMI(\%prediction, \%truth, $minCluster);
}

exit;

sub getCommandLineOptions {

    my @modes = qw(mcode mcl truthonly);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	ref => '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOComplexes.min4',
	mincluster => 3,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;
    
    my $usage = "usage: $0 -in input < $modeString $defaultString ".
	" -clusterref >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "ref=s", "mincluster=i", 
	       "clusterref");
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

sub NMI {
    my ($X, $Y, $minCluster) = @_;

    my $Xinverse = invertCluster($X);
    my $Yinverse = invertCluster($Y);

    my ($XFile, $YFile, $OUT);
    ($OUT, $XFile) = tempfile();
    say $OUT join "\t", @$_ for values %$Xinverse;
    close $OUT;
    ($OUT, $YFile) = tempfile();
    say $OUT join "\t", @$_ for values %$Yinverse;
    close $OUT;

    my $oNMI = '/home/glocke/DPiM/cpp/oNMI/onmi';
    my $cmd = "$oNMI $XFile $YFile";
    #say $cmd;
    my @read = `$cmd`;
    my ($max, $lfk, $sum);
    for (@read) {
	my $line = $_;
	chomp;
	my @spl = split;
	if ($line =~ /Max/) {
	    $max = $spl[-1];
	} elsif ($line =~ /lfk/) {
	    $lfk = $spl[-1];
	} elsif ($line =~ /Sum/) {
	    $sum = $spl[-1];
	}
    }
    return ($max, $lfk, $sum);
}


# repeats the test in the paper
sub test1 {
    my (@true, @tests);

    my $nClust = 20;
    for my $i (0..($nClust-1)) {
	my $start = 10*$i;
	push @true, [$start..($start+9)];
	my @test;
	push @test, $_ for @true;
	push @tests, \@test;
    }
    my %true = listsToClusterFormat(\@true);
    my %universe = map {$_=> 1} 0..(10*$nClust-1);
    
    my @NMI;
    for my $testArray (@tests) {
	my %test = listsToClusterFormat($testArray);
	push @NMI, [NMI(\%test, \%true)];
    }
    dumpNMI(\@NMI);
}

# like the test in the paper, but the test case has an extra (unmatched) cluster
sub test2 {
    my (@true, @tests);

    my @badClust = ('a'..'j');
    
    my $nClust = 20;
    for my $i (0..($nClust-1)) {
	my $start = 10*$i;
	push @true, [$start..($start+9)];
	my @test = (\@badClust);
	push @test, $_ for @true;
	push @tests, \@test;
    }
    my %true = listsToClusterFormat(\@true);
    my %universe = map {$_=> 1} 0..(10*$nClust);
    $universe{$_} = 1 for @badClust;
    
    my @NMI;
    for my $testArray (@tests) {
	my %test = listsToClusterFormat($testArray);
	push @NMI, NMI(\%test, \%true, \%universe, 'max');
    }

    dumpNMI(\@NMI);

}

# same size of clusters, but increasing degree of shuffling
sub test3 {
    my (@true, @tests);

    my $nClust = 20;
    my $size = 10;
    @true = map { my $start = $_*$size; 
		  my $end = $start+$size-1;
		  [$start..$end] } 0..($nClust-1);

    # half the space is not in the true set
    my %universe = map {$_ => 1} 0..(2*$size*$nClust -1);

    my $nTest = 10;

    for my $i (1..$nTest) {
	my $fidelity = $i/$nTest;
	my $nCorrect = int($fidelity * $size);
	my %elems = %universe;

	my @test;	
	for my $j (0..($nClust-1)) {
	    my @clust;
	    # add correct elements
	    for my $k (0..($nCorrect-1)) {
		my $elem = $j*$size + $k;
		push @clust, $elem;
		delete $elems{$elem} if exists $elems{$elem};
	    }
	    # add incorrect elements
	    my @k = shuffle keys %elems;
	    while (@clust < $size) {
		my $elem = 0+ pop @k;
		push @clust, $elem;
	    }
	    push @test, \@clust;
	}
	push @tests, \@test;
    }

    my %true = listsToClusterFormat(\@true);

    my @NMI;
    for my $testArray (@tests) {
	my %test = listsToClusterFormat($testArray);
	push @NMI, [NMI(\%test, \%true)];
    }
    dumpNMI(\@NMI);
    #die Dumper(\$tests[-2]);
}
# a test with greater overlap
# truth puts qw(a b c) into every cluster
# plus sequential numbers unique to each cluster
sub test4 {
    my (@true, @tests);
   
    my $nClust = 20;
    my $sequen = 10;
    my @core = qw(a b c d e f g);
    @true = map { my $start = $_*$sequen; 
		  my $end = $start+$sequen-1;
		  [$start..$end] } 0..($nClust-1);
    unshift @$_, @core for @true;
    
    # half the space is not in the true set
    my %pickFrom = map {$_ => 1} 0..(2*$sequen*$nClust -1);

    my $size = $sequen + @core;
    my $nTest = $size;
    for my $i (1..$nTest) {
	my $fidelity = $i/$nTest;
	my $nCorrect = int($fidelity * $size);
	my %elems = %pickFrom;

	my @test;	
	for my $j (0..($nClust-1)) {
	    # add correct elements
	    my @thisTrue = @{ $true[$j] };
	    my @clust = @thisTrue[0..($nCorrect-1)];

	    # add incorrect elements
	    my @k = shuffle keys %elems;
	    while (@clust < $size) {
		my $elem = 0+ pop @k;
		push @clust, $elem;
	    }
	    push @test, \@clust;
	}
	push @tests, \@test;
    }

    my %true = listsToClusterFormat(\@true);

    my @NMI;
    for my $testArray (@tests) {
	my %test = listsToClusterFormat($testArray);
	push @NMI, [NMI(\%test, \%true)];
    }
    dumpNMI(\@NMI);    
}

# for debugging
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

sub dumpNMI {
    my $nmiArray = shift;
    say "NMI<Max>\tlfkNMI\tNMI<Sum>";
    for my $row (@$nmiArray) {
	printf("%.3f\t%.3f\t%.3f\n", @$row);
    }
}

