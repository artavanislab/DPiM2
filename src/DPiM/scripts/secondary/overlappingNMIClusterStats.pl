#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min max sum shuffle);
use Statistics::Basic qw(median);
use HomeBrew::IO qw(checkExist readColRef readList readHeader readColsRef writeCols);
use DpimLib qw(networkHashFromEdgeList);

# calculate the normalized mutual information between known complexes and
#  putative clusters *given that the known complexes overlap*
# The formulae are (from McDaid, Greene, and Hurley 2013) more or less  
#   straightforward, but too lengthy to write in this head.

# reference:
# McDaid, Greene, and Hurley - Normalized Mutual Information to evaluate 
#   overlapping community finding algorithms
#   http://arxiv.org/abs/1110.2515

my $LOG2 = 1/log(2);
#test1();
#test2();
#test3();
#test4();
#exit;
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
    my $normMode= $opts{norm};

    die "truthonly and clusterref cannot be used simultaneously"
	if exists $opts{clusterref} && $mode eq 'truthonly';
    
    my %universe;
    networkHashFromEdgeList(\%universe, $netFile, 'symmetric');

    my (%truth, $trueClusters);
    unless (exists $opts{clusterref}) {
	$trueClusters = clusterFromCommaCol(\%truth, $refFile, 'proteins');
    }
    
    my ($nClusters, %prediction);
    if ($mode eq 'mcode') {
	$nClusters = clusterFromNetworkList(\%prediction, $in);
	if (exists $opts{clusterref}) {
	    $trueClusters = clusterFromNetworkList(\%truth, $refFile);
	}
    } elsif ($mode eq 'mcl') {
	$nClusters = clusterFromEachLine(\%prediction, $in);
	if (exists $opts{clusterref}) {
	    $trueClusters = clusterFromEachLine(\%truth, $refFile);
	}
    } elsif ($mode eq 'truthonly') {
	$nClusters = $trueClusters;
	%prediction = %truth;
    }
    #my ($min, $max, $median) = clusterStats(\%prediction);
    
    say NMI(\%prediction, \%truth, \%universe, $normMode);
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

    my @norms = qw(avg joint max);
    my %norms = map {$_ => 1} @norms;
    @arr = @norms;
    $arr[0] = "*".$arr[0]."*";
    my $normString = "-norm ".join("/", @arr);

    my %defaults = (
	ref => '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOComplexes.min4',
	net => '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;
    
    my $usage = "usage: $0 -in input < $modeString $defaultString $normString ".
	" -clusterref -printmatches file.out>\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "ref=s", "net=i", "norm=s",
	"clusterref", "printmatches=s");
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

    die "you must select one of the following normalizations: ", 
        join(", ",keys(%norms)), "\n" 
	if (exists $opts{norm} && ! exists $norms{$opts{norm}});
    $opts{norm} //= $norms[0];

    return %opts;
}

# $ret is a ref to hash in cluster format
sub clusterFromCommaCol {
    my ($ret, $in, $col) = @_;
    
    my @read;
    my $idCol = (readHeader($in))[0];
    readColsRef(\@read, $in, [$idCol, 'proteins'], undef, "\t");
    for my $row (@read) {
	my $id = $row->{$idCol};
	my @cluster = split /,/, $row->{proteins};
	$ret->{$_}{$id} = 1 for @cluster;
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

  ##
 # Normalized Mutual Information 
##

# find the normalized mutual information between two sets of clusters X and Y
# $X and $Y are in cluster format
# $n is the number of nodes in the whole network
# if $maxFlag is set, use max(H(X), H(Y)) instead of H(X,Y)
sub NMI {
    my ($X, $Y, $universe, $normMode) = @_;
    
    #my $n = 0+ keys %$universe;

    my $ignorance = 'ignorance'; # unknown cluster
    my %partialUniverse;
    #if ($ignoranceMode) {
    if (0) {
	addIgnorance($X, $universe, $ignorance);
	addIgnorance($Y, $universe, $ignorance);
	%partialUniverse = %$universe;
    } else {
	# find the subset of proteins that are found in X and/or Y
	%partialUniverse = map { $_ => 1 } (keys %$X, keys %$Y);
    }
    my $n = 0+ keys %partialUniverse;
    #say "size of universe vs partial universe = ", (0+ keys %$universe)
    #	, " vs $n";
    
    #die Dumper([sort {$a <=> $b} keys %partialUniverse]);
    
    my $Xinverse = invertCluster($X);
    my $Yinverse = invertCluster($Y);
    my $Xids = [keys %$Xinverse]; # cluster id's
    my $Yids = [keys %$Yinverse]; #

    my ($a, $b, $c, $d) = firstStats($X, $Y, $Xids, $Yids, 
				     \%partialUniverse);
    my ($H_XiYj, $H_YjXi) = ({}, {}); # equation 2
    computeClusterH($H_XiYj, $H_YjXi, $Xids, $Yids, $a, $b, $c, $d, $n);
    if (exists $opts{printmatches}) {
	printMatches($H_XiYj, $H_YjXi, $Xinverse, $Yinverse, $Xids, $Yids,
		     $n, $opts{printmatches});
    }
    
    my ($H_XY, $H_YX) = HXlessY($H_XiYj, $H_YjXi, $Xids, $Yids); # eqn 4

    my $H_X = selfEntropy($Xinverse, $n);
    my $H_Y = selfEntropy($Yinverse, $n);

    my $MI= 0.5 * ($H_X - $H_XY + $H_Y - $H_YX);
    #die "MI= 0.5 * ($H_X - $H_XY + $H_Y - $H_YX)";

    #say "H(X)   =? I(X:Y)+ H(X|Y) => $H_X =? ", ($MI + $H_XY);
    #say "H(Y)   =? I(X:Y)+ H(Y|X) => $H_Y =? ", ($MI + $H_YX);
    #say "H(X,Y) =? ", ($H_X + $H_YX), " = H(X) + H(Y|X)";
    #say "H(X,Y) =? ", ($H_Y + $H_XY), " = H(Y) + H(Y|X)";
    #say "H(X,Y) =? ", $H_XY + $MI + $H_YX, " = H(X|Y) + I(X:Y) + H(Y|X)";
    #say "H_X, H_Y, H_XY, H_YX, MI = ", join ", ", map {sprintf("%.1f", $_)} $H_X, $H_Y, $H_XY, $H_YX, $MI;
    if ($normMode eq 'max') {
	#say "denom = ", max($H_X, $H_Y);
	return $MI / max($H_X, $H_Y);
    } elsif ($normMode eq 'joint') {
	# no need to average the three computations for H(X,Y) because the value
	# of I(X:Y) is already computed based on the same average
	# H(X,Y) = H(X) + H(Y|X) = H(Y) + H(X|Y) = H(X|Y) + I(X:Y) + H(Y|X)
	
	my $totalH = ( ($H_XY + $MI + $H_YX) +
		       ($H_X + $H_YX) + ($H_Y + $H_XY) ) / 3;
	#say "denom = $totalH";
	return $MI / $totalH;
    } else {
	#say "denom = ", 0.5*($H_X + $H_Y);
	return 2 * $MI/($H_X + $H_Y);
    }
}

# put any fbgn without a cluster into the ignorance cluster
sub addIgnorance {
    my ($C, $universe, $ignorance) = @_;

    my $noCluster = { $ignorance => 1 };
    for my $p (keys %$universe) {
	#$noCluster;
	$C->{$p} //= {"$p$ignorance" => 1}
    }

    return;
}

# as defined prior to eqn 1 in McDaid, Greene, and Hurley
# i refers to a cluster in X, j refers to a cluster in Y
# a_ij = Sum{m in universe} ( m not in X_i && m not in Y_j)
# b_ij = Sum{m in universe} ( m not in X_i && m     in Y_j)
# c_ij = Sum{m in universe} ( m     in X_i && m not in Y_j)
# d_ij = Sum{m in universe} ( m     in X_i && m     in Y_j)
# return format for the above is a hash ref
# $x->{$i}{$j} := x_ij
sub firstStats {
    my ($X, $Y, $Xids, $Yids, $universe) = @_;

    my (%a, %b, %c, %d);
    for my $i (@$Xids) {
	for my $j (@$Yids) {
	    $a{$i}{$j} = $b{$i}{$j} = $c{$i}{$j} = $d{$i}{$j} = 0;
	}
    }
    for my $p (keys %$universe) {
	my $x = $X->{$p};
	my $y = $Y->{$p};
	for my $i (@$Xids) {
	    for my $j (@$Yids) {
		if (exists $x->{$i}) {
		    if (exists $y->{$j}) {
			# d_ij = ( m     in X_i && m     in Y_j)
			$d{$i}{$j}++;
		    } else {
			# c_ij = ( m     in X_i && m not in Y_j)
			$c{$i}{$j}++;
		    }
		} else {
		    if (exists $y->{$j}) {
			# b_ij = ( m not in X_i && m     in Y_j)
			$b{$i}{$j}++;
		    } else {
			# a_ij = ( m not in X_i && m not in Y_j)
			$a{$i}{$j}++;
		    }		    
		}
	    }
	}
    }

    return (\%a, \%b, \%c, \%d);
}

# the lack of information between two vectors
# equation 2, written as H*(X_i, Y_j) in the paper (H star)
sub computeClusterH {
    my ($H_XiYj, $H_YjXi, $Xids, $Yids, $a, $b, $c, $d, $n) = @_;

    for my $i (@$Xids) {
	for my $j (@$Yids) {
	    my $A = $a->{$i}{$j};
	    my $B = $b->{$i}{$j};
	    my $C = $c->{$i}{$j};
	    my $D = $d->{$i}{$j};
	    my $hA = h($A, $n);
	    my $hB = h($B, $n);
	    my $hC = h($C, $n);
	    my $hD = h($D, $n);
	    if ($hA + $hD >= $hB + $hC) {
		$H_XiYj->{$i}{$j} = $hA + $hB + $hC + $hD - 
		    h($B+$D, $n) - h($A+$C, $n);
		$H_YjXi->{$i}{$j} = $hA + $hB + $hC + $hD - 
		    h($C+$D, $n) - h($A+$B, $n);
	    } else {
		$H_XiYj->{$i}{$j} = h($A+$B, $n) + h($C+$D, $n);
		$H_YjXi->{$i}{$j} = h($A+$C, $n) + h($B+$D, $n);
	    }
	}
    }
    return;
}

# h(w,n) = -wlog_2(w/n)
sub h {
    my ($w, $n) = @_;
    return 0 if $w == 0;
    return -$w * $LOG2 * log($w/$n);
}

# compute equation 4 for H(X|Y) and H(Y|X)
# return ($H_XY, $H_YX) := ( H(X|Y), H(Y|X) )
sub HXlessY {
    my ($H_XiYj, $H_YjXi, $Xids, $Yids) = @_;

    # first compute equation 3
    my (%HXiY, %HYjX);
    for my $i (@$Xids) {
	for my $j (@$Yids) {
	    my $hij = $H_XiYj->{$i}{$j};
	    $HXiY{$i} //= $hij;
	    $HXiY{$i} = $hij if $hij < $HXiY{$i};
	}
    }
    for my $j (@$Yids) {
	for my $i (@$Xids) {
	    my $hij = $H_YjXi->{$i}{$j};
	    $HYjX{$j} //= $hij;
	    $HYjX{$j} = $hij if $hij < $HYjX{$j};
	}
    }
    #say "H_XiY{$_} = ", sprintf("%.3f", $HXiY{$_}) for sort {$a <=> $b} @$Xids;
    #say "H_YjX{$_} = ", sprintf("%.3f", $HYjX{$_}) for sort {$a <=> $b} @$Yids;

    my ($H_XY, $H_YX) = (0, 0);
    for my $i (@$Xids) {
	$H_XY += $HXiY{$i};
    }
    for my $j (@$Yids) {
	$H_YX += $HYjX{$j};
    }

    return ($H_XY, $H_YX);
}

# unconditional entropy for a cover
sub selfEntropy {
    my ($inverse, $n) = @_;

    my $H = 0;
    for my $cluster (values %$inverse) {
	my $size = 0+ @$cluster;
	$H += h($size, $n) + h($n-$size, $n);
    }
    
    return $H;
}

# report which clusters match which other clusters
sub printMatches {
    my ($H_XiYj, $H_YjXi, $Xinverse, $Yinverse, $Xids, $Yids, $n, $baseOut)
	= @_;

    my ($first, $last);
    {
	my @spl = split/\./, $baseOut;
	$last = pop @spl;
	$first = join ".", @spl;
    }

    my @cols = qw(from to H missingH nmi membersFrom membersTo);
    my $header = join "\t", @cols;
    my $format = join "\t", qw(%s %s %.3f %.3f %.3f %s %s);
	
    {
	my @byRow;
	for my $i (@$Xids) {
	    # find the Y cluster with the least missing information
	    my $bestJ = $Yids->[0];
	    my $bestH;
	    for my $j (@$Yids) {
		my $hij = $H_XiYj->{$i}{$j};
		$bestH //= $hij;
		if ($hij < $bestH) {
		    $bestJ = $j;
		    $bestH = $hij;
		}
	    }
	    my %row = (from => $i, to=> $bestJ, missingH => $bestH);
	    my $size = 0+ @{ $Xinverse->{$i} };
	    $row{H} = h($size, $n) + h($n - $size, $n);
	    $row{nmi} = ($row{H} - $bestH) / $row{H};
	    $row{membersFrom} = join ",", @{ $Xinverse->{$i} };
	    $row{membersTo} = join ",", @{ $Yinverse->{$bestJ} };
	    push @byRow, \%row;
	}
	my @data;
	for my $c (@cols) {
	    push @data, [map {$_->{$c}} @byRow];
	}
	my $out = "$first.predToRef.$last";
	writeCols($out, \@data, $header, $format);
    }

    {
	my @byRow;
	for my $j (@$Yids) {
	    # find the Y cluster with the least missing information
	    my $bestI = $Xids->[0];
	    my $bestH;
	    for my $i (@$Xids) {
		my $hij = $H_YjXi->{$i}{$j};
		$bestH //= $hij;
		if ($hij < $bestH) {
		    $bestI = $i;
		    $bestH = $hij;
		}
	    }
	    my %row = (from => $j, to=> $bestI, missingH => $bestH);
	    my $size = 0+ @{ $Yinverse->{$j} };
	    $row{H} = h($size, $n) + h($n - $size, $n);
	    $row{nmi} = ($row{H} - $bestH) / $row{H};
	    $row{membersFrom} = join ",", @{ $Yinverse->{$j} };
	    $row{membersTo} = join ",", @{ $Xinverse->{$bestI} };
	    push @byRow, \%row;
	}
	my @data;
	for my $c (@cols) {
	    push @data, [map {$_->{$c}} @byRow];
	}
	my $out = "$first.refToPred.$last";
	writeCols($out, \@data, $header, $format);
    }
}

# # #       # # #
## ## DEBUG ## ##
# # # X # X # # #


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
	push @NMI, NMI(\%test, \%true, \%universe, 'max');
    }
    say Dumper('test1', \@NMI);
    #die Dumper(\@NMI);
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

    say Dumper('test2', \@NMI);

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

    say "avg\tjoint\tmax";
    for my $testArray (@tests) {
	my %test = listsToClusterFormat($testArray);
	printf("%.3f\t%.3f\t%.3f\n", 
	       NMI(\%test, \%true, \%universe, 'avg'),
	       NMI(\%test, \%true, \%universe, 'joint'),
	       NMI(\%test, \%true, \%universe, 'max')
	    );
    }
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

    say "avg\tjoint\tmax";
    for my $testArray (@tests) {
	my %test = listsToClusterFormat($testArray);
	my %universe = map { $_ => 1 } (keys %test, keys %true);
	printf("%.3f\t%.3f\t%.3f\n", 
	       NMI(\%test, \%true, \%universe, 'avg'),
	       NMI(\%test, \%true, \%universe, 'joint'),
	       NMI(\%test, \%true, \%universe, 'max')
	    );
    }
    #die Dumper(\$tests[-2]);
    
}

