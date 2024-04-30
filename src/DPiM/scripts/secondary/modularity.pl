#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(sum);
use HomeBrew::IO qw(checkExist readCols readList readHeader);
use DpimLib qw(networkHashFromEdgeList);

# apply the modularity measure defined in 
#   A. Lázár, D. Ábel, and T. Vicsek, European Physics Letters 2010
#   http://arxiv.org/pdf/0910.5072v1.pdf

#test();
#exit;

my %opts = getCommandLineOptions();

{
    my $netFile = $opts{net};
    my $moduleFile = $opts{mod};
    my $out = $opts{out};
    my $mode = $opts{mode};
    my $minCluster = $opts{mincluster};
    
    my @modules;
    my @moduleData;
    if ($mode eq 'reflist') {
	@moduleData = readCols($moduleFile, ['proteins'], 'line', "\t");
	@modules = map { [ split /,/, $_->{proteins} ] } @moduleData;
    } elsif ($mode eq 'perline') {
	my @rawList = readList($moduleFile);
	@modules = map { [ split /,/ ] } @rawList;
	$moduleData[$_]{line} = $rawList[$_] for 0..$#rawList;
    } elsif ($mode eq 'mcl') {
	readMCL(\@modules, $moduleFile, $minCluster);
	if ($mode eq 'iter') {
	    shift @$_ for @modules;
	}
	for my $i (0..$#modules) {
	    my $size = 0+ @{ $modules[$i] };
	    $moduleData[$i]{line} = 'clust'.($i+1)."_size$size";
	}
    } elsif($mode eq 'iter') {
	readMCL(\@modules, $moduleFile, $minCluster);
	my @names = map { shift @$_ } @modules;

	for my $i (0..$#modules) {
	    my $size = 0+ @{ $modules[$i] };
	    $moduleData[$i]{line} = 'clust'.$names[$i]."_size$size";
	}
    } else {
	die "unknown mode '$mode'";
    }
    my %net;
    networkHashFromEdgeList(\%net, $netFile, 'symmetric', undef, undef, $opts{human});
    #die Dumper(\@modules);
    
    my @modularity = modularity(\@modules, \%net);
    my @dense = map {density($_, \%net)} @modules;
    my $meanMod = sum(@modularity) / @modularity;
    say $meanMod;

    my $meanConnected; # calculate the mean only among clusters that have
    #                  # at least one edge outside the cluster
    # this requires that the connected component for this cluster contain
    # at least two clusters
    {
	my ($sum, $n) = (0,0);
	for my $i (0..$#modularity) {
	    if ($dense[$i] != $modularity[$i]) {
		$n++;
		$sum+=$modularity[$i];
	    }
	}
	if ($n == 0) {
	    $meanConnected = -666;
	} else {
	    $meanConnected = $sum / $n;
	}
    }
    
    my @indices = sort {$modularity[$b] <=> $modularity[$a] } 0..$#modularity;
    if ($opts{nosort}) {
	@indices = 0..$#modularity;
    }
    
    my @cols = readHeader($moduleFile);
    @cols = ('cluster') if $mode eq 'perline' || $mode eq 'mcl' || 
	$mode eq 'iter';
    push @cols, 'modularity', 'density';

    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# $0 computed modularity using $netFile and $moduleFile";
    say $OUT "# mean modularity $meanMod";
    say $OUT "# mean connected modularity $meanConnected";
    say $OUT join "\t", @cols;
    for my $i (@indices) {
	if (! defined $moduleData[$i]{line}) {
	    $moduleData[$i]{line} = join "\t", -1, "null", "FBgn0000000";
	    $moduleData[$i]{line} = 'null' if $mode eq 'perline';
	    $moduleData[$i]{line} = 'clustNull_size'.(0+ @{ $modules[$i] });
	}
	chomp $moduleData[$i]{line};
	say $OUT join "\t", $moduleData[$i]{line}, $modularity[$i], $dense[$i];
    }
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    # mode:
    my @modes = qw(reflist perline mcl iter);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	mincluster => 3,
	);
    my $defaultString = join " ", map { "-$_ $defaults{$_}" } keys %defaults;
    my $usage = "usage: $0 -net input -mod module.list -out output ".
	"< $modeString $defaultString -nosort -human >\n";

    my %opts = ();
    GetOptions(\%opts, "net=s", "mod=s", "out=s", "mode=s", "mincluster=i", 
	       "nosort", "human");
    die $usage unless exists $opts{net} && exists $opts{mod} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});
    checkExist('f', $opts{mod});

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    return %opts;
}

sub readMCL {
    my ($ret, $in, $minCluster) = @_;
    open my $IN, "<", $in or die "can't read from $in. $!";

    while (<$IN>) {
	next if /^#/;
	next if /^protein/;
	chomp;
	my $line = $_;
	my @cluster = split /\t/;
	next if @cluster < $minCluster;
	push @$ret, \@cluster;
    }
    return;
}

sub modularity {
    my ($modules, $net) = @_;

    # lump all non-annotated nodes into Russell's Paradox
    my $null = nullModule($modules, $net);
    push @$modules, $null if @$null > 1;
    
    # degree of each node
    my %degree = map { $_ => 0+ keys %{$net->{$_}} } keys %$net;
    ##die Dumper(\%degree);
    
    # how many modules does each node belong to?
    my %membership;
    calcMembership(\%membership, $modules, $net);
    
    my @mod; # modularity of each module
    for my $i (0..$#$modules) {
	push @mod, oneModule($modules->[$i], $i, $net, \%degree, \%membership);
    }
    #@mod = sort { $a<=> $b} @mod;
    #print Dumper(\@mod);
    #my $netMod = sum(@mod) / @mod;

    return @mod;
}

# make a list of every node not in a module
sub nullModule {
    my ($modules, $net, $membership) = @_;

    my %already;
    for my $mod (@$modules) {
	$already{$_} = 1 for @$mod;
    }

    my @nullMod = ();
    for my $node (keys %$net) {
	push @nullMod, $node if ! exists $already{$node};
    }

    return \@nullMod;
}

# for each node, make a hash of all modules it's in
# reverseModules{$node}{$module} = 1 if $node is in $module 
#   (where $module denotes an array index in the @modules list)
sub calcMembership {
    my ($ret, $modules) = @_;
    
    for my $i (0..$#$modules) {
	my $mod = $modules->[$i];
	for my $node (@$mod) {
	    $ret->{$node}{$i} = 1;
	}
    }

    return;
}

# find the modularity of this module
sub oneModule {
    my ($module, $modIndex, $net, $degree, $membership) = @_;

    my $ret=0;
    for my $member (@$module) {
	my $thisNode = 0;
	for my $neighbor (keys %{ $net->{$member} }) {
	    if (exists $membership->{$neighbor}{$modIndex} ) {
		$thisNode++;
	    } else {
		$thisNode--;
	    }
	}
	if (! exists $degree->{$member}) {
	    warn "can't find degree->{$member}\n";
	    next;
	}
	if (0 == (0+ keys %{ $membership->{$member} })) {
	    die Dumper($member, $membership->{$member});
	}
	$thisNode /= ($degree->{$member} * 
		      (0+ keys %{ $membership->{$member} }));
	#say "\t$member\t$thisNode";
	$ret+= $thisNode;
    }
    $ret *= density($module, $net) / @$module;

    return $ret;
}

sub density {
    my ($module, $net) = @_;

    # density is the number of edges among cluster members divided by the 
    #   number of possible edges
    my $edgeCount=0;
    for my $i (0..($#$module-1)) { # don't double-count
	my $nodeI = $module->[$i];
	for my $nodeJ (@$module[($i+1)..$#$module]) {
	    $edgeCount++ if exists $net->{$nodeI}{$nodeJ};
	}
    }
    
    my $maxEdges = @$module * $#$module / 2;

    die "maxEdges = 0\n".Dumper($module) if $maxEdges == 0;
    
    return $edgeCount / $maxEdges;
}

sub test {

    # two disconnected powersets
    my $size1 = 4;
    my $size2 = 5;

    my $prevNodes = 0;
    my (@correct, @incomplete, @swap, %net);

    @correct = ([], []);
    @incomplete = ([], []);
    @swap = ([], []);
    for my $i (1..$size1) {
	my $nodeI = $i + $prevNodes;
	#push @{$correct[0]}, $nodeI;
	for my $j (1..$size1) {
	    my $nodeJ = $j + $prevNodes;
	    next if $i == $j;
	    $net{$nodeI}{$nodeJ} = 1;
	}
    }

    $prevNodes = $size1;
    for my $i (1..$size2) {
	my $nodeI = $i + $prevNodes;
	for my $j (1..$size2) {
	    my $nodeJ = $j + $prevNodes;
	    next if $i == $j;
	    $net{$nodeI}{$nodeJ} = 1;
	}
    }

    
    @correct = ([1..$size1], [($size1+1)..($size1+$size2)]);
    my @modules = @correct;
    say "correct modularity is ", modularity(\@modules, \%net);

    # swap
    ($modules[0][0], $modules[1][0]) = ($modules[1][0], $modules[0][0]);
    say "swap-one modularity is ", modularity(\@modules, \%net);

    # drop one
    shift @{ $modules[0] };
    shift @{ $modules[1] };
    say "drop-one modularity is ", modularity(\@modules, \%net);
}

