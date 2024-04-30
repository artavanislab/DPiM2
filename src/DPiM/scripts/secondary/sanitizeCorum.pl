#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min max);
use Statistics::Basic qw(median stddev);
use Graph;
use HomeBrew::IO qw(checkExist readColsRef writeCols);
use DpimLib qw(networkHashFromEdgeList);

# remove complexs less than 4 members
# remove duplicate complexes
# collect similar complexes into bigger complexes

my %opts = getCommandLineOptions();

{
    my $netFile = $opts{net};
    my $out = $opts{out};
    my $corumFile = $opts{corum};
    my $maxJaccard = $opts{maxjaccard};
    my $minComplex = $opts{mincomplex};
    my $subsetCutoff = $opts{subsetcutoff};

    my %net;
    networkHashFromEdgeList(\%net, $netFile, 'symmetric');

    my @cols = qw(id name proteins);
    my @complexes;
    readColsRef(\@complexes, $corumFile, \@cols, undef, "\t");
    say "found ", (0+ @complexes), " complexes";

    $_->{complex} = [split ',', $_->{proteins} // die Dumper($_)] 
	for @complexes;
    $_->{complex} = pruneComplex($_->{complex}, \%net) for @complexes;
    $_->{size} = 0+ @{$_->{complex}} for @complexes;
    @complexes = grep { $_->{size} >= $minComplex } @complexes;

    say "", (0+ @complexes), " complexes after requiring footpriqnt";
    
    # remove complexes that are similar to one another
    
    # first remove those that are contained in each other
    @complexes = removeSubsets(\@complexes, $subsetCutoff);
    say "", (0+ @complexes), " complexes after removing subsets";

    @complexes = collectSimilar(\@complexes, $maxJaccard, $opts{modjaccard});
    say "", (0+ @complexes), " complexes after removing similar complexes";

    my @data = (
	[ map { $_->{id} } @complexes],
	[ map { $_->{name} } @complexes],
	[ map { join ",", @{ $_->{complex} } } @complexes],
	);
    my $header = join "\t", qw(id name proteins);
    my $format = join "\t", qw(%s %s %s);
    my $preComments = "# projected $corumFile onto $netFile\n";
    $preComments   .= "# requiring at least $minComplex members\n";
    $preComments   .= "# requiring a maximum Jaccard similarity of $maxJaccard\n";
    writeCols($out, \@data, $header, $format, $preComments);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	corum=>'/home/glocke/DPiM/corum/complexList_allOrthoMin3.tab',
	maxjaccard => 0.5, # ~0.5% of pairs match at this level
	mincomplex => 4,
	subsetcutoff => 0.8,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -net network -out output < ".
	"$defaultString -modjaccard >\n";

    my %opts = ();
    GetOptions(\%opts, "net=s", "out=s", "corum=s", "maxjaccard=f", 
	       "mincomplex=i", "subsetcutoff=f", "modjaccard");
    die $usage unless exists $opts{net} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});

    return %opts;
}

# remove all elements that are not present in $net
sub pruneComplex {
    my ($complex, $net) = @_;

    my @ret;
    for my $member (@$complex) {
	push @ret, $member if exists $net->{$member};
    }

    return \@ret;
}

# remove complexes that are contained in other complexes
sub removeSubsets {
    my ($complexes, $subsetCutoff) = @_;

    my %subsets; # corum id's for complexes contained in other complexes
    for my $c1 (@$complexes) {
	for my $c2 (@$complexes) {
	    next unless $c1->{size} < $c2->{size};
	    my $jac = modifiedJaccard($c1->{complex}, $c2->{complex});
	    die "can't get modJac ", Dumper($c1, $c2) if ! defined $jac;
	    if ($jac > $subsetCutoff)
	    {
		$subsets{$c1->{id}} = 1;
		last;
	    }
	}
    }

    return grep {! exists $subsets{$_->{id}} } @$complexes;
}

# intersection of elements divided by size of smaller set
sub modifiedJaccard {
    my ($set1, $set2) = @_;

    my $denom = min ((0+ @$set1), (0+ @$set2));
    
    my %set2 = map { $_ => 1 } @$set2;
    my %intersect;
    for my $k (@$set1) {
	$intersect{$k} = 1 if exists $set2{$k}
    }

    return (0+ keys %intersect) / $denom;
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



# find jaccard similarity between all complexes
# draw an edge between any two complexes with jaccard > maxJaccard
# collect members of all connected components into a super-complex
# return the result
sub collectSimilar {
    my ($complexes, $maxJaccard, $modJaccard) = @_;

    my $similarity = \&jaccard;
    $similarity = \&modifiedJaccard if defined $modJaccard;

    my @jaccard;
    my %matches;
    for my $i (0..($#$complexes - 1)) {
	for my $j (($i+1)..$#$complexes) {
	    my $jac = $similarity->($complexes->[$i]{complex}, 
				    $complexes->[$j]{complex});
	    $jaccard[$i][$j] = $jaccard[$j][$i] = $jac;
	    if ($jac > $maxJaccard) {
		$matches{$i}{$j} = $matches{$j}{$i} = $jac;
		#say "$i\t$j";
	    }
	}
	$jaccard[$i][$i] = 1;
    }
    $jaccard[$#$complexes][$#$complexes] = 1;
    #warn((join ",", @{ $jaccard[$_] }), "\n") for 0..$#jaccard;
    #exit;

    my %drop = map {$_ => 1} keys %matches;

    my $g = Graph->new( undirected => 1 );

    for my $src ( keys %matches ) {
	for my $tgt ( keys %{ $matches{$src} } ) {
	    $g->add_edge($src, $tgt);
	}
    }
    
    my @connected = $g->connected_components; # node lists
    @connected = sort { @$a <=> @$b } @connected;
    
    my @newComplexes;
    for my $comp (@connected) {
	push @newComplexes, collect(map {$complexes->[$_]} @$comp);
    }
    say "max new cluster = ", max(map {$_->{size}} @newComplexes);
	    
    
    # note that a cluster is defined by a group of complexes where each member
    #  of the cluster is jaccard-similar to at least one other complex
    #my %drop = map { $_ => 1 } findBiggestInCluster(\@complexes, \%matches, 
#						    \@jaccard);
    my @keep = grep { ! exists $drop{$_} } 0..$#$complexes;
    
    my @ret = @$complexes[@keep];
    push @ret, @newComplexes;
    return @ret;

}

sub collect {
    my @complexes = @_;

    my @order = sort 0..$#complexes;
    @complexes = @complexes[@order];
    my %members;
    for my $c (@complexes) {
	$members{$_} = 1 for @{$c->{complex}};
    }

    my %ret;
    $ret{id}   = join "_", (map { $_->{id} }   @complexes);
$ret{name} = join "_", (map { $_->{name} } @complexes);
    $ret{complex} = [ sort keys %members ],
    $ret{size} = 0+ keys %members;

    return \%ret;
}

# make a graph where two complexes are connected if they are jaccard similar
# find the connected components of this graph
# for each component, find the complex with the most members
# return a list of those complexes that are not the biggest in their component
sub findBiggestInCluster {
    my ($complexes, $matches, $jaccard) = @_;
    
    my $g = Graph->new( undirected => 1 );
    $g->add_vertex($_) for 0..$#$complexes;

    for my $src ( keys %$matches ) {
	for my $tgt ( keys %{ $matches->{$src} } ) {
	    $g->add_edge($src, $tgt);
	}
    }
    
    my @subgraphs = $g->connected_components; # node lists
    @subgraphs = sort { @$a <=> @$b } @subgraphs;
    for my $sg (@subgraphs) {
	say 0+ @$sg if 1 < 0+ @$sg;
	warn "$_\n" for @$sg;
    }
    exit;

    my (@jMatch, @jNoMatch); # collect the jaccard similarities for cluster
    # members that match and don't match
    # note that jMatch > $maxJaccard and jNoMatch <= $maxJaccard by construction
    
    my @nonMaxMembers; # these are the complexes to be removed, 

    
    for my $component (@subgraphs) {
	my $biggestSize = 0;
	my $biggestComplex = undef;
	for my $i (0..$#$component) {
	    my $c1 = $component->[$i];
	    my $size = $complexes->[$c1]{size};
	    if ($size > $biggestSize) {
		$biggestSize = $size;
		$biggestComplex = $c1;
	    }
	    last if $i == $#$component;
	    for my $j (($i+1)..$#$component) {
		my $c2 = $component->[$j];
		if (exists $matches->{$c1}{$c2}) {
		    push @jMatch, $jaccard->[$c1][$c2];
		} else {
		    push @jNoMatch, $jaccard->[$c1][$c2];
		}
	    }
	}
	push @nonMaxMembers, grep { $_ != $biggestComplex } @$component;
    }

    say "median, sd of jaccard for intra-cluster = "
	, (median(@jMatch, @jNoMatch)), ", ", (stddev(@jMatch, @jNoMatch));
    say "median, sd of jaccard for matching within cluster = "
	, (median(@jMatch)), ", ", (stddev(@jMatch));
    say "median, sd of jaccard for non-matching within cluster = "
	, (median(@jNoMatch)), ", ", (stddev(@jNoMatch));

    return @nonMaxMembers;
}