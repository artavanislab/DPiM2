#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(networkHashFromEdgeList);

# look for nodes present in target but not in source (X) and vice versa (Y)
# for every x in X, find the y in Y with max(jacc(x,y))
#   where jacc(x,y) is the jaccard score for overlap of neighbors

my %opts = getCommandLineOptions();

{
    my $sourceFile = $opts{source};
    my $targetFile = $opts{target};
    my $out = $opts{out};

    my (%source, %target);
    networkHashFromEdgeList(\%source, $sourceFile, 'sym');
    networkHashFromEdgeList(\%target, $targetFile, 'sym');

    my (@sNoT, @tNoS);
    for my $node (keys %source) {
	push @sNoT, $node if ! exists $target{$node};
    }
    for my $node (keys %target) {
	push @tNoS, $node if ! exists $source{$node};
    }

    my (%match, %jacc);
    for my $s (@sNoT) {
	my ($bestMatch, $bestJ, $bestI) = ('none', 0, 0);
	for my $t (@tNoS) {
	    my ($j, $intersect) = jaccard($source{$s}, $target{$t});
	    $jacc{$s}{$t} = $j;
	    if ($j > $bestJ) {
		$bestMatch = $t;
		$bestJ = $j;
		$bestI = $intersect;
	    }
	}
	$match{$s} = { target  => $bestMatch, 
		       jaccard => sprintf("%.4f", $bestJ), 
		       degree  => 0+ keys %{$source{$s}},
		       intersect => $bestI
	};
    }

    my @order = sort {$match{$b}{degree} <=> $match{$a}{degree}} @sNoT;
    
    
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 tried to find matches for nodes present in $sourceFile not present $targetFile among the converse set";
    my @cols = qw(source target jaccard degree intersect);
    say $OUT join "\t", @cols;
    shift @cols;
    for my $s (@order) {
	say $OUT join "\t", $s, map {$match{$s}{$_}} @cols;
    }
    close $OUT; 
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0  -source doThese.net -target appearHere?net -out ".
	"targetList.tsv\n";

    my %opts = ();
    GetOptions(\%opts, "source=s", "target=s", "out=s");
    die $usage unless exists $opts{source} && exists $opts{target} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{source});
    checkExist('f', $opts{target});

    return %opts;
}

sub jaccard {
    my ($set1, $set2) = @_;

    my @k1 = keys %$set1;
    my @k2 = keys %$set2;
    my $debug = @k2 == 1 && $k2[0] eq 'FBgn0035471';
    $debug = undef;

    my ($union, $intersect) = (@k1 + @k2, 0);
    for my $k (@k1) {
	if (exists $set2->{$k}) {
	    say "intersect $k" if $debug;
	    $intersect++;
	    $union--;
	} else {
	    say "no-union  $k" if $debug;
	}
    }
    die "$intersect/$union" if $debug;
    return (0+ $intersect/$union, $intersect);
}
