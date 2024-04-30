#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);
use DpimLib qw(networkHashFromEdgeList);

# find the minimum score for all edges for all interactions

my %opts = getCommandLineOptions();

{
    my $inList = $opts{in};
    my $out = $opts{out};

    my @memberFiles = readList($inList);
    checkExist('f', $_) for @memberFiles;

    my %minScore; # minScore{fb1}{fb2} = min score for this this edge
    networkHashFromEdgeList(\%minScore, shift @memberFiles, undef, 'keepScore', 
			    'sort');
    #
    for my $f (@memberFiles) {
	say $f;
	my %net;
	networkHashFromEdgeList(\%net, $f, undef, 'keepScore', 'sort');
	my @delete; # delete[i] = [n1, n2] or n1
	# if $f doesn't contain n1, n2, delete it
	# if $f doesn't contain have any edges for n1, delete all
	for my $n1 (keys %minScore) {
	    if (! exists $net{$n1}) {
		push @delete, $n1;
		next;
	    }
	    for my $n2 (keys %{ $minScore{$n1} }) {
		if (! exists $net{$n1}{$n2}) {
		    push @delete, [$n1, $n2];
		} elsif ($minScore{$n1}{$n2} > $net{$n1}{$n2}) {
		    $minScore{$n1}{$n2} = $net{$n1}{$n2};
		}
	    }
	}
	for my $del (@delete) {
	    if (ref($del)) {
		my ($n1, $n2) = @$del;
		delete $minScore{$n1}{$n2};
	        delete $minScore{$n1} if 0 == keys %{ $minScore{$n1} };
	    } else {
		delete $minScore{$del};
	    }
	}
    }

    my @sort = sortByScore(\%minScore);

    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0";
    say $OUT "# found min scores among files listed in $inList";
    say $OUT join "\t", qw(protein1 protein2 score);
    for my $pair (@sort) {
	my ($n1, $n2) = @$pair;
	say $OUT join "\t", $n1, $n2, $minScore{$n1}{$n2};
    }
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

    my $usage = "usage: $0 -in real.list -out out.net\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

# ret = [[n1, n2], [n1, n2],...] 
#   where $net{$ret[i][0]}{$ret[i][1]} >= $net{$ret[i+1][0]}{$ret[i+1][1]}
sub sortByScore {
    my $net = shift;

    my %sortable;
    for my $n1 (keys %$net) {
	for my $n2 (keys %{ $net->{$n1} }) {
	    my $key = "$n1-$n2";
	    $sortable{$key} = $net->{$n1}{$n2};
	}
    }
    
    my @doubleKeys = sort {$sortable{$b} <=> $sortable{$a}} sort keys %sortable;

    return map { [split /-/, $_] } @doubleKeys;
}
