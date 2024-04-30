#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(sum);

use lib '/home/glocke/perl5/lib/perl5';
use Math::GSL::CDF qw(gsl_cdf_hypergeometric_P);

use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist readColsRef writeCols);

# given so many rejected edges (black) and so many accepted edges (white)
# what is the probability to find a given number of rejected/accepted edges
# on a given day under a hypergeomtric expectation?


my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    my @batches;
    readColsRef(\@batches, $in, [qw(enrich total date)]);

    my $whiteTot = 0;
    my $allTot = 0;
    my %days;
    my @inCols = qw(enrich total);
    for my $row (@batches) {
	my $date = $row->{date};
	$days{$date} //= { map {$_ => 0} @inCols };
	$days{$date}{$_} += $row->{$_} for @inCols;
	$whiteTot += $row->{enrich};
	$allTot += $row->{total};
    }

    if (exists($opts{month})) {
	%days = aggregateMonths(\%days);
    }
    
    my $blackTot = $allTot - $whiteTot;
    for my $d (keys %days) {
	my $acc = $days{$d}{accepted} = $days{$d}{enrich};
	my $rej = $days{$d}{rejected} = $days{$d}{total} - $days{$d}{accepted};
	$days{$d}{pbin} = 
	    gsl_cdf_hypergeometric_P($acc, $whiteTot, $blackTot, $acc+$rej);
    }

    
    my @colNames = qw(accepted rejected pbin);
    my @data;
    my @k = sort keys %days;
    push @data, \@k;
    for my $c (@colNames) {
	push @data, [map {$days{$_}{$c}} @k]; 
    }
    unshift @colNames, "date";
    my $header = (join "\t", @colNames)."\n";
    my $format = (join "\t", qw(%s %d %d %.4e))."\n";
    my $preComments = "# binomial probabilities assessed for a total of $whiteTot enriched edges and $blackTot non-enriched edges";
    writeCols($out, \@data, $header, $format, $preComments);
    
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

    my $usage = "usage: $0 -in enrich.tab -out output < -month >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "month");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

sub aggregateMonths {
    my ($byDay) = @_;

    my %byMonth; # return value
    my @cols = qw(enrich total);
    for my $d (keys %$byDay) {
	my $first = $d;
	$first =~ s/\d\d$/01/;
	if (exists $byMonth{$first}) {
	    $byMonth{$first}{$_}+= $byDay->{$d}{$_} for @cols;
	} else {
	    $byMonth{$first} = $byDay->{$d};
	}
    }

    return %byMonth;
}
