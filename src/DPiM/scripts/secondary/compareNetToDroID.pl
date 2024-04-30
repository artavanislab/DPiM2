#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readCols readHeader);
use DpimLib qw(networkHashFromEdgeList);

# for all edges with score at least X, find the fraction of edges that have
#  support in the DroID set

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $droidFile = $opts{droid};
    
    my %net;
    networkHashFromEdgeList(\%net, $in, undef, 'keepScore');
    my @edges; # edges[i] = {n1, n2, score};
    edges(\@edges, \%net);
    @edges = sort {$b->{score} <=> $a->{score} } @edges;

    my %droid;
    networkHashFromEdgeList(\%droid, $droidFile, 'symm', 'keepScore');

    my %droidSupport; # droidSupport{score} = { 1=>sum1, 2=> sum2, 3=>sum3...}
    # where sumX is the number of edges up to and including i that have X
    #   support in droid
    my ($sum1, $sum2, $sum3, $sum4) = (0,0,0,0);
    my $n = 0;
    for my $e (@edges) {
	my ($n1, $n2, $score) = ($e->{n1}, $e->{n2}, $e->{score});
	my $support = $droid{$n1}{$n2} // 0;

	if ($support > 3) {
	    $sum1++;
	    $sum2++;
	    $sum3++;
	    $sum4++;
	} elsif ($support > 2) {
	    $sum1++;
	    $sum2++;
	    $sum3++;
	} elsif ($support > 1) {
	    $sum1++;
	    $sum2++;
	} elsif ($support > 0) {
	    $sum1++;
	}
	$n++;
	$droidSupport{$score} = { 1=> $sum1, 2=> $sum2, 3=>$sum3, 4=>$sum4, 
				  n=> $n};
    }

    my $scoreCol = (readHeader($in))[2];
    $scoreCol =~ s/HGSC/HGSc/; ## some file says HGSCore
    
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# found droid support of $in according to $droidFile";
    my @cols = (1..4);
    say $OUT join "\t", $scoreCol, qw(frac1 frac2 frac3 frac4 sum1 sum2 sum3 sum4);
    for my $s (sort {$a <=> $b} keys %droidSupport) {
	my $ds = $droidSupport{$s};
	my @sums = map { $ds->{$_} } @cols;
	my @rates = map { sprintf("%.3e", $ds->{$_} / $ds->{n}) } @cols;
	say $OUT join "\t", $s, @rates, @sums;
    }
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	droid => '/home/glocke/DPiM/droid/support.updateFBgn.net',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

# ret[i] = {n1, n2, score};
sub edges {
    my ($ret, $net) = @_;

    for my $n1 (keys %$net) {
	for my $n2 (keys %{ $net->{$n1} }) {
	    push @$ret, {n1 => $n1, n2 => $n2, score => $net->{$n1}{$n2}};
	}
    }

}
