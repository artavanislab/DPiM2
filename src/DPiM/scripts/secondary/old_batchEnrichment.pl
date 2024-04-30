#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use lib '/home/glocke/perl5/lib/perl5';
use Math::GSL::CDF qw(gsl_cdf_hypergeometric_P);
    
use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist writeCols);

use lib '/home/glocke/DPiM/scripts/lib';
use DpimLib qw(readAPMS readHS);

# count which batches are enriched for listed edges (e.g. discarded)
# for each batch: find full size of batch, find number of discarded edges
# find global average of discarded edges (per batch? per edge?)
# apply 

my %opts = getCommandLineOptions();

{
    my $enrichFile = $opts{enrich};
    my $apmsFile = $opts{apms};
    my $out = $opts{out};

    my @apms;
    say "about to read $apmsFile";
    readAPMS(\@apms, $apmsFile);

    my $total = 0+@apms;
    
    my %batch; # batch{tap_id} = { total, enrich, pbin }
    # pbin is the probability of getting {enrich} given {total} pulls under
    # a binomial expectation  

    # find batch sizes
    for my $edge (@apms) {
	my $tid = $edge->{tap_id};
	$batch{$tid}{total}++;
    }

    say "parsing enriched edges";
    my $enrichTotal = 0;
    my $preyPreyCnt = 0;
    open my $ENR, "<", $enrichFile or die "can't read from $enrichFile. $!";
    while (my $line = readHS($ENR)) {
	$_ = $line;
	chomp;
	my ($n1, $n2, @spl) = split;
	my @find = findEdgeBatches($n1, $n2, \@apms);
	if (! defined $find[0]) {
	    $preyPreyCnt++;
	    next;
	    #die "can't find batch for '$n1' -> '$n2'" if 
	}
	$enrichTotal++;
	$batch{$_}{enrich}++ for @find;
	$total = $total - (@find - 1);
    }
    say "prey-prey count is $preyPreyCnt";
    #my $p = $enrichTotal / $total;
    my $white = $enrichTotal;
    my $black = $total - $enrichTotal;

    # hypergeometric test
    say "statistical tests";
    for my $b (keys %batch) {
	my $found = $batch{$b}{enrich};
	my $outOf = $batch{$b}{total};
	$batch{$b}{pbin} = 
	    gsl_cdf_hypergeometric_P($found, $white, $black, $outOf);
    }

    my @colNames = qw(enrich total pbin);
    my @data;
    my @k = sort keys %batch;
    push @data, \@k;
    for my $c (@colNames) {
	push @data, [map {$batch{$_}{$c}} @k]; 
    }
    unshift @colNames, "tap_id";

    my $header = (join "\t", @colNames)."\n";
    my $format = (join "\t", qw(%s %d %d %.4e))."\n";
    my $preComments = "# binomial probabilities assessed for a total of $white enriched edges and $black non-enriched edges";
    writeCols($out, \@data, $header, $format, $preComments);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	apms => '/home/glocke/DPiM/data/dpim3.12232014.nrBait',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -enrich edge.list -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "enrich=s", "out=s", "apms=s");
    die $usage unless exists $opts{enrich} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{enrich});
    checkExist('f', $opts{apms});

    return %opts;
}

# given a particular PPI, find which batch(es) it comes from
sub findEdgeBatches {
    my ($n1, $n2, $batch) = @_;

    my @find;
    push @find, 
    grep { $_->{bait_ref} eq $n1 && $_->{prey_ref} eq $n2 ||
	       $_->{prey_ref} eq $n1 && $_->{bait_ref} eq $n2 } @$batch;

    if (0+ @find) {
	return map { $_->{tap_id} } @find;
    } else {
	return undef;
    }
}
