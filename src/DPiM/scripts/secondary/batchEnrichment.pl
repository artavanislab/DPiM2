#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use lib '/home/glocke/perl5/lib/perl5';
use Math::GSL::CDF qw(gsl_cdf_hypergeometric_P);
    
use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist readColsRef writeCols);

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

    say "make enrichMap";
    my %enrichMap;
    my $enrichTotal = makeEnrichMap(\%enrichMap, $enrichFile);
    
    my @apms;
    say "read $apmsFile";
    readAPMS(\@apms, $apmsFile);

    my $total = 0+@apms;

    say "parse batches";
    my %batch; # batch{search_id} = { total, enrich, pbin }
    # pbin is the probability of getting {enrich} given {total} pulls under
    # a binomial expectation  
    my $cnt++;
    for my $edge (@apms) {
	my $tid = $edge->{search_id};
	$batch{$tid} //= {date => $edge->{sample_date}, 
			  bait_ref => $edge->{bait_ref}};
	$batch{$tid}{total}++;
	my $n1 = $edge->{bait_ref};
	my $n2 = $edge->{prey_ref};
	$batch{$tid}{enrich}++ if exists $enrichMap{$n1}{$n2};
	say "$cnt..." if ($cnt % 10000) == 0;
    }

    
    my $white = $enrichTotal;
    my $black = $total - $enrichTotal;

    say "hypergeometric tests";
    for my $b (keys %batch) {
	$batch{$b}{enrich} //= 0;
	my $found = $batch{$b}{enrich};
	my $outOf = $batch{$b}{total};
	$batch{$b}{pbin} = 
	    gsl_cdf_hypergeometric_P($found, $white, $black, $outOf);
    }

    my @colNames = qw(enrich total pbin date bait_ref);
    my @data;
    my @k = sort {$a <=> $b} keys %batch;
    push @data, \@k;
    for my $c (@colNames) {
	push @data, [map {$batch{$_}{$c}} @k]; 
    }
    unshift @colNames, "search_id";

    if (exists $opts{enrichedonly}) {
	# only print out info on batches with at least one enriched edge
	my @ind = grep { $data[1][$_] > 0 } 0..$#{$data[1]};
	my @newData;
	for my $col (0..$#data) {
	    my @sel = map { $data[$col][$_] } @ind;
	    push @newData, \@sel;
	}
	@data = @newData;
    }
    
    my $header = (join "\t", @colNames)."\n";
    my $format = (join "\t", qw(%s %d %d %.4e %s %s))."\n";
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

    my $usage = "usage: $0 -enrich edge.list -out output < $defaultString ".
	"-enrichedonly >\n";

    my %opts = ();
    GetOptions(\%opts, "enrich=s", "out=s", "apms=s", "enrichedonly");
    die $usage unless exists $opts{enrich} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{enrich});
    checkExist('f', $opts{apms});

    return %opts;
}

sub makeEnrichMap {
    my ($ret, $inFile) = @_;
        
    my @enrich;
    readColsRef(\@enrich, $inFile, [qw(node1 node2)], 'line');

    for my $row (@enrich) {
	my $n1 = $row->{node1};
	my $n2 = $row->{node2};
	say Dumper("undef", $row) if ! (defined $n1 && defined $n2);
	$ret->{$n1}{$n2} = $ret->{$n2}{$n1} = 1
    }

    return(0+ @enrich);
}
