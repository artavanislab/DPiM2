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

use lib '/home/glocke/DPiM/scripts/lib';
use DpimLib qw(getLineAPMS);

# find the distribution of # of bad batches per day

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $apmsFile = $opts{apms};
    my $pbinMax = $opts{pbinmax};
    
    my %date = mapBatchToDate($apmsFile);
    my %qual = mapBatchToQual($in, $pbinMax);

    my %perDay = dayDistribution(\%date, \%qual);
    hyperGeom(\%perDay, \%qual);
    
    writer($out, \%perDay);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	pbinmax => 0.01,
	apms => '/home/glocke/DPiM/data/dpim3.12232014.nrBait',	
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "pbinmax=f", "apms=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{apms});

    return %opts;
}

# ret{$search_id} = sample_date (string)
sub mapBatchToDate {
    my ($apmsFile) = @_;

    my %ret;
    
    open my $IN, "<", $apmsFile or die "Can't read from $apmsFile. $!";
    my %row;
    while (getLineAPMS(\%row, $IN)) {
	my $id = $row{search_id};
	my $date = $row{sample_date};
	if (exists $ret{$row{search_id}}) {
	    die "different dates. id = $id.  old = '$ret{$id}', new = '$date'."
		if $ret{$id} ne $date;
	} else {
	    $ret{$id} = $date;
	}
    }
    close $IN;

    #die Dumper(\%ret);
    
    return %ret;
}

# ret{search_id} = 0 for a bad run (less than max), 1 for good run (more than max)
sub mapBatchToQual {
    my ($in, $pbinMax) = @_;

    my @cols = qw(search_id pbin);
    my @read;
    readColsRef(\@read, $in, \@cols);

    my %ret = map { $_->{search_id} => (($_->{pbin} > $pbinMax)?1:0) } @read;

    return %ret;
}

# ret{date} = {good = #, bad = #}
sub dayDistribution {
    my ($date, $qual) = @_;

    my %qualMap = (0 => 'bad', 1 => 'good');
    my %ret;
    for my $id (keys %$date) {
	my $date = $date->{$id};
	$ret{$date} //= {good => 0, bad => 0};
	my $q = $qualMap{$qual->{$id}};
	$ret{$date}{$q}++;
    }

    return %ret;
}

sub hyperGeom {
    my ($perDay, $qual) = @_;

    my ($goodTot, $badTot) = qualStats($qual);

    for my $d (keys %$perDay) {
	my $found = $perDay->{$d}{good};
	my $outOf = $perDay->{$d}{bad} + $perDay->{$d}{good};
	$perDay->{$d}{pbin} = 
	    gsl_cdf_hypergeometric_P($found, $goodTot, $badTot, $outOf);
    }
}

# return (number that meet criteria, number that don't meet criteria)
sub qualStats {
    my ($qual) = @_;

    my $goodTotal = sum values %$qual;
    my $all = 0+ keys %$qual;
    
    return ($goodTotal, $all - $goodTotal);
}

sub writer {
    my ($out, $perDay) = @_;

    my @cols = qw(date good bad pbin);

    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT join "\t", @cols;
    for my $date (sort keys %$perDay) {
	say $OUT (join "\t", $date, $perDay->{$date}{good},
		  $perDay->{$date}{bad}, $perDay->{$date}{pbin});
    }
    close $OUT;
    return;
}
