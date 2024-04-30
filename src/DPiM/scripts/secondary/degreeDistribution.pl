#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist writeCols);

use lib '/home/glocke/DPiM/scripts/lib';
use DpimLib qw(degreeDistribution);

# find and print the degree distributinon

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $in2 = $opts{in2};

    my %dist = degreeDistribution($in);
    my @k = sort {$dist{$b} <=> $dist{$a}} keys %dist;
    my @d = map {$dist{$_}} @k;

    
    my $header = join "\t", qw(rank protein degree);
    my $format = join "\t", qw(%d %s %d);
    my $preComments = "# found degree distribution from $in\n";
    my @data = (\@k, \@d);

    if (defined $in2) {
	my %dist2 = degreeDistribution($in2);
	my @d2 = map {$dist2{$_}} @k;
	$header.="\tdegree2";
	$format.="\t%d";
	$preComments.= "# degree2 refers to $in2";
	push @data, \@d2;
    }
    writeCols($out, \@data, $header, $format, $preComments, 1)
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

    my $usage = "usage: $0 -in input -out output < -in2 another.net >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "in2=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{in2}) if exists $opts{in2};

    return %opts;
}
