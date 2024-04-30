#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);

# describe the purpose for this script

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $dir = $opts{dir};

    my @l = readList($in);
    checkExist('f', @l) for @l;

    for my $f (@l) {
	system("cp $f $dir/.");
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	);

    my $usage = "usage: $0 -in list.in -dir copyToHere\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "dir=s");
    die $usage unless exists $opts{in} && exists $opts{dir};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

