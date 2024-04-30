#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);
use feature ':5.10'; 

# given a list of scripts, locate within each script the "out" argument.
#   if the output file doesn't exist, write the name of the script to a new list

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    my @scripts = readList($in);
    my @outFiles = map { findOut($_) } @scripts;

    open my $OUT, ">", $out or die "can't write to $out. $!";
    for my $i (0..$#scripts) {
	next if -e $outFiles[$i];
	next if -s $outFiles[$i];
	say $OUT $scripts[$i];
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

    my $usage = "usage: $0 -in bash.list -out failed.list\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

# find the -out argument and report what it is
sub findOut {
    my $script = shift;

    open my $IN, "<", $script or die "can't read from $script. $!";

    my $ret;
    while (my $line = <$IN>) { 
	next unless $line =~ /-out\s+(\S+)\s/; # note that \s catches "\n"
	$ret = $1;
    }

    die "can't find output argument in $script.  quitting." unless $ret;

    return $ret;
}
