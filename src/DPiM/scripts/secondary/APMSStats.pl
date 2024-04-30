#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineAPMS);

# find the total # of interactions and the number of runs in an APMS file

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};

    my %sID;
    my $count;

    open my $IN, "<", $in or die "can't read from $in. $!";
    my %row;
    while (getLineAPMS(\%row, $IN)) {
	$count++;
	$sID{$row{search_id}} = 1;
    }
    my $runs = 0+ keys %sID;
    
    say "$count interactions";
    say "$runs batches";
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

    my $usage = "usage: $0 -in apms.file\n";

    my %opts = ();
    GetOptions(\%opts, "in=s");
    die $usage unless exists $opts{in};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
