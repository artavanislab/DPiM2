#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use lib '/home/glocke/DPiM/scripts/lib';
use DpimLib qw(getLineAPMS);

use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist);

# list all the proteins that have exactly one TSC across all AP-MS data

my %opts = getCommandLineOptions();

{
    my $apmsFile = $opts{apms};
    my $out = $opts{out};

    my %TSC;
    getTSC(\%TSC, $apmsFile);

    open my $OUT, ">", $out or die "Can't write t. $out . $!";
    say $OUT $_ for sort grep { $TSC{$_} == 1 } keys %TSC;
    close $OUT;

    if (exists $opts{all}) {
	$out = $opts{all};
	open my $OUT, ">", $out or die "Can't write t. $out . $!";
	say $OUT $_ for sort keys %TSC;
	close $OUT;
	
    }
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

    my $usage = "usage: $0 -out output < $defaultString -all allProt.list >\n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "apms=s", "all=s");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{apms});

    return %opts;
}

sub getTSC {
    my ($ret, $apmsFile) = @_;

    open my $IN, "<", $apmsFile or die "Can't read from $apmsFile. $!";

    my %row;
    while (getLineAPMS(\%row, $IN)) {
	for my $p ($row{bait_ref}, $row{prey_ref}) {
	    $ret->{$p} += $row{total_peptides};
	}
    }

    return;
}
