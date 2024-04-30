#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineDP4APMS);

# find experiments without any bait

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    my $freq = $ENV{DPSCR}."/calcCommonContam_KJ_GL.pl";
    my $stats = $ENV{DPSCR}."/baitPullDown.pl";
    my $select = $ENV{DPSCR}."/selectNoBait.pl";
    my $rebrand = $ENV{DPSCR}."/rebrandBait_rare.pl";
    checkExist('f', $_) for ($stats, $select, $rebrand);
    
    my $freqFile = "$out.working.input.freq";
    my $statsFile = "$out.working.input.stats";
    my $withBaitFile = "$out.working.withBait";
    my $noBaitFile = "$out.working.noBait";

    {
	my $cmd = "$freq -in $in -out $freqFile";
	say $cmd;
	system($cmd);
    }

    {
	my $cmd = "$stats -in $in -out $statsFile";
	say $cmd;
	system($cmd);
    }

    {
	my $cmd = "$select -apms $in -stats $statsFile -nobait $noBaitFile -clean $withBaitFile";
	say $cmd;
	system($cmd);
    }

    {
	my $cmd = "$rebrand -nobait $noBaitFile -withbait $withBaitFile -freq $freqFile -out $out";
	say $cmd;
	system($cmd);
    }
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

    my $usage = "usage: $0 -in input -out output\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
