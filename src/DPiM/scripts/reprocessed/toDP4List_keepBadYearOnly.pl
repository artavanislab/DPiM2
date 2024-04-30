#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Temp qw(tempfile);
use HomeBrew::IO qw(checkExist readList);
#use DpimLib qw(getLineDP4APMS);

# convert each file to dp4 format and then concatenate them

my %opts = getCommandLineOptions();

{
    my $listFile = $opts{in};
    my $out = $opts{out};
    my $toDP4 = $opts{todp4};
    my $catAPMS = $opts{catapms};

    my @inFiles = readList($listFile);
    checkExist('f', $_) for @inFiles;

    my @outFiles;
    for my $f (@inFiles) {
	my ($FH, $out) = tempfile();
	close $FH;
	push @outFiles, $out;
	my $cmd = "$toDP4 -in $f -out $out";
	say $cmd;
	system($cmd);
    }
    checkExist('f', $_) for @outFiles; # check that the above ran correctly

    my ($TMP, $tmpFile) = tempfile();
    say $TMP $_ for @outFiles;
    close $TMP;

    my $cmd = "$catAPMS -in $tmpFile -out $out";
    say $cmd;
    system($cmd);

	# removing badyear experiment (2009-04-01 till 2010-06-01) Added by Kejie, 08082018
	#$cmd = "/home/kli3/DPiM/scripts/removeBadDate.pl -in $out -out $tmpFile; mv $out $out.beforeRmBadYear; mv $tmpFile $out";
	$cmd = "/home/kli3/DPiM/scripts/keepBadDate.pl -in $out -out $tmpFile; mv $out $out.beforekeepBadYear; mv $tmpFile $out";
	system($cmd);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	todp4 => "$ENV{DPSCR}/reprocessed/toDP4.pl",
	catapms => "$ENV{DPSCR}/reprocessed/catAPMS.pl",
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "todp4=s", "catapms=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

