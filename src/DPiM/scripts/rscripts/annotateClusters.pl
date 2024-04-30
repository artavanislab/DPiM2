#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);
#use DpimLib qw(getLineDP4APMS);

# run annotateClusters.R on a list of files

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $ext = $opts{ext};
    my $scr = $opts{scr};
    my $mode = $opts{mode};

    my @files = readList($in);
    checkExist('f', $_) for @files;

    ## usage: clustPerLine.txt annotation.out < *fly*,human exclusions >
    for my $f (@files) {
	my $cmd = "$scr $f $f.$ext $mode";
	say $cmd;
	system($cmd);
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   


sub getCommandLineOptions {
    # mode:
    my @modes = qw(fly human);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	scr => $ENV{DPSCR}."/rscripts/annotateClusters.R",
	);
    my $defaultString = join " ", map { "-$_ $defaults{$_}" } keys %defaults;
    my $usage = "usage: $0 -in input -ext output < $modeString ".
	"$defaultString>\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "ext=s", "mode=s", "scr=s");
    die $usage unless exists $opts{in} && exists $opts{ext};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{mode} //= $modes[0];
    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" if ! exists $modes{$opts{mode}};

    checkExist('f', $opts{in});
    checkExist('f', $opts{scr});
    
    return %opts;
}

