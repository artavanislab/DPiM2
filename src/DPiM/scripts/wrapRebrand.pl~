#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineDP4APMS);

# describe the purpose for this script

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    
    open my $IN, "<", $in or die "can't read $in. $!";
    close $IN; 
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 describe this";
    close $OUT; 
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

sub getCommandLineOptions {
    # mode:
    my @modes = qw(mode1 mode2);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	);
    my $defaultString = join " ", map { "-$_ $defaults{$_}" } keys %defaults;
    my $usage = "usage: $0 -in input -out output < $modeString>\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{mode} //= $modes[0];
    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" if ! exists $modes{$opts{mode}};

    checkExist('f', $opts{in});

    return %opts;
}

