#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readHeader);
use DpimLib qw(getLineDP4APMS getLineRawAPMS getLineHyperspecAPMS);

# When a prey fbgn appears more than once as in the same experiment, sum the 
# peptides so that it appears just once


my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};

    my %expBySID; # expBySID{sID}{prey} = \%row;

    my @cols = readHeader($in);
    my $reader = \&getLineDP4APMS;
    $reader = \&getLineRawAPMS if $mode eq 'raw';
    $reader = \&getLineHyperspecAPMS if $mode eq 'fourcols';
    
    
    open my $IN,  "<", $in or die "Cannot open $in. $!";
    my %row;
    while($reader->(\%row, $IN, 'line')) {
	my $sID = $row{search_id};
	my $prey = $row{prey_ref};
	if (exists $expBySID{$sID}{$prey}) {
	    $expBySID{$sID}{$prey}{total_peptides}+= $row{total_peptides};
	} else {
	    $expBySID{$sID}{$prey} = { map { $_ => $row{$_} } @cols };
	}
    }


    open my $OUT, ">", $out or die "Can't write to $out. $!";

    say $OUT join "\t", @cols;
    for my $sID (sort {$a <=> $b} keys %expBySID) {
	for my $prey (sort {$expBySID{$sID}{$b}{total_peptides} <=> 
				$expBySID{$sID}{$a}{total_peptides} } 
		      keys %{$expBySID{$sID}}) 
	{
	    say $OUT join "\t", map { $expBySID{$sID}{$prey}{$_} } @cols;
	}
    }

    exit;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my @modes = qw(dp4 raw fourcols);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $modeString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    checkExist('f', $opts{in});

    return %opts;
}


