#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readCols);
#use DpimLib qw(getLineDP4APMS);

# update old fbgn's to new ones

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $refMap = $opts{ref};

    my %update = readRefMap($refMap);
    say "\%update constructed";

    open my $IN, "<", $in or die "Can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";

    my $cnt = 0;
    while (my $line = <$IN>) {
	for my $old (keys %update) {
	    if (exists $opts{findkeys}) {
		my $flag;
		for my $k (keys %update) {
		    $flag = 1 if $line =~ /$k/;
		}
		print $OUT $line if $flag;
	    } else {
		while ($line =~ s/$old/$update{$old}/g) {
		    $cnt++;
		}
	    }
	}
	print $OUT $line;
    }
    close $IN;
    close $OUT;

    say "made $cnt substitutions";
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	ref=>'/home/glocke/DPiM/prevDPIM/FlyBase_IDs_FB2016_01.txt',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString -findkeys>\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "ref=s", "findkeys");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{ref});

    return %opts;
}

sub readRefMap {
    my ($inFile) = @_;

    my @read = readCols($inFile, [qw(submitted_id current_id converted_id)]);
    my %ret;
    for my $row (@read) {
	die "$row->{current_id} doesn't match $row->{converted_id}"
	    unless $row->{current_id} eq $row->{converted_id};
	warn "already have a map for $row->{submitted_id}" 
	    if exists $ret{$row->{submitted_id}};
	next if $row->{submitted_id} eq $row->{current_id};
	$ret{$row->{submitted_id}} = $row->{current_id};
    }

    return %ret;
}
