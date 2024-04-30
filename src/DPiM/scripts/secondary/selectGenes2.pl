#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist readHeader);

# select edges that have certain genes in them

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $protString = $opts{prot};
    my $out = $opts{out};
    my $colString = $opts{cols};

    my %lookFor = map { $_ => 1 } split /,/, $protString;
    my @cols = split /,/, $colString;

    my $reader;
    {
	open my $IN, "<", $in or die "can't read from $in. $!";
	$reader = sub {
	    if ($_ = <$IN>) {
		chomp;
		return [split];
	    } else {
		close $IN;
		return undef;
	    }
	}
    }

    open my $OUT, ">", $out or die "Can't write to $out. $!";

    say $OUT "# seeking proteins $protString in columns $colString";
    say $OUT "# selecting lines listed in $in";

    while (my $row = $reader->()) {
	my ($p1, $p2) = ($row->[$cols[0]], $row->[$cols[1]]);
	next unless exists $lookFor{$p1} && exists $lookFor{$p2};
	say $OUT join "\t", @$row;
    }
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	cols => '0,1',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in selectFrom network -prot fbgn1,fbgn2,... ".
	"-out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "prot=s", "out=s", "cols=s");
    die $usage unless exists $opts{in} && exists $opts{prot} &&
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

