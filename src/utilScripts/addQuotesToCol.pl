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
    my $colN = $opts{col};
    my $out = $opts{out};

    open my $IN, "<", $in or die "can't open $in. $!";
    open my $OUT, ">", $out or die "can't open $out. $!";

    while (<$IN>) {
	if (/^#/) { 
	    print $OUT $_; 
	    next; 
	}

	chomp;
	my @spl = split /\t/; 
	$spl[$colN] = qq("$spl[$colN]"); 
	say $OUT join "\t", @spl;
    }

    close $IN;
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

    my $usage = "usage: $0 -in input -col n -out output\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "col=i", "out=s");
    die $usage unless exists $opts{in} && exists $opts{col} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
