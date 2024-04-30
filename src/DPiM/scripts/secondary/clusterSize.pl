#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineDP4APMS);

# given a bunch of covers/cluster maps, find the number of nodes/edges in each

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    
    open my $IN, "<", $in or die "can't read $in. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    close $IN; 
    close $OUT; 
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	paramregex = 
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in clust.list -out output < -net findEdges.in.net ".
	">\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "net=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
