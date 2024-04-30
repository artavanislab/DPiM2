#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist writeCols);
use DpimLib qw(networkHashFromEdgeList);

# * take a network with string ids and convert it to integer ids
# * report the result without a header, and provide log for converting int
#   back to string
#
# symmetrize the network (report a b x and b a x) unless -nosym is set

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $noSym = $opts{nosym};
    
    open my $OUT, ">", $out or die "can't write to $out. $!";

    my %network;
    networkHashFromEdgeList(\%network, $in, undef, 'keep score');

    my %idmap;
    my $nextID = 1;
    for my $string1 (keys %network) {
	my $int1;
	if (exists $idmap{$string1}) {
	    $int1 = $idmap{$string1};
	} else {
	    $int1 = $idmap{$string1} = $nextID;
	    $nextID++;
	}
	for my $string2 (keys %{$network{$string1}}) {
	    my $int2;
	    if (exists $idmap{$string2}) {
		$int2 = $idmap{$string2};
	    } else {
		$int2 = $idmap{$string2} = $nextID;
		$nextID++;
	    }
	    if ($opts{noweight}) {
		say $OUT "$int1 $int2";
		say $OUT "$int2 $int1" unless $noSym;
	    } else {
		say $OUT "$int1 $int2 $network{$string1}{$string2}";
		say $OUT "$int2 $int1 $network{$string1}{$string2}" 
		    unless $noSym;
	    }
	}
    }
    close $OUT;

    my $log = "$out.log";
    my @k = sort {$idmap{$a} <=> $idmap{$b}} keys %idmap;
    my @num = map {$idmap{$_}} @k;
    my $header = "nodeString nodeInt";
    my $format = "%s\t%d\n";
    my $preComments = "# converted node names from $in";
    writeCols($log, [\@k, \@num], $header, $format, $preComments);
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

    my $usage = "usage: $0 -in net.in -out output < -nosym -noweight > \n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "nosym", "noweight");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
