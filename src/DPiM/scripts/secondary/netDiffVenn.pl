#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(networkHashFromEdgeList);

# find the number of shared edges among (for now) precisely three networks

my %opts = getCommandLineOptions();

{
    my $netString = $opts{nets};
    my $out = $opts{out};

    my @files = split /,/, $netString;
    checkExist('f', $_) for @files;

    die "hackish implementation" if @files != 3;

    my @nets;
    for my $f (@files) {
	my $net = {};
	fillNet($net, $f);
	push @nets, $net;
    }

    my @alone = map { 0+ keys %$_ } @nets;
    
    my $center = sizeOfUnion3(@nets);

    say "center = $center";
    $_ -= $center for @alone;
    
    my @pairWise;
    for my $i (0..($#nets-1)) {
	for my $j (($i+1)..$#nets) {
	    my $pairwise = sizeOfUnion($nets[$i], $nets[$j]) - $center;
	    print $pairwise;
	    print "\t";
	    $alone[$i] -= $pairwise;
	    $alone[$j] -= $pairwise;
	}
	print "\n";
    }

    say $_ for @alone;
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

    my $usage = "usage: $0 -nets net1,net2,... -out tab\n";

    my %opts = ();
    GetOptions(\%opts, "nets=s", "out=s");
    die $usage unless exists $opts{nets} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }


    return %opts;
}

sub fillNet {
    my ($ret, $in) = @_;

    my %net;
    networkHashFromEdgeList(\%net, $in, 0, 1, 1);
    for my $k1 (keys %net) {
	for my $k2 (keys %{$net{$k1}}) {
	    my $edge = "$k1$k2";
	    $ret->{$edge} = 1;
	}
    }

    return;
}

sub sizeOfUnion {
    my ($hash1, $hash2) = @_;

    my $cnt = 0;
    for my $k (keys %$hash1) {
	$cnt++ if exists $hash2->{$k};
    }

    return $cnt;
}

sub sizeOfUnion3 {
    my ($hash1, $hash2, $hash3) = @_;

    my $cnt = 0;
    for my $k (keys %$hash1) {
	$cnt++ if exists $hash2->{$k} && exists $hash3->{$k};
    }

    return $cnt;
}
