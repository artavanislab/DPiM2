#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(networkHashFromEdgeList);

# given edges in one order, repeate the order using edges from another network

my %opts = getCommandLineOptions();

{
    my $netFile = $opts{net};
    my $edgeFile = $opts{edge};
    my $out = $opts{out};
    my $newCol = $opts{newcol};
    my $null = $opts{notfound};

    $null = '' if $null eq '""' || $null eq "''";
    
    my %net;
    networkHashFromEdgeList(\%net, $netFile, 'symmetric', 'keepScore');
    
    my @edges;
    open my $IN, "<", $edgeFile or die "can't read $edgeFile. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    while (my $line = <$IN>) {
	chomp $line;
	if ($line =~ /^#/ || $line !~ /FBgn\d+/) {
	    $line.= "\t$newCol" if $line !~ /^#/;
	    print $OUT $line;
	}
	$_ = $line;
	my ($n1, $n2) = split;
	say $OUT join "\t", $line, $net{$n1}{$n2} // $null;
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
	newcol=>'support',
	notfound => '0',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -net net.in -edge edge.in -out output < ".
	"$defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "net=s", "edge=s", "out=s", "newcol=s", "notfound=s");
    die $usage unless exists $opts{net} && exists $opts{edge} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});
    checkExist('f', $opts{edge});

    return %opts;
}
