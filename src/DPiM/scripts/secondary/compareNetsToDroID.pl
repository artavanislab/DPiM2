#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(max);
use List::MoreUtils qw(each_array);
use HomeBrew::IO qw(checkExist readCols readHeader readList);
use DpimLib qw(networkHashFromEdgeList readHS);

# for the networks given, find the number of DroID edges captured at 
## varying levels of support

# TODO: find precision, recall, and f-statistic

my %opts = getCommandLineOptions();

{
    my $netList = $opts{in};
    my $out = $opts{out};
    my $droidFile = $opts{droid};
    my @names = split /,/, $opts{names};
    
    my @netFiles = readList($netList);
    checkExist('f', $_) for @netFiles;
    if ($opts{names} eq 'default') {
	@names = map {"net$_"} 1..@netFiles;
    }
    
    my %droid;
    networkHashFromEdgeList(\%droid, $droidFile, undef, 'keepScore');
    my %tot;
    my $maxS = 0;
    for my $prot1 (keys %droid) {
	for my $prot2 (keys %{$droid{$prot1}}) {
	    my $support = $droid{$prot1}{$prot2};
	    $tot{$_}++ for 1..$support;
	    $maxS = $support if $support > $maxS;
	}
    }


    my %droidSupport; # droidSupport{net} = { 1=>sum1, 2=> sum2, 3=>sum3...}
    my %size; # size{net} = number of edges in network
    
    my $it = each_array(@netFiles, @names);
    while (my ($netFile, $name) = $it->()) {
	my %net;
	networkHashFromEdgeList(\%net, $netFile, 'symm', 'keepScore');
	for my $p1 (keys %droid) {
	    next unless exists $net{$p1};
	    for my $p2 (keys %{ $droid{$p1} }) {
		next unless exists $net{$p1}{$p2};
		my $support = $droid{$p1}{$p2};
		$droidSupport{$name}{$_}++ for 1..$support;
	    }
	}
	$size{$name} = countEdges($netFile);
    }	

    my @supp = (1..$maxS);

    my (%recall, %precision, %F);
    # $X{net}{1..5} = statistic for this true set
    for my $net (@names) {
	my $edges = $size{$net};
	for my $s (@supp) {
	    my $r = $droidSupport{$net}{$s}/$tot{$s};
	    my $p = $droidSupport{$net}{$s}/$edges;
	    $recall{$net}{$s} = $r;
	    $precision{$net}{$s} = $p;
	    $F{$net}{$s} = 2 * $p*$r / ($p+$r);
	}
    }
      
    

    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 found droid support of ".(join ", ", @netFiles).
	" according to $droidFile";
    say $OUT join "\t", 'support', (
	map { ("recall-$_", "precision-$_", "F-$_") } @names), "tot";
    for my $s (@supp) {
	print $OUT $s;
	for my $net (@names) {
	    printf $OUT "\t%.3e\t%.3e\t%.3e", $recall{$net}{$s}, 
	        $precision{$net}{$s}, $F{$net}{$s};
	}
	print $OUT "\t$tot{$s}";
	print $OUT "\n";
    }
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	droid => '/home/glocke/DPiM/droid/support.updateFBgn.net',
	names => 'net1,net2,...',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "droid=s", "names=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{names} = 'default' if $opts{names} eq $defaults{names};
    
    checkExist('f', $opts{in});

    return %opts;
}

sub countEdges {
    my ($inFile) = @_;

    my $ret = 0;
    open my $IN, "<", $inFile or die "Can't read from $inFile: $!";
    while (my $line = readHS($IN)) {
	$ret++;
    }

    return $ret;
}
