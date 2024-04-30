#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readCol);
#use DpimLib qw(getLineDP4APMS);

# given a list of nodes, find all edges containing at least one of these nodes

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $fbgnString = $opts{fbgn};
    my $nodeFile = $opts{nodefile};
    my $n = $opts{n};

    my %fbgn; 
    if (defined $fbgnString) {
	%fbgn = map {$_ => 1} split /,/, $fbgnString;
    } elsif (defined $nodeFile) {
	my @nodes = readCol($nodeFile, 'protein');
	if (defined $n) {
	    @nodes = @nodes[0..($n-1)];
	}
	%fbgn = map {$_ => 1} @nodes;
	$fbgnString = join ",", @nodes;
    }

    my $testEdge = \&inclusiveTest;
    $testEdge = \&strictTest if exists $opts{strict};
    
    open my $IN, "<", $in or die "can't read $in. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 pulling all edges including nodes $fbgnString from $in"
	unless exists $opts{abcformat};
    my $line;
    if (! exists $opts{abcformat}) {
	do {
	    $line = <$IN>;
	    print $OUT $line;
	} until (eof($IN) || $line =~ /^protein1/);
    }

    while (<$IN>) {
	#next if $opts{abcformat} && $_ !~ /^FBgn/;
	next if $opts{abcformat} && $_=~/^Bait/; #human file header must start with Bait
	my @spl = split;
	
	print $OUT $_ if $testEdge->($spl[0], $spl[1], \%fbgn);
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

    my $usage = "usage: $0 -in in.net -out output [-fbgn FBgn1,FBgn2,...] or [-nodefile -n] < -abcformat -strict > \n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "fbgn=s", "nodefile=s", "n=i", 
	       "abcformat", "strict");
    die $usage unless exists $opts{in} && exists $opts{out} && 
	(exists $opts{fbgn} || exists $opts{nodefile});

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{nodefile}) if exists $opts{nodefile};

    return %opts;
}

# return true if this edge includes a sought node
sub inclusiveTest {
    my ($node1, $node2, $seek) = @_;
    return exists $seek->{$node1} || exists $seek->{$node2};
}

# return true if this edge connects two sought nodes
sub strictTest {
    my ($node1, $node2, $seek) = @_;
    return exists $seek->{$node1} && exists $seek->{$node2};
}
