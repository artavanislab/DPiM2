#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineAPMS);

# given a list of node names, select all edges connecting those nodes from a 
# parent network

my %opts = getCommandLineOptions();

{
    my $nodeList = $opts{node};
    my $networkFile = $opts{net};
    my $out = $opts{out};
    my $mode = $opts{mode};

    my @clusters; # cluster[i] = { node1=> 1, node2=> 1, ...}
    if ($mode eq 'edgefiles') {
	readEdgeFiles(\@clusters, $nodeList);
    } elsif ($mode eq 'nodesbyline') {
	readNodesByLine(\@clusters, $nodeList);
    } else {
	die "Unknown mode '$mode'";
    }
    
    
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

    # mode:
    my @modes = qw(edgefiles nodesbyline);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my $usage = "usage: $0 -node in.list -net in.net -out output  ".
	"< $modeString>\n";

    my %opts = ();
    GetOptions(\%opts, "node=s", "net=s", "out=s", "mode=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

sub readEdgeFiles {
    my ($clusters, $in) = @_;

    my @files = readList($in);

    for my $f (@files) {
	my %this;

	open my $IN, "<", $f or die "Can't read from $f. $!";
	while (<$IN>) {
	    while (/(FBgn\d+)/g) {
		$this{$1} = 1;
	    }
	}
	close $IN;

	push @$clusters, \%this if 0 < keys %this;
    }

    return;
}

sub readNodesByLine {
    my ($clusters, $in) = @_;
    
    open my $IN, "<", $in or die "can't read from $in. $!";

    while (<$IN>) {
	my %this;
	while (/(FBgn\d+)/g) {
	    $this{$1} = 1;
	}
	push @$clusters, \%this if 0 < keys %this;
    }

    return;
}
