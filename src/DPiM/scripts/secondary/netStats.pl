#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(max);
use HomeBrew::IO qw(checkExist readHeader readColsRef);
use Graph;
use Graph::Undirected;


#use DpimLib qw(networkHashFromEdgeList);

# find the following informations
# * number of edges
# * number of nodes
# * size of large connected component as a fraction of total nodes

my %opts = getCommandLineOptions();

{
    my $netFile = $opts{net};
    my $out = $opts{out};

    my $graph = Graph::Undirected->new; 
    {
	my @cols = readHeader($netFile);
	my @read;
	readColsRef(\@read, $netFile, \@cols);
	for my $row (@read) {
	    $graph->add_edge($row->{$cols[0]}, $row->{$cols[1]});
	}
    }

    my $size = 0+ $graph->vertices();
    my @cc = $graph->connected_components();
    my $LCCSize = max( map { 0+ @$_ } @cc);
    
    say "$size nodes";
    say "".(0+ $graph->edges())." edges";
    printf "%.4f LCCFraction\n", $LCCSize / $size;
    say "$LCCSize nodes in LCC";
    
    if (defined $out) {
	open my $OUT, ">", $out or die "can't write to $out. $!";
	close $OUT;
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

    my $usage = "usage: $0 -net in.net < -out out.txt > \n";

    my %opts = ();
    GetOptions(\%opts, "net=s", "out=s");
    die $usage unless exists $opts{net};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});

    return %opts;
}
