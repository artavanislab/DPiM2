#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist writeCols);
#use Graph::GML; # doesn't capture attributes !?  v0.01....
#use DpimLib qw(getLineAPMS);

# convert a gml file into an edge list
#
# not maximally functional.  just serves a single purpose, namely writing a
# undirected graph with node 'names' specified in a gml file

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    my (%nodes, %edges);
    readGML($in, \%nodes, \%edges);

    my @edgeList;
    my %cols;
    for my $n1 (sort keys %edges) {
	my $node1 = $nodes{$n1}{name} // die "can't find nodes{$n1}{name}";
	for my $n2 (sort keys %{ $edges{$n1} }) {
	    my $node2 = $nodes{$n2}{name} // die "can't find nodes{$n2}{name}";
	    my %row = (node1 => $node1, node2 => $node2);
	    for my $k (keys %{ $edges{$n1}{$n2} }) {
		next if $k eq 'source' || $k eq 'target';
		$row{$k} = $edges{$n1}{$n2}{$k};
		$cols{$k} = 1;
	    }
	    push @edgeList, \%row;
	}
    }
    #die Dunmp
    my @cols = ('node1', 'node2', sort keys %cols);
    my @data;
    for my $c (@cols) {
	push @data, [ map {$_->{$c}} @edgeList ];
    }
    my $header = join "\t", @cols;
    my $format = join "\t", ('%s') x @cols;
    my $preComments = "# made edgelist from $in";
    writeCols($out, \@data, $header, $format, $preComments);
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

    my $usage = "usage: $0 -in in.gml -out out.edgelist\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

# assume undirected
# nodes{id} = {id => id, attribute1 => a1, attr2 => a2,....}
# edges{id1}{id2} = { source, target, attributes... }
sub readGML {
    my ($in, $nodes, $edges) = @_;

    open my $IN, "<", $in or die "Can't read from $in. $!";

    my %elemTypes = ( node => 'node', edge => 'edge' );
    my $started = undef;
    my ($elemType, $elemOpen, %elem);
    while (my $line = <$IN>) {
	if (! $started) {
	    $started = 1 if $line =~ /^\[/;
	}
	$_ = $line;
	chomp;
	my @spl = split;
	if (! $elemOpen) {
	    if (exists $elemTypes{$spl[0]}) { # open an element
		$elemType = $spl[0];
		$elemOpen = 1;
	    } 
	} else {
	    if (@spl > 1) {
		# parse attributes
		my $key = shift @spl;
		my $val = join " ", @spl;
		$val =~ s/['"\]]//g;
		$elem{$key} = $val;
	    }
	    if ($line =~ /\]/) {
		# close element
		if ($elemType eq 'node') {
		    my $id = $elem{id} 
		        // die "can't find node id. \n", Dumper(\%elem);
		    $nodes->{$id} = { map { $_ => $elem{$_} } keys %elem };
		} elsif ($elemType eq 'edge') {
		    die "can't find source.\n", Dumper(\%elem) 
			unless exists $elem{source};
		    die "can't find target.\n", Dumper(\%elem) 
			unless exists $elem{target};
		    ($elem{source}, $elem{target}) = 
			sort ($elem{source}, $elem{target});
		    
		    my $source = $elem{source};
		    my $target = $elem{target};
		    
		    $edges->{$source}{$target} = {
			map { $_ => $elem{$_} } keys %elem 
		    };
		} 
		    
		%elem = ();
		$elemType = undef;
		$elemOpen = undef;
	    }
	}
    }
}
