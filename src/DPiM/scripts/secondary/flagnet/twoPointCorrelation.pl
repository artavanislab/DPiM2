#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Storable;
use Math::GSL::CDF qw(gsl_cdf_hypergeometric_P);
#use File::stat;

use HomeBrew::IO qw(checkExist readColsRef writeCols);
use DpimLib qw(readGoDB networkHashFromEdgeList);

# How much more likely is my neighbor to have a flag if I a flag?

# FlagNet samples edges, which is two GO lookups at a time.
# The question here is whether these lookups are correlated.

# Given that the question we're asking is whether our sampling is a good
#   approximation to hypergeometric sampling, it should come as no surprise that
#   a hypergeomtric test will come in handy!

# What is the null model?
#   Ultimately, the null model is that the flagged nodes are not clustered.
#   In the context of edge-wise sampling, this means that the conditional 
#   probability that one end is flagged when we know the other is flagged is 
#   equal to the probability of any edge-end being flagged.
# In this case, we're pulling edge-ends, and we don't pull the same edge-end 
#   twice.  This is hyper-geometric sampling more or less equivalent to the 
#   FlagNet test itself.

#            flagged  nonflagged 
#          +---------+---------+
# adjacent | ##      |         |
#          +---------+---------+
# the rest |         |         |
#          +---------+---------+


# algorithm: 
# For every flagged node, count the flagged and unflagged neighbors
# Perform HG test.  report the likelihood ratio and the p-value.

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $refHistFile = $opts{refhist};
    my $goDBFile = $opts{godb};
    my $goNameFile = $opts{goname};

    my %goDB;
    say "parsing $goDBFile";
    readGoDB(\%goDB, $goDBFile);

    my %invGoDB;
    invertGoDB(\%invGoDB, \%goDB);

    my %network;
    networkHashFromEdgeList(\%network, $in);
    for my $k1 (keys %network) {
	for my $k2 (keys %{ $network{$k1} }) {
	    $network{$k2}{$k1} = $network{$k1}{$k2};
	}
    }
    
    

    # for every term, seek every flagged node
    # examine every adjacent edge-end, counting which are/not flagged

    my %neigh; # neigh{term}{flags} = number of edge-ends flagged with term
    #          # neigh{term}{pulls}   = number of edge-ends tested
    my @terms;
    for my $term (keys %invGoDB) {
	say $term;
	for my $flagged (keys %{ $invGoDB{$term} }) {
	    say "\t", $flagged;
	    $neigh{$term}{flags} = 0;
	    $neigh{$term}{pulls} = 0;
	    for my $neighbor (keys %{ $network{$flagged} }) {
		#say "\t$neighbor";
		$neigh{$term}{pulls}++;
		$neigh{$term}{flags}++ if exists $invGoDB{$term}{$neighbor};
	    }
	}
	exit;
	push @terms, $term if $neigh{$term}{pulls} > 0;
	# keep only terms when there are annotated nodes in the network
    }
    
    my %refHist = readHist($refHistFile);
    my $refPulls = countProteins($refHistFile);

    my $goNameMap = retrieve($goNameFile);

    my @results;
    for my $t (@terms) {
	my $found = $neigh{$t}{flags};
	my $pulls = $neigh{$t}{pulls};
	my $white = $refHist{$t} // die Dumper($neigh{$t});
	my $black = $refPulls - $white;
	my $expect = $pulls * $white / $refPulls;

	my $name = $goNameMap->{$t}{name} // "UNKNOWN";

	my %row = (term => $t, found => $found, expect => $expect, 
		   white => $white, black => $black, pulls => $pulls, 
		   name => '"'.$name.'"');
	
	if ($found < $expect) {
	    $row{test} = '-'; 
	} else {
	    $row{test} = '+';
	    ($white, $black) = ($black, $white);
	    $found = $pulls - $found;
	}
	say "$t\t$name:\tgsl_cdf_hypergeometric_P($found, $white, $black, $pulls)";
	my $p = gsl_cdf_hypergeometric_P($found, $white, $black, $pulls);
	$row{p} = $p;

	push @results, \%row;
    }
    

    @results = sort {$a->{p} <=> $b->{p}} @results;
    
    my @cols = qw( term p test found expect white black pulls name);
    my @data;
    for my $c (@cols) {
	push @data, [map {$_->{$c}} @results];
    }
    die join ", ", map { $data[$_][0] } 0..$#cols;
    
    my $header = join "\t", @cols;
    my $format = join "\t", qw( %s %.3e %s %d %d %.2f %d %d %s );
    my $preComments = "# test = '-' is a test for extreme ".
	"underrepresentation\n";
    $preComments .= "# test = '+' is a test for extreme overrepresentation\n";
    $preComments .= "# q = p * number of tests\n";
    
    writeCols($out, \@data, $header, $format, $preComments);

}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	godb => '/home/glocke/DPiM/flybase/gene_association.fb',
	goname => '/home/glocke/DPiM/flybase/goMap.storable',
	refhist => '/home/glocke/DPiM/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOHist',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in full.network -out output < $defaultString > \n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "godb=s", "goname=s", "refhist=s", );
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

sub readHist {
    my ($in) = @_;
    my @cols = qw(term count);
    my @read;
    readColsRef(\@read, $in, \@cols);
    return map { $_->{term} => $_->{count} } @read;
}
sub countProteins {
    my ($in) = @_;

    open my $IN, "<", $in or die "Can't read from $in. $!";
    my $line;
    do {
	$line = <$IN>;
    } until ($line =~ /proteins found = (\d+)/);
    die "can't count proteins" unless defined $1 && $1 > 0;

    return $1;
}

# invGoDB{term}{prot} is defined iff this protein is flagged by term
sub invertGoDB {
    my ($invGoDB, $goDB) = @_;

    for my $prot (keys %$goDB) {
	for my $term (@{ $goDB->{$prot} }) {
	    $invGoDB->{$term}{$prot} = 1;
	}
    }

    return;
}
