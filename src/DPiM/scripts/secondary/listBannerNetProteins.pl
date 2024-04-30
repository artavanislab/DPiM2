#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readCols readColsHash readHeader);
use DpimLib qw(networkHashFromEdgeList);

# append a column to flagnet output
# this column lists every protein in the subnetwork associated with this term
#  - along with degree

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $netFile = $opts{net};
    my $goListFile = $opts{golist};
    my $out = $opts{out};

    my %net;
    networkHashFromEdgeList(\%net, $netFile, 'symmetric');
    my %goList; # goList{GO:012345} = [protein1, protein2,...]
    #                                 list of all proteins with this term
    {
	for (readCols($goListFile, [qw(term proteins)], undef, "\t")) {
	    my @proteins = split /,/, $_->{proteins};
	    $goList{$_->{term}} = \@proteins;
	}
    }

    
    my @in = readCols($in, ['term'], 'line', "\t");
    my @cols = readHeader($in);
    push @cols, 'protein_degree';
    
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# added protein_degree column to $in";
    say $OUT "# degree counts only banner edges"; 
    say $OUT join "\t", @cols;
    for my $row (@in) {
	my $term = $row->{term};
	my $proteins = $goList{$term};
	my %bannerProteins; # these proteins have at least one banner edge
	for my $i (0..($#$proteins-1)) {
	    my $p1 = $proteins->[$i];
	    for my $j (($i+1)..$#$proteins) {
		my $p2 = $proteins->[$j];
		if (exists $net{$p1}{$p2}) {
		    $bannerProteins{$p1}++;
		    $bannerProteins{$p2}++;
		}
	    }
	}
	my $protString = proteinString($proteins, \%bannerProteins);
	chomp $row->{line};
	say $OUT $row->{line}, "\t", $protString;
    }
    close $OUT;
    
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my @modes = qw(edge node banner);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in base.flagnet -deg degree.dist -golist ".
	"golist.tsv -out output < -mode $modeString >\n";
    $usage .= "see also $ENV{DPSCR}/secondary/degreeDistribution.pl , $ENV{DPSCR}/secondary/compareGOHists.pl\n" if exists $ENV{DPSCR};
    
    my %opts = ();
    GetOptions(\%opts, "in=s", "deg=s", "golist=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{deg} && 
	exists $opts{golist} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{deg});
    checkExist('f', $opts{golist});

    return %opts;
}

# return fbgn01234_deg where deg is the degree of this protein (an integer)
sub proteinString {
    my ($proteins, $degree) = @_;

    my %degHash;
    for my $p (@$proteins) {
	my $deg = $degree->{$p};
	$degHash{$p} = $deg;
    }
    my @sorted = sort {$degHash{$b} <=> $degHash{$a}} @$proteins;
    my @elems = map { $_.'_'.$degHash{$_} } @sorted;
    return join ",", @elems;
}
