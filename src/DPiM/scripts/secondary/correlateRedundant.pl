#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use Statistics::Basic qw(correlation);
use HomeBrew::IO qw(checkExist readList readCols writeCols);
use DpimLib qw(readTable);

# make a table showing the correlation between redundant runs
# correlate redundant runs against the "best" replicate, defined as having the
# highest bait TSC 

my %opts = getCommandLineOptions();

{
    my $in  = $opts{in};
    my $out = $opts{out};
    my $degreeFile = $opts{degree};
    
    my @cols = qw( cor bait search_id replicateRank degree baitCount 
                   bestBaitCount bestSearchID );

    my %degrees = getDegree($degreeFile);
    
    my @output; # output[] = { $cols[] => '', }
    for my $file (readList($in)) {
	my @colNames;
	my @table = readTable($file, \@colNames, 'prey');
	my $baitRow = shift @table;
	my $bait = $baitRow->{prey};
	my $degree = $degrees{$bait} // 0;
	
	my @search_id = @colNames;

	my $bestSearchID = shift @search_id;
	my $bestBaitCount = $baitRow->{$bestSearchID};

	my %baitCount = map {$_ => $baitRow->{$_}} @search_id;
	
	my @bestTSC = map { $_->{$bestSearchID} } @table;

	my %baseRow = ( cor => 0, bait => $bait, search_id => 0, 
			replicateRank => 0, degree => $degree, 
			baitCount => 0,
			bestBaitCount => $bestBaitCount, 
			bestSearchID => $bestSearchID,
	    );
	my $replicateRank = 2;
	for my $sID (@search_id) {
	    my @TSC = map { $_->{$sID} } @table;
	    my $cor = 0+ correlation(\@bestTSC, \@TSC);

	    my %row = map { $_ => $baseRow{$_} } @cols;
	    $row{cor} = $cor;
	    $row{search_id} = $sID;
	    $row{replicateRank} = $replicateRank;
	    $row{baitCount} = $baitCount{$sID};
	    
	    push @output, \%row;

	    $replicateRank++;				     
	}
    }

    my @data;
    for my $c (@cols) {
	my @co = map { $_->{$c} // die "can't find $c ??" } @output;
	push @data, \@co;
    }

    #my @cols = qw( cor bait search_id replicateRank degree baitCount 
     #              bestBaitCount bestSearchID );
    my $header = join "\t", @cols;
    my $format = join "\t", qw(%.4f %s %d %d %d %d %d %d);
    my $preComments = "# found correlations among files listed in $in";
    writeCols($out, \@data, $header, $format, $preComments);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	degree => '/home/glocke/DPiM/data/dpim3.degreeDistribution.tab',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in tab.list -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "degree=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{degree});

    return %opts;
}

sub getDegree {
    my ($in) = @_;

    my @cols = qw(protein degree);
    my @read = readCols($in, \@cols);

    return map { $_->{protein} => $_->{degree} } @read;
}
