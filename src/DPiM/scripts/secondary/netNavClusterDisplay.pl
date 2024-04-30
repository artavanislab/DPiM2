#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempfile);
use Data::Dumper;
#use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineAPMS);

# make a csv file that netnav can read in that displays mcl clustering results
# show only edges within clusters
# output format:
# protein1 protein2 cutoff mcl_parameter

#my %opts = getCommandLineOptions();

{
    my $baseDir = '/home/glocke/DPiM/dpim4/witInstr';
    my $outFile = "$baseDir/netnav.clusterDisplay.csv";
    my $clusterDir = 'mcl2/mcl.clusters';
    my $pull = '/home/glocke/DPiM/scripts/secondary/pullComplexMembers.pl';
    my $translate = '/home/glocke/DPiM/scripts/secondary/fbgnToGeneSymbol.pl';
    my $minSize = 2; # do not write edges for clusters smaller than this
    # at 2, this will report all clusters, as "clusters" with 1 member HAVE NO EDGES
	
    my @subDirs = qw(
percentileMin
percentile30
percentile34
percentile40
percentile50
percentile60
);
    
    my %cutoff = qw(
percentile30 992
percentile34 950
percentile40 790
percentile50 497
percentile60 288
percentileMin 999);

    my %netFiles = qw(
    percentile30 cons1_01-30-2016_30Percent.tsv
percentile34 cons1_01-30-2016_34Percent.tsv
percentile40 cons1_01-30-2016_40Percent.tsv
percentile50 cons1_01-30-2016_50Percent.tsv
percentile60 cons1_01-30-2016_60Percent.tsv
percentileMin cons1_01-30-2016_MinPercent.tsv);

    my ($FH, $tmpFile) = tempfile();
    close $FH;

    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";

    my $clusterID = 0;
    for my $dir (@subDirs) {
	my $netFile = "$baseDir/$dir/$netFiles{$dir}";
	my $cutoff = $cutoff{$dir};
	my $clDir = "$baseDir/$dir/$clusterDir";
	my @mcl = `ls $clDir/*txt`;
	chomp @mcl;
	for my $mclFile (@mcl) {
	    my $mclParam = findMCLParam($mclFile);
	    open my $MCLIN, "<", $mclFile or die "Can't read from $mclFile. $!";
	    while (<$MCLIN>) {
		chomp;
		my @clusterMembers = split;
		next if @clusterMembers < $minSize;
		my $cmd = "$pull -net $netFile -out $tmpFile -nodes ".
		    (join ",", @clusterMembers);
		say $cmd;
		system($cmd);
		my @nodePairs = readOneClust($tmpFile);
		for my $pair (@nodePairs) {
		    # protein1 protein2 cutoff mcl_parameter
		    say $OUT join ",", $pair->[0], $pair->[1], $cutoff, 
		        $mclParam, $clusterID;
		}
		$clusterID++;
	    }
	    close $MCLIN;
	}
    }
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

    my $usage = "usage: $0 -in input -out output\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    #checkExist('f', $opts{in});

    return %opts;
}

sub findMCLParam {
    my ($mclFile) = @_;
    
    $mclFile =~ m/format.i([\d\.]+).txt/ or die "can't parse $mclFile";

    return $1;
}


sub readOneClust {
    my $inFile = shift;

    open my $IN, "<", $inFile	or die "can't read from $inFile. $!";

    my @ret;
    while (<$IN>) {
	my @spl = split;
	push @ret, [$spl[0], $spl[1]];
    }
    return @ret;
}
