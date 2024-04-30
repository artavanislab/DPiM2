#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Storable qw(retrieve);
use HomeBrew::IO qw(checkExist readColsHashRef);
use DpimLib qw(readHS);

## in 5 source networks, find the cluster that is most enriched for this term
## generate table with 6+ columns:
##   node1 node2 symbol1 symbol2 score source 
## additional columns are the annotated term in different MCL parameters
## source networks: DPIM1 allRep bestRep CompPASS humanHG_allRep

my %opts = getCommandLineOptions();

{
    my $netFile = $opts{net};
    my $mclDir = $opts{mcldir};
    my $out = $opts{out};
    my $source = $opts{source};
    my $goNameFile = $opts{goname};
    my $onlyI = $opts{onlyi};
    
    my @i = qw(1.2 1.4 1.6 1.8 2 3 4 6 10 15 20);
    @i = ($onlyI) if defined $onlyI;
    
    my @cols = qw(protein1 protein2 symbol1 symbol2 score source mclparam name);
    
    my %symbolMap;
    readColsHashRef(\%symbolMap, $opts{humansymbol}, [qw(entrez symbol)]);
    makeMap(\%symbolMap, $opts{flysymbol});
    
    say "retrieving goNameMap";
    my $goNameMap = retrieve($goNameFile);
    # $goNameMap->{UNKNOWN} = { name => '"unknown"' };
    
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT join "\t", @cols;
    
    my %clusters; # cluster{i}{gene}{id}= 1 iff gene is in cluster #id
    my %clusterGO; # clusterGO{$i}
    for my $i (@i) {
	$clusters{$i} = {};
	$clusterGO{$i} = {};
	my @f = glob("$mclDir/*i$i.txt");
	die Dumper($netFile, $i, \@f) if 1 != @f;
	my $in = $f[0];
	clusterFromEachLine($clusters{$i}, $in);
	
	my $goTestFile = "$in.max5.GOTest";
	my @cols = qw(i term);
	@cols =  qw(i bestTerm) if testBestTerm($goTestFile);
	readColsHashRef($clusterGO{$i}, $goTestFile, \@cols);
    }

    open my $IN, "<", $netFile or die "Can't read from $netFile: $!";
    while (my $line = readHS($IN, $opts{human})) {
	chomp $line;
	my ($p1, $p2, $score) = split /\s+/, $line;
	my $s1 = $symbolMap{$p1} // die "can't find symbol for $p1";
	my $s2 = $symbolMap{$p2} // die "can't find symbol for $p2";
	#my @clusterNames = 
	#map { nameCluster($p1, $p2, $clusters{$_}, $clusterGO{$_},
	#$goNameMap) } @i;
	for my $i (@i) {
	    my $clusterName = nameCluster($p1, $p2, $clusters{$i}, 
					  $clusterGO{$i}, $goNameMap);
	    
	    say $OUT join "\t", $p1, $p2, $s1, $s2, $score, $source, $i
		, $clusterName;
	}
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
	goname => '/home/glocke/DPiM/flybase/goMap.full.storable',
	flysymbol => '/home/glocke/DPiM/flybase/fbgn_annotation_ID_fb_2015_04.tsv',
	humansymbol => '/home/glocke/DPiM/human/nsaf/entrez2symbol.concise.tsv',
	source => 'nrBait',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -net nrBait.net -mcldir pathTo/[mcl.txt] -out output < $defaultString -onlyi 3(optional) -human >\n";

    my %opts = ();
    GetOptions(\%opts, "net=s", "mcldir=s", "out=s", "nrBait=s", "goname=s", 
	       "flysymbol=s", "humansymbol=s", "source=s", "onlyi=i", "human");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});
    checkExist('d', $opts{mcldir});
    return %opts;
}

# make map from fbgn to gene symbol
sub makeMap {
    my ($ret, $mapFile) = @_;

    open my $IN, "<", $mapFile or die "Can't read from $mapFile. $!";
    while (<$IN>) {
	next if /^#/;
	next if length($_) < 2;
	chomp;
	my @spl = split;
	$ret->{$spl[1]} = $spl[0];
	if (@spl > 2) {
	    my @spl2 = split ',', $spl[2];
	    $ret->{$_} = $spl[0] for @spl2;
	}
	#die Dumper($ret);
    }
    return;
}

sub testBestTerm {
    my ($goFile) = @_;

    my @gr = `grep "^i" $goFile | grep bestTerm`;
    ##die Dumper(\@gr);
    return @gr > 0 && length($gr[0]) > 2;
}

# assume each line is a cluster, with members in a tab-delimited list
# modified to use 1-based indexing
sub clusterFromEachLine {
    my ($ret, $in, $minCluster) = @_;
    $minCluster //= 3;

    open my $IN, "<", $in or die "can't read from $in. $!";
    my $i=1;
    while (<$IN>) {
	next if /^#/;
	chomp;
	my @cluster = split /\t/;
	next if @cluster < $minCluster;
	$ret->{$_}{$i} = 1 for @cluster;
	$i++;
    }

    return $i;
}

# figure out which cluster this edge is in, then give the name for that GO term
sub nameCluster {
    my ($p1, $p2, $clusters, $clusterGO, $goNameMap) = @_;

    my @cluster1 = keys %{ $clusters->{$p1} };
    my @matches;
    for my $cid (@cluster1) {
	push @matches, $cid if exists $clusters->{$p2}{$cid};
    }
    die "too many matches for $p1 $p2" if @matches > 1;
    if (0 == @matches) {
	return "-1";
    }

    my $cid = pop @matches;
    my $go = $clusterGO->{$cid};
    my $name = $go;
    if ($go ne 'none') {
	$name = $goNameMap->{$go}{name} // 
	    die "can't find goName for '$go' ($p1 $p2 $cid)";
    }
    return "$cid $name";
}
