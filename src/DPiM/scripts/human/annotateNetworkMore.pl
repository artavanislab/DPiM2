#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Storable qw(retrieve);
use HomeBrew::IO qw(checkExist readColsRef);
use DpimLib qw(readHS networkHashFromEdgeList);

## generate table with 6+ columns:
##   node1 node2 symbol1 symbol2 score source 
## additional columns are the annotated term in different MCL parameters
## source networks: DPIM1 allRep bestRep CompPASS humanHG_allRep

my %opts = getCommandLineOptions();
$opts{human} = 1;

{
    my $netFile = $opts{net};
    my $annotFile = $opts{annot};
    my $out = $opts{out};
    my $prevNetFile = $opts{prevnet};
    
    # $annot{$fbgn} = {
    #   flySymbol=>X, humanSymbol=>X, clusterI=>$i, clusterName=>$clusterName,
    #   go { goID1 => goName1, goID2=> goName2,...}
    #   reactome { reID1 => reName1, reID2=> reName2,...}
    #   pfam { goID1 => goName1, goID2=> goName2,...}
    #   fb { goID1 => goName1, goID2=> goName2,...}
    #   ot { goID1 => goName1, goID2=> goName2,...}
    #   mb { goID1 => goName1, goID2=> goName2,...}
    # };
    my %annot;
    readAnnot(\%annot, $annotFile);
    
    my %prevNet;
    networkHashFromEdgeList(\%prevNet, $prevNetFile, undef, 'symmetric', undef, 1);

    my @cols = qw( protein1 protein2 score BioPlex2.0
                   humanSymbol1 humanSymbol2 clusterI clusterName
                   sharedGO sharedReactome sharedPfam sharedFlyBaseDisease sharedOpenTargetsDisease sharedMetaBaseDisease);
    
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT join "\t", @cols;

    open my $IN, "<", $netFile or die "Can't read from $netFile: $!";
    while (my $line = readHS($IN, $opts{human})) {
	chomp $line;
	my ($fbgn1, $fbgn2) = (split /\s+/, $line)[0..1];

	die "can't find $fbgn1" if ! exists $annot{$fbgn1};
	die "can't find $fbgn2" if ! exists $annot{$fbgn2};
	
	#my $flySymbol1 = $annot{$fbgn1}{Symbol};
	#my $flySymbol2 = $annot{$fbgn2}{Symbol};

	my $prevNet = isTrue(exists $prevNet{$fbgn1}{$fbgn2});

	my $humanSymbol1 = $annot{$fbgn1}{Symbol};
	my $humanSymbol2 = $annot{$fbgn2}{Symbol};
	
	my ($clusterI, $clusterName) = qw(NA NA);
	if ($annot{$fbgn1}{clusterI} eq $annot{$fbgn2}{clusterI}) {
	    $clusterI = $annot{$fbgn1}{clusterI};
	    $clusterName = $annot{$fbgn1}{clusterName};
	}

	my $goString = termString($annot{$fbgn1}{go}, $annot{$fbgn2}{go},);
	my $reString = termString($annot{$fbgn1}{reactome}, 
				  $annot{$fbgn2}{reactome});
	my $pfString = termString($annot{$fbgn1}{pfam}, 
				  $annot{$fbgn2}{pfam});
    my $fbString = termString($annot{$fbgn1}{fb},
                  $annot{$fbgn2}{fb});
    my $otString = termString($annot{$fbgn1}{ot},
                  $annot{$fbgn2}{ot});
    my $mbString = termString($annot{$fbgn1}{mb},
                  $annot{$fbgn2}{mb});
	
	$line.="\t";
	$line.= join("\t", $prevNet, $humanSymbol1, 
		     $humanSymbol2, $clusterI, $clusterName, $goString, 
		     $reString, $pfString, $fbString, $otString, $mbString);
	say $OUT $line;
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
	prevnet => '/home/kli3/proj/Interactome/data/BioPlex2.0Nature/BioPlex_interactionList_v4a_slim.tsv',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -net nrBait.net -annot nodeAnnot.tsv -out output < ".
	"$defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "net=s", "annot=s", "out=s", "prevnet=s");
    die $usage unless exists $opts{net} && exists $opts{annot} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});
    checkExist('f', $opts{annot});
    checkExist('f', $opts{prevnet});

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

# $annot{$fbgn} = {
#   flySymbol=>X, humanSymbol=>X, clusterI=>$i, clusterName=>$clusterName,
#   go { goID1 => goName1, goID2=> goName2,...}
#   reactome { reID1 => reName1, reID2=> reName2,...}
#   pfam { goID1 => goName1, goID2=> goName2,...}
#   fb { goID1 => goName1, goID2=> goName2,...}
#   ot { goID1 => goName1, goID2=> goName2,...}
#   mb { goID1 => goName1, goID2=> goName2,...}
# };
sub readAnnot {
    my ($ret, $in) = @_;

    my @cols = qw(EntrezID Symbol clusterI clusterName Symbol
                  goID goName reactomeID reactomeName pfamID  pfamName flybaseDOID flybaseDOIDName OpenTargetsDiseaseID OpenTargetsDiseaseName MetaBaseDiseaseID MetaBaseDiseaseName);

    my @read;
    readColsRef(\@read, $in, \@cols, undef, "\t");

    for my $row (@read) {
	my ($EntrezID, $Symbol, $clusterI, $clusterName, $humanSymbol, 
	    $goIDString, $goNameString, $reIDString, $reNameString, 
	    $pfIDString, $pfNameString, $fbIDString, $fbNameString, $otIDString, $otNameString, $mbIDString, $mbNameString) =
		map { $row->{$_} } @cols;
	my %go = termHash($goIDString, $goNameString);
	my %reactome = termHash($reIDString, $reNameString);
	#my %pfam =  termHash($pfIDString, $pfNameString);
	#my %fb =  termHash($fbIDString, $fbNameString);
	my %ot =  termHash($otIDString, $otNameString);
	my %mb =  termHash($mbIDString, $mbNameString);
	

	#$ret->{$EntrezID} = { EntrezID=>$EntrezID, humanSymbol=>$Symbol,
	#		  clusterI=>$clusterI, clusterName=>$clusterName,
	#		  reactome => \%reactome, go => \%go, pfam => \%pfam, fb => \%fb, ot => \%ot, mb => \%mb};
    #}
	$ret->{$EntrezID} = { Symbol=>$Symbol, humanSymbol=>$Symbol,
              clusterI=>$clusterI, clusterName=>$clusterName,
              reactome => \%reactome, go => \%go, ot => \%ot, mb => \%mb};
    }

    return;
}

sub termHash {
    my ($idString, $nameString) = @_;

    my @keys = split /,/, $idString;
    my @values = split /','/, $nameString;
    $values[ 0] =~ s/^'//;
    $values[-1] =~ s/'$//;
    die "different number of keys (".(0+ @keys).") and values (".
	(0+ @values).") ".Dumper(\@keys, \@values) if @keys != @values;

    my %ret = map { $keys[$_] => $values[$_] } 0..$#keys;
    delete $ret{NA} if exists $ret{NA};
    return %ret;
}

sub termString {
    my ($termHash1, $termHash2) = @_;
    
    my @sharedTerms;
    for my $term (keys %$termHash1) {
	if (exists $termHash2->{$term}) {
	    push @sharedTerms, $termHash1->{$term};
	}
    }

    if (0 == @sharedTerms) {
	return "NA";
    } else {
	return join ",", map { "\"$_\"" } @sharedTerms;
    }
}

sub isTrue {
    my $arg = shift;
    return "1" if $arg;
    return "0";
}
