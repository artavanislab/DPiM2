#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HTTP::Tiny;
use JSON;
use HomeBrew::IO qw(checkExist readColsRef writeCols);
#use DpimLib qw(getLineAPMS);

# make a map from (human) entrez ids to FBgn


# turn corum clusters listed in allComplexes.csv into something I can use

my %opts = getCommandLineOptions();

#Complex id;Complex name;Synonyms;organism;subunits (UniProt IDs);subunits (Entrez IDs);protein complex purification method;PubMed id;FunCat categories;"functional comment";"disease comment";"subunit comment";
{
    my $corumFile = $opts{corum};
    my $flyFile = $opts{flymap};
    my $out = $opts{out};
    my $minScore = $opts{minscore};

    my @corumIDs;
    parseCorum(\@corumIDs, $corumFile);
    @corumIDs = sort {$a <=> $b} @corumIDs;

    my %flyMap;
    readDIOPT(\%flyMap, $flyFile);
    my ($totalUniq, $totalAmbig) = topScoreOnly(\%flyMap) 
	unless $opts{allmatches};
    
    my @cols = qw(fbgn score flySymbol humanSymbol source);
    my @orthos; # orthos[i] = {entrez=> $entrez, map {$col => $val}}
    my ($uniq, $ambig, $noMatch) = (0, 0, 0);
    for my $entrez (@corumIDs) {
	next if $entrez==0;
	my $matches = $flyMap{$entrez};
	if (! defined $matches) {
	    my %row = map { $_ => 'NA' } @cols;
	    $row{entrez} = $entrez;
	    $noMatch++;
	    next;
	} elsif( @$matches == 1) {
	    $uniq++;
	} else {
	    $ambig++;
	}
	for my $match (@$matches) {
	    next if $match->{score} < $minScore;
	    my %row = %$match;
	    $row{entrez} = $entrez;
	    push @orthos, \%row;
	}
    }

    unshift @cols, 'entrez';
    my @data;
    for my $col (@cols) {
	push @data, [ map { $_->{$col} } @orthos ];
    }
    
    my $header = join "\t", @cols;
    my $format = join "\t", qw(%d %s %f %s %s %s);
    my $preComments = "# $uniq unique matches and $ambig ambiguous matches\n";
    $preComments .= "# $noMatch genes had no match\n";
    writeCols($out, \@data, $header, $format, $preComments);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	flymap => '/home/glocke/DPiM/corum/DIOPT_Translations_for_George.txt',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -corum corumClusters -out output < $defaultString ".
	" -allmatches -minscore none>\n";

    my %opts = ();
    GetOptions(\%opts, "corum=s", "out=s", "flymap=s", "allmatches", 
	       "minscore=f");
    die $usage unless exists $opts{corum} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{corum});
    checkExist('f', $opts{flymap});

    return %opts;
}


# It turns out, we don't have to do this map because the DIOPT table
# already has human entrez id's
# parse the corum file to learn which id's to translate
# then use rest.ensembl.org to map entrez id's to gene symbol
# use values obtained from the EntrezGene database
sub mapEntrez2Symbol {
    my ($ret, $corumFile) = @_;

    my %allEntrez; # allEntrez{$entrezID} = 1
    
    parseCorum(\%allEntrez, $corumFile);

    ## debug null response $allEntrez{123456789} = 1;

    my $submit = submitter();
    for my $id (keys %allEntrez) {
	my %hash = $submit->($id);
	#die Dumper(\%hash);
	die "null response for $id" unless exists $hash{EntrezGene};
	$ret->{$id} = { symbol => $hash{EntrezGene}{display_id},
			synonyms => $hash{EntrezGene}{synonyms} };
    }

    return;
}

# just find the entrez id's in the file
sub parseCorum {
    my ($ret, $inFile) = @_;

    open my $IN, "<", $inFile or die "Can't read from $inFile. $!";
    
    my %allEntrez;

    my $speciesCol = 3;
    my $entrezCol = 5;
    <$IN>; # ignore header
    while (<$IN>) {
	my @spl = split /;/;
	next unless $spl[$speciesCol] =~ /uman/ || 
	    $spl[$speciesCol] =~ /Mammalia/;
	#next if $spl[$entrezCol] =~ /\(/;
	my @ids = split /,/, $spl[$entrezCol];
	$_ =~ s/\(//g for @ids;
	$_ =~ s/\)//g for @ids;
	$allEntrez{$_}=1 for @ids;
	#last if 3 < keys %allEntrez;
    }
    
    push @$ret, $_ for keys %allEntrez;

    return;
}


# return a closure that interfaces with rest.ensembl.org
# the only reason to use a closure is just to wrap up the various $http, 
#   $server, etc. arguments, and keep them from cluttering the main body of code
sub submitter {
    my $http = HTTP::Tiny->new();
    my $server = 'http://rest.ensembl.org';
    my $ext = '/xrefs/name/human/';
    my %options = (headers => { 'Content-type' => 'application/json',
				'external_db' => 'EntrezGene'});

    return sub {
	my ($id) = @_;
	my $response = $http->get($server.$ext.$id.'?', \%options);	
	die "Failed to access $id http\n" unless $response->{success};
	
	my $arr = decode_json($response->{content});
	return map { $_->{dbname} => $_ } @$arr;
	if(length $response->{content}) {
	    my $hash = decode_json($response->{content});
	    local $Data::Dumper::Terse = 1;
	    local $Data::Dumper::Indent = 1;
	    print Dumper $hash;
	    print "\n";
	}
    }
}
 
# populate the flyMap hash
# flyMap{entrez} = [{ fbgn=>$fbgn, score=>$score, flySymbol=>$tla1, 
#                     humanSymbol => $TLA1, source=$databases }, {}, {}...]
sub readDIOPT {
    my ($flyMap, $flyFile) = @_;
    
    my @cols = ('FBgn Search Term', 
		'FlyBaseID',
		'Human GeneID',
		'Human Symbol',
		'Fly Symbol',
		'Weighted Score',
		'Prediction Derived From');
    my @read;
    readColsRef(\@read, $flyFile, \@cols, undef, "\t");
    
    for my $row (@read) {
	next unless length($row->{'Human GeneID'});
	my $id = $row->{'Human GeneID'};
	$flyMap->{$id} //= [];
	push @{$flyMap->{$id}}, { fbgn => $row->{'FlyBaseID'},
				  altfbgn=> $row->{'FBgn Search Term'}, 
				  humanSymbol => $row->{'Human Symbol'},
				  flySymbol => $row->{'Fly Symbol'},
				  score=> $row->{'Weighted Score'}, 
				  source=> $row->{'Prediction Derived From'}, 
	};
    }
    
    for my $list (values %$flyMap) {
	$list = [sort {$b->{score} <=> $a->{score}} @$list];
    }

    return;
}

# remove any match that has less than the top score for that gene
sub topScoreOnly {
    my ($flyMap) = @_;

    my ($uniq, $ambig) = (0, 0);
    for my $match (values %$flyMap) {
	my $topScore = $match->[0]{score};
	pop @$match while $match->[-1]{score} != $topScore;
	if (@$match == 1) {
	    $uniq++;
	} else {
	    $ambig++;
	}
    }
    
    return($uniq, $ambig);
}
