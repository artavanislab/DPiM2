#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(max);
use File::Temp qw(tempfile);
use Data::Dumper;
use File::Slurp;
use HomeBrew::IO qw(checkExist readColsRef readList2 writeCols);
#use DpimLib qw(getLineDP4APMS);

# for all nodes in a network, annotate the following:
# * gene symbol
# * whether it was in DPIM
# * cluster membership
# * human orthologs (HGNC symbols)
# * all terms

my %opts = getCommandLineOptions();

{
    my $clustFile = $opts{in};
    my $clustAnnotFile = $opts{clusterannot};
    my $out = $opts{out};
    my $dioptFile = $opts{diopt};
    my $freqFile = $opts{freq};
    my $termFile = $opts{terms};
    my $termScr = $opts{termscr};
    my $transFile = $opts{trans};

    # $cluster{$fbgn} = { clustID => cluster number, term => mostEnrichedTerm,
    #                     sig => true/false, name => '$clustID. ! $term' }
    #  the "!" in the name is set if sig is false         ^
    my %clusters;
    if (exists $opts{mcl}) {
	readMCLCluster(\%clusters, $clustFile, $clustAnnotFile);
    } else {
	readCluster(\%clusters, $clustFile, $clustAnnotFile);
    }
    my %diopt;
    readDiopt(\%diopt, $dioptFile, \%clusters);

    # freq{$fbgn} = { freq, avgPep };
    my %freq;
    readFreq(\%freq, $freqFile);
    
    if (! defined $termFile) {
	my $FH;
	($FH, $termFile) = tempfile();
	close $FH;
	my $cmd = "$termScr $clustFile $termFile";
	say $cmd;
	system($cmd);
    }
    # $ret->{$fbgn}{GO}{names => "'name1','name2',.. ids => "GO:123,GO:124,.."}
    # $ret->{$fbgn}{Reactome}{names => "'name1','name2',.. ids => "GO:123,GO:124,.."}
    my %terms;
    readTerms(\%terms, $termFile, \%clusters);

    my %fbgnMap; ## fbgn2genesymbol
    say "parsing $transFile...";
    makeMap(\%fbgnMap, $transFile);

    # * gene symbol
    # * cluster membership
    # * human orthologs (HGNC symbols)
    # * all terms
    my @fbgn = sort keys %clusters;
    # KJ commented this die line to allow no match found for FBgn999_ cases 01/31/2018
    #my @symbol = map { $fbgnMap{$_} // die "can't get symbol for $_" } @fbgn;
    my @symbol = map { $fbgnMap{$_} } @fbgn;
    my @freq = map { $freq{$_}{freq} } @fbgn;
    my @avgPep = map { $freq{$_}{avgPep} } @fbgn;
    my @clusterI = map { $clusters{$_}{clustID} } @fbgn;
    my @clusterName = map { $clusters{$_}{name} } @fbgn;
    my @orthos = map { $diopt{$_} } @fbgn;
    my @goIDs = map { $terms{$_}{GO}{ids} } @fbgn;
    my @goNames = map { $terms{$_}{GO}{names} } @fbgn;
    my @reIDs = map { $terms{$_}{Reactome}{ids} } @fbgn;
    my @reNames = map { $terms{$_}{Reactome}{names} } @fbgn;
    my @pfIDs = map { $terms{$_}{Pfam}{ids} } @fbgn;
    my @pfNames = map { $terms{$_}{Pfam}{names} } @fbgn;
    my @fbIDs = map { $terms{$_}{FlyBase}{ids} } @fbgn;
    my @fbNames = map { $terms{$_}{FlyBase}{names} } @fbgn;
    my @otIDs = map { $terms{$_}{OpenTargets}{ids} } @fbgn;
    my @otNames = map { $terms{$_}{OpenTargets}{names} } @fbgn;
    my @mbIDs = map { $terms{$_}{MetaBase}{ids} } @fbgn;
    my @mbNames = map { $terms{$_}{MetaBase}{names} } @fbgn;


    my $format = join "\t", ("%s") x 19;
    my $header = join "\t", qw(fbgn flySymbol freq avgPep clusterI clusterName 
                               humanSymbol goID goName reactomeID reactomeName 
                               pfamID pfamName flybaseDOID flybaseDOIDName OpenTargetsDiseaseID OpenTargetsDiseaseName MetaBaseDiseaseID MetaBaseDiseaseName );
    my @data = (\@fbgn, \@symbol, \@freq, \@avgPep, \@clusterI, \@clusterName, 
		\@orthos, \@goIDs, \@goNames, \@reIDs, \@reNames, \@pfIDs, 
		\@pfNames, \@fbIDs, \@fbNames, \@otIDs, \@otNames, \@mbIDs, \@mbNames);
    my $preComments = "# $0 annotated all nodes in $clustFile";
    writeCols($out, \@data, $header, $format, $preComments);
}

exit;
    
   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	clusterannot => '/home/glocke/DPiM/augRemap/nrBait11-08-2016/mcl/mcl.clusters/mcl.clusters.nrBait.ppi.abc-format.i4.txt.R_notEASE.GOTest',
	diopt => '/home/glocke/DPiM/corum/DIOPT_Translations_for_George_110916.txt',
	freq => '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.rmDupes.pepFDR.sumIso.applyLC.untrans.tscCut.rebrand.protFreq',
	termscr => $ENV{DPSCR}."/rscripts/annotateNodes.R",
	trans => '/home/glocke/DPiM/flybase/fbgn_annotation_ID_fb_2015_04.tsv',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in in.mcl.txt -out output < -terms ".
	"annotateNodes.R.out $defaultString -mcl >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "terms=s", "clusterannot=s", "diopt=s",
	       "freq=s", "mindiopt=f", "termscr=s", "trans=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{clusterannot});
    checkExist('f', $opts{diopt});
    checkExist('f', $opts{trans});
    checkExist('f', $opts{terms}) if exists $opts{terms};
    checkExist('f', $opts{termscr}) if ! exists $opts{terms};

    return %opts;
}

# $cluster{$fbgn} = { clustID => cluster number, term => mostEnrichedTerm,
#                     sig => true/false, name => '$clustID. ! $term' }
#  the "!" in the name is set if sig is false         ^
sub readCluster {
    my ($ret, $clustFile, $clustAnnotFile) = @_;

    my @rawClust = readList2($clustFile);
    my %clust2Fbgn; ## clustID{$clustID}{$fbgn}=1 iff this cluster contains $fbgn
    for my $clust (@rawClust) {
	my $clustID = shift @$clust;
	$clust2Fbgn{$clustID}{$_}=1 for @$clust;
    }

    my @cols = qw(name sig clustID);
    my @read;
    readColsRef(\@read, $clustAnnotFile, \@cols, undef, "\t");
    for my $row (@read) {
	my $clustID = $row->{clustID};
	my $term = $row->{name};
	$term =~ s/"//g;
	my $sig = $row->{sig};
	if ($sig eq 'FALSE') {
	    $sig = undef;
	} else {
	    $sig = 1;
	}
	my $name = "'$clustID. ";
	$name.= "!" if !$sig;
	$name.= $term."'";
	
	my $val = {clustID=>$clustID, term=>$term, sig=>$sig, name=>$name};
	$ret->{$_} = $val for keys %{ $clust2Fbgn{$clustID} };
    }

    return;
}

# $cluster{$fbgn} = { i => cluster number, term => mostEnrichedTerm,
#                     sig => true/false, name => '$clustID. ! $term' }
#  the "!" in the name is set if sig is false         ^
sub readMCLCluster {
    my ($ret, $clustFile, $clustAnnotFile) = @_;

    my @clustStrings = read_file($clustFile);
    my %clust2Fbgn; ## clustID{$clustID}{$fbgn}=1 iff this cluster contains $fbgn
    for my $i (0..$#clustStrings) {
	my @fbgns = split /\s+/, $clustStrings[$i];
	$clust2Fbgn{$i+1}{$_}=1 for @fbgns;
    }

    my @cols = qw(name sig i);
    my @read;
    #readColsRef(\@read, $clustAnnotFile, \@cols, undef, "\t", "checkQuotes");
    readColsRef(\@read, $clustAnnotFile, \@cols, undef, "\t");
    for my $row (@read) {
	my $clustID = $row->{i};
	my $term = $row->{name};
	$term =~ s/"//g;
	my $sig = $row->{sig};
	if ($sig eq 'FALSE') {
	    $sig = undef;
	} else {
	    $sig = 1;
	}
	my $name = "'$clustID. ";
	$name.= "!" if !$sig;
	$name.= $term."'";
	
	my $val = {clustID=>$clustID, term=>$term, sig=>$sig, name=>$name};
	$ret->{$_} = $val for keys %{ $clust2Fbgn{$clustID} };
    }

    return;
}

# $ret->{$fbgn} = "TLA1,TLA2,..."
sub readDiopt {
    my ($ret, $in, $nodes) = @_;

    my %ranks = ( high => 2, moderate => 1 );
    
    my @cols = ("FlyBaseID", "Fly GeneID", "Human Symbol", "Rank");
    my @read;
    readColsRef(\@read, $in, \@cols, undef, "\t");
    my %allOrthos;
    for my $row (@read) {
	my $fbgn = $row->{"FlyBaseID"};
	next unless exists $nodes->{$fbgn};
	my $rnk = $ranks{$row->{"Rank"}} // next;
	
	$allOrthos{$fbgn}{$row->{"Human Symbol"}} //= $rnk;
	$allOrthos{$fbgn}{$row->{"Human Symbol"}} = 
	    max($rnk, $allOrthos{$fbgn}{$row->{"Human Symbol"}});
    }

    for my $fbgn (keys %$nodes) {
	if (exists $allOrthos{$fbgn}) {
	    $ret->{$fbgn} = join ","
		, sort {$allOrthos{$fbgn}{$b} <=> $allOrthos{$fbgn}{$a}}
	    sort keys %{ $allOrthos{$fbgn} };
	} else {
	    $ret->{$fbgn} = 'NoMatch'
	}
    }
	
    return;
}

# ret{$fbgn} = { freq, avgPep };
sub readFreq {
    my ($ret, $freqFile) = @_;

    my @cols = qw(fbgn    Fraction        avg_tot_pep);
    my @read;
    readColsRef(\@read, $freqFile, \@cols);
    for my $row (@read) {
	$ret->{$row->{fbgn}} = {
	    freq => $row->{Fraction},
	    avgPep => $row->{avg_tot_pep},
	};
    }

    return;
}

# $ret->{$fbgn}{GO}{names => "'name1','name2',.. ids => "GO:123,GO:124,.."}
# $ret->{$fbgn}{Reactome}{names => "'name1','name2',.. ids => "GO:123,GO:124,.."}
sub readTerms {
    my ($ret, $termFile, $nodes) = @_;

    my @cols = qw(FLYBASE ID name);
    my @read;
    readColsRef(\@read, $termFile, \@cols, undef, "\t");
    my @DBs = qw(GO Reactome Pfam FlyBase OpenTargets MetaBase);
    my $termDB;
    my %terms;
    for my $row (@read) {
	my $fbgn = $row->{FLYBASE};
	next unless exists $nodes->{$fbgn};
	if ($row->{ID} =~ /^GO/) {
	    $termDB = 'GO';
	} elsif ($row->{ID} =~ /^RE/) {
	    $termDB = 'Reactome';
	} elsif ($row->{ID} =~ /^PF/) {
	    $termDB = 'Pfam';
	} elsif ($row->{ID} =~ /^DOID/) {
		$termDB = 'FlyBase';
	} elsif ($row->{ID} =~ /^OpenTargets/) {
        $termDB = 'OpenTargets';
	} elsif ($row->{ID} =~ /^MetabaseID/) {
        $termDB = 'MetaBase';
	} else {
	    die "can't detect what database this term is from\n".Dumper($row);
	}
	$terms{$fbgn}{$termDB}{$row->{ID}} = $row->{name};
    }

    
    my $null = { names => 'NA', ids => 'NA' };
    for my $fbgn (keys %$nodes) {
	for my $termDB (@DBs) {
	    if (exists $terms{$fbgn}{$termDB}) {
		my @ids = keys %{ $terms{$fbgn}{$termDB} };
		my @names = map { "'$_'" } values %{ $terms{$fbgn}{$termDB} };
		$ret->{$fbgn}{$termDB} = {
		    names => (join ",", @names),
		    ids => (join ",", @ids),
		};
	    } else {
		$ret->{$fbgn}{$termDB} = $null;
	    }
	}
    }
    
    return;
}

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
