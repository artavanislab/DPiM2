#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min);
use Data::Dumper;
use Storable qw(retrieve);
use HomeBrew::IO qw(checkExist readColsHash readColsRef writeCols);
use Statistics::R;
use DpimLib qw(readGoDB);

# test if clusters have enriched GO terms
# what to report?
#  - for each cluster, min hypgeo p value, best go term
#  - for each cluster, what GO terms?
#
# something confusing: **minP is actually minQ**

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $goHistFile = $opts{gohist};
    my $out = $opts{out};
    my $mode = $opts{mode};
    my $alpha = $opts{alpha};
    my $excludeString = $opts{exclude};
    my $goDBFile = $opts{godb};
    my $goNameFile = $opts{goname};
    my $maxPercent = $opts{maxpercent};
    my $minCluster = $opts{mincluster};
    
    die "mcl is the only mode implemented" if $mode ne 'mcl';
    
    my $nProtein = findNProtein($goHistFile);
    my %globalGO = readColsHash($goHistFile, [qw(term count)]);
    # prune go terms with a single count
    {
	my @singletons = grep { 1 == $globalGO{$_} } keys %globalGO;
	#die Dumper(\@singletons);
	delete $globalGO{$_} for @singletons;
    }
    my %droppedTerms;
    if ($maxPercent > 0) {
	# remove terms that cover massive portions of the network as they are
	#   generally meaningless for cluster identity
	my $max = $nProtein * $maxPercent / 100;
	my @vague = grep { $max < $globalGO{$_} } keys %globalGO;
	delete $globalGO{$_} for @vague;
	$droppedTerms{$_} = 1 for @vague;
    }
    if ($excludeString ne 'none') {
	my @excl = split /,/, $excludeString;
	$droppedTerms{$_} = 1 for @excl;
    }
    my $nTerms = 0+ keys %globalGO;
    #say $nTerms;
    
    my %goDB;
    # goDB{proteinID} => [ term1, term2,... ]
    if (exists $opts{human}) {
	readGoTab(\%goDB, $goDBFile);
    } else {
	readGoDB(\%goDB, $goDBFile);
    }
    
    my $goNameMap = retrieve($goNameFile);
    ##$goNameMap->{none} = { name => '"NA"' };

    # LATER: switch on modes to read other cluster formats
    my @clusters = clusterFromEachLine($in, $minCluster);
    
    my @stats = map { findTerms($_, \%goDB, $goNameMap, \%droppedTerms) } @clusters;
    $stats[$_]{clusterID} = $_+1 for 0..$#stats;
    #stats[i] = {size=> size of cluster, 
    #            terms => {term1=> cnt, term2=> cnt} 
    #            nterms => [number of terms with cnt>1]
    #            clusterID => numerical cluster id}
    
    my $R = Statistics::R->new();
    $R->startR;

    # Main Loop:
    # add stats[i]{minP=1, term='none', test => {term1=>p, term2=>p}}
    my $sigSum = 0; # how many have minP < 0
    my @cols = qw(term minP bestCnt nWayTie sig);
    for my $cl (@stats) {
	if ($cl->{nterms} == 0) {
	    $cl->{term} = 'none';
	    $cl->{minP} = 1;
	    $cl->{bestCnt} = 0;
	    $cl->{nWayTie} = 0;
	    $cl->{sig} = 0;
	    $cl->{name} = 'NA';
	    next;
	}

	my %goCnts;
	my $pulls = $cl->{size};
	for my $term (keys %{ $cl->{terms} }) {
	    my $found = $cl->{terms}{$term};
	    ## exclude term unless more than 1 node in cluster is annotated
	    next if 1 == $found;
	    my $white = $globalGO{$term} // die "can't find $term in global??";
	    my $black = $nProtein - $white;
	    $goCnts{$term} = { found => $found, white => $white, 
			       black=> $black };
	    ##$tests{$term} = 1 - gsl_cdf_hypergeometric_P($found, $white, 
	    ##						 $black, $pulls);
	}
	my @terms = enrich($R, \%goCnts, $pulls, $nTerms);
	$cl->{$cols[$_]} = $terms[$_] for  0..$#terms;

	my $sig = ($cl->{minP} <= $alpha)?1:0;
	$sigSum+=$sig;
	$cl->{sig} = $sig;
	$cl->{name} = '"'.$goNameMap->{$cl->{term}}{name}.'"' // 'NA';
    }
    my $sigPercent = sprintf("%.1f", 100 * $sigSum / @stats);
    
    my @reportCols = qw(minP term sig size bestCnt nWayTie name);
    my @data = ();
    for my $col (@reportCols) {
	push @data, [map { $_->{$col} } @stats];
    }
    unshift @reportCols, 'i';
    my $header = join "\t", @reportCols;
    my $format = join "\t", qw(%d %.3e %s %d %d %d %d %s);
    my $preComments = "# compared $in against $goHistFile using $goDBFile and $goNameFile\n";
    $preComments.= "# $sigPercent% of clusters beat alpha = $alpha\n";
    $preComments.= "# performed p-value corrections for $nTerms tests\n";
    $preComments.= "# excluded terms covering more than $maxPercent of nodes\n" 
	if $maxPercent > 0;
    $preComments.= "# following goids were excluded: ".
	(join ", ", keys %droppedTerms)."\n" if 0 < keys %droppedTerms;
    writeCols($out, \@data, $header, $format, $preComments, 1);
    say $sigPercent;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    # mode:
    my @modes = qw(mcl mode2);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my @excludedTerms = qw( GO:0003674 GO:0005515 GO:0005575 GO:0005634 
GO:0005737 GO:0005739 GO:0005829 GO:0008150 GO:0005829 GO:0008150 );
    my $exDefault = join ",", @excludedTerms;
    ## that's molecular_function, protein binding, cellular component, nucleus, 
    ## cytoplasm, mitochondrion, cytosol, biological_process
    
    my %defaults = (
	alpha => 0.05,
	exclude => $exDefault,
	godb => '/home/glocke/DPiM/flybase/gene_association.fb',
	goname => '/home/glocke/DPiM/flybase/goMap.full.storable',
	maxpercent => -1,
	mincluster => 3,
	);
    my $defaultString = join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;
    my $usage = "usage: $0 -in cluster.list -gohist wholeNet.nodewise.gohist ".
	"-out output < $modeString $defaultString -human -reportall >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "gohist=s", "out=s", "mode=s", "alpha=f",
	       "exclude=s", "godb=s", "goname=s", "maxpercent=f", 
	       "mincluster=i", "human");
    die $usage unless exists $opts{in} &&  exists $opts{gohist} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{gohist});
    checkExist('f', $opts{godb});
    checkExist('f', $opts{goname});

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    return %opts;
}

# ret->{proteinID} => [ term1, term2,... ]
sub readGoTab {
    my ($ret, $in, $ignore) = @_;

    $ignore //= {};
    
    my %goHash;

    my @read;
    readColsRef(\@read, $in, [qw(entrez term)]);
    
    for my $row (@read) {
	my ($prot, $term) = ($row->{entrez}, $row->{term});
	next if exists $ignore->{$term};
	$goHash{$prot}{$term} = 1;
    }

    for my $p (keys %goHash) {
	$ret->{$p} = [ keys %{ $goHash{$p} } ];
    }
    
}

sub findNProtein {
    my ($in) = @_;

    open my $IN, "<", $in or die "Can't read from $in. $!";

    do { $_ = <$IN> } until ( eof($IN) || /total proteins found = (\d+)/ ||
			      /total tests = (\d+)/);
    die "can't find nProtein in $in." if eof($IN) || !defined $1 || $1 < 1;

    return 0+ $1;
}

# ret[i] = all fbgns in i'th row
sub clusterFromEachLine {
    my ($in, $minCluster) = @_;
    $minCluster //= 3;

    my @ret;
    
    open my $IN, "<", $in or die "can't read from $in. $!";
    while (<$IN>) {
	next if /^#/;
	chomp;
	my @cluster = split /\t/;
	next if @cluster < $minCluster;
	push @ret, \@cluster;
    }

    return @ret;
}

# ret = {size=> size of cluster, terms => {term1=> cnt, term2=> cnt}, nterms=>2}
sub findTerms {
    my ($cluster, $goDB, $goNameMap, $droppedTerms) = @_;

    my %ret;
    $ret{size} = 0+ @$cluster;

    my %terms;
    for my $prot (@$cluster) {
	next unless exists $goDB->{$prot};
	for my $term (@{ $goDB->{$prot} }) {
	    next if exists $droppedTerms->{$term};
	    $terms{$term}++;
	    #say "$prot -> GO:0031532" if $term eq 'GO:0031532';
	}
    }

    # standardize GO terms
    my %oldNew;
    for my $go (keys %terms) {
	if (exists $goNameMap->{$go}{standard_id} &&
	    $goNameMap->{$go}{standard_id} ne $go) {
	    $oldNew{$go} = $goNameMap->{$go}{standard_id};
	    die "standard equals self" if $oldNew{$go} eq $go;
	}
    }
    for my $oldGo (keys %oldNew) {
	my $newGo = $oldNew{$oldGo};
	$terms{$newGo} += $terms{$oldGo};
	delete $terms{$oldGo};
    }

    $ret{terms} = \%terms;
    $ret{nterms} = 0+ grep { 1 < $_ } values %terms;
    return \%ret;
}

# perform hypergeomtric tests on the given terms
sub enrich {
    my ($R, $goCnts, $pulls, $nTerms) = @_;
    
    my @terms = keys %$goCnts;
    my (@found, @white, @black);
    for my $t (@terms) {
	push @found, $goCnts->{$t}{found};
	push @white, $goCnts->{$t}{white};
	push @black, $goCnts->{$t}{black};
    }
    {
	# sort lists according to @found
	my @i = sort {$found[$a] <=> $found[$b]} 0..$#found;
	@terms = @terms[@i];
	@found = @found[@i];
	@white = @white[@i];
	@black = @black[@i];
    }
    my $foundString = join ",", @found;
    my $whiteString = join ",", @white;
    my $blackString = join ",", @black;

    $R->send(qq'found <- c($foundString)');
    $R->send(qq'white <- c($whiteString)');
    $R->send(qq'black <- c($blackString)');
    $R->send(qq'hg <- phyper(found-1, white, black, $pulls, lower.tail=F)');
    $R->send(qq'paste(hg, collapse=",")');
    $R->send(qq'pp <- p.adjust(hg, method="holm", n=$nTerms)');
    # note that setting n=$nTerms is equivalent to adding 1's untl length(hg)==$nTerms
    $R->send(qq'paste(pp, collapse=",")');

    my $pString = $R->read;
    my $x = $pString;
    $pString =~ s/\[1\] //;
    $pString =~ s/"//g;
    my @q = split ",", $pString;
    my @i = sort { $q[$a] <=> $q[$b] } 0..$#terms;
    @q = @q[@i];
    
    ##die Dumper(\@q, "nTerms = ".(0+ @black)."; size of q = ".(0+ @q));
    
    my @cols = qw(term minP bestCnt nWayTie sig);
    my $bestI = minIndex(\@q);
    my $term = $terms[$bestI];
    my $minP = $q[$bestI];
    my $bestCnt = $goCnts->{$term};
    my $nWayTie = 0+ grep { $_ == $minP } @q;
    if ($nWayTie > 1) {
	my @allBest = grep { $q[$_] == $minP } 0..$#q;
	for my $i (@allBest) {
	    next if ! defined $terms[$i];
	    my $t = $terms[$i];
	    if ($bestCnt < $goCnts->{$t}) {
		$term = $t;
		$bestCnt = $goCnts->{$t};
	    }
	}
    }

    #if ($nWayTie > 1) {
    if (0) { ## debug
	say $foundString;
	say $whiteString;
	say $blackString;
	say $nTerms;
	say 0+ @terms;
	say join ", ", map { "'$_'" } @terms;
	say $pulls;
	my @ret = ($bestI, $term, $minP, $bestCnt, $nWayTie);
	die Dumper(\@ret);
    }
    
    return($term, $minP, $bestCnt->{found}, $nWayTie);
}

# return the index of the minimum element
sub minIndex {
    my( $aref, $ret ) = ( shift, 0 );
    my $min = $aref->[$ret];
    for my $i (1 .. $#{$aref}) {
	if ($aref->[$i] < $min) {
	    $ret = $i;
	    $min = $aref->[$i];
	}
    }
    return $ret;
}
