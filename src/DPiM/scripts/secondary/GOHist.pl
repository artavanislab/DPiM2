#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Storable;

use HomeBrew::IO qw(checkExist writeCols readColsRef);
use DpimLib qw(readGoDB networkHashFromEdgeList);

# make a histogram of the go terms associated with every protein (FBgn) found
# in the input file

my $MAXCNT = 0;
my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};
    my $goDBFile = $opts{godb};
    my $goNameFile = $opts{goname};
    my $proteinCols = $opts{proteincols};
    my $countCol = $opts{countcol};
    $MAXCNT = $opts{dumpafter} if exists $opts{dumpafter};
    
    my %goDB;
    # goDB{proteinID} => [ term1, term2,... ]
    say "parsing $goDBFile";
    if (exists $opts{human}) {
	readGoTab(\%goDB, $goDBFile);
    } else {
	readGoDB(\%goDB, $goDBFile);
    }
    #die Dumper($goDB{3952}, $goDB{4843}, $goDB{10257});

    my ($nodeFlag, $testType);
    if ($mode eq 'edge') {
	$nodeFlag = 0;
	$testType = 'stubs';
    } elsif ($mode eq 'node') {
	$nodeFlag = 1;
	$testType = 'nodes';
    } elsif ($mode eq 'banner') {
	$nodeFlag = -1;
	$testType = 'edges annotated on both ends';
    }

    say "parsing $in";
    my ($protCount, %goHist);

    if ($proteinCols eq 'auto') {
	$protCount = countTerms(\%goHist, $in, \%goDB, $nodeFlag);
    } else {
	my @pCols = split ",", $proteinCols;
	$protCount = countTermsByCol(\%goHist, $in, \%goDB, \@pCols, $countCol);
    }
    
    say "retrieving goNameMap";
    my $goNameMap = retrieve($goNameFile);
    $goNameMap->{UNKNOWN} = { name => 'unknown' };

    {
	# standardize go terms
	my %oldNew;
	for my $go (keys %goHist) {
	    if (exists $goNameMap->{$go}{standard_id} &&
		$goNameMap->{$go}{standard_id} ne $go) {
		$oldNew{$go} = $goNameMap->{$go}{standard_id};
	    }
	}
	for my $oldGo (keys %oldNew) {
	    my $newGo = $oldNew{$oldGo};
	    $goHist{$newGo} += $goHist{$oldGo};
	    delete $goHist{$oldGo};
	}
    }
    
    my @terms = sort keys %goHist;
    #die Dumper(\@terms);
    

    my @goNames = map { '"'.($goNameMap->{$_}{name} // 
			     die "can't find goNames->{$_}{name}").'"'} @terms;

    my @data = (\@terms);
    push @data, [map {$goHist{$_}} @terms];
    push @data, \@goNames;
    
    my $preComments = "# found go terms associated with proteins listed in $in\n";
    $preComments .= "# total tests = $protCount \n";
    $preComments .= "# testing all $testType";
    my $header = join "\t", qw(term count name);
    my $format = join "\t", qw(%s %d %s);
    writeCols($out, \@data, $header, $format, $preComments);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	godb => '/home/glocke/DPiM/flybase/gene_association.fb',
	goname => '/home/glocke/DPiM/flybase/goMap.full.storable',
	proteincols => 'auto',
	countcol => 'auto'
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my @modes = qw(edge node banner);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my $usage = "usage: $0 -in input -out output < $modeString $defaultString ".
	" -dumpafter n -human >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "godb=s", "goname=s", 
	       "proteincols=s", "countcol=s", "dumpafter=i", "human");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{godb});

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

# for every protein in every edge, count GO terms
# if $byNode == 1, do not take degree into account
# if $byNode == -1, only count terms annotated on both stubs
sub countTerms {
    my ($goHist, $in, $goDB, $byNode) = @_;

    my %network;
    networkHashFromEdgeList(\%network, $in, undef, undef, undef, $opts{human});
    if ($byNode == -1) {
	return bannerNet($goHist, \%network, $goDB);
    }
    for my $k1 (keys %network) {
	for my $k2 (keys %{ $network{$k1} }) {
	    $network{$k2}{$k1} = $network{$k1}{$k2};
	}
    }

    my $protCount = 0;
    for my $k1 (keys %network) {
	my $deg = 0+ keys %{ $network{$k1} }; # degree of node
	$deg = 1 if $byNode;
	$protCount += $deg;
	
	my $go = $goDB->{$k1};
	#warn ">\tno GO terms for $k1\n" unless defined $go;
	$go //= ["UNKNOWN"];
	for my $term (@$go) {
	    die "NOT protein = $k1\n", Dumper($go) if $term =~ /NOT/;
	    $goHist->{$term} += $deg;
	}
	die Dumper($goHist) if $MAXCNT && $MAXCNT <= keys %$goHist;
    }
    
    # DEBUG open my $OUT, ">", 'totalEdgeLog';
    #say $OUT $_ for sort keys %test;
    
    return $protCount;
}

# count only terms shared on both ends of the network
sub bannerNet {
    my ($goHist, $network, $goDB) = @_; 

    my $edgeCount = 0;
    for my $k1 (keys %$network) {
	my %go1 = map { $_ => 1 } @{ $goDB->{$k1} };
	#die Dumper($k1, \%go1);
	for my $k2 (keys %{ $network->{$k1} }) {
	    $edgeCount++;
	    for my $term2 (@{ $goDB->{$k2} }) {
		if (exists $go1{$term2}) {
		    $goHist->{$term2}++;
		    #die join "\t", $k1, $k2, $term2;
		}
	    }
	}
    }
    return $edgeCount;
}

# read (at least one) protein name column and a count column
# if count column is "auto", interpret that to mean we're just counting each
#   protein once
sub countTermsByCol {
    my ($goHist, $in, $goDB, $pCols, $countCol) = @_;

    my $autoFlag = ($countCol eq 'auto');

    my @cols = @$pCols;
    push @cols, $countCol unless $autoFlag;

    my $protCount;

    my @read;
    readColsRef(\@read, $in, \@cols);
    for my $row (@read) {
	my $count = ($autoFlag)?1:$row->{$countCol};
	$protCount+=$count;
	for my $col (@$pCols) { 
	    my $go = $goDB->{$row->{$col}}; # all GO terms for this protein
	    warn ">\tno GO terms for $row->{$col}\n" unless defined $go;
	    $go //= ["UNKNOWN"];
	    for my $term (@$go) {
		die "NOT protein = $row->{$col}\n"
		    , Dumper($go) if $term =~ /NOT/;
		$goHist->{$term} += $count;
	    }
	}
	#die dumpHash($goHist) if $protCount > $count;
    }
    

    return $protCount;
}

sub dumpHash {
    my $h = shift;

    for my $k (sort keys %$h) {
	print "$k => ", Dumper($h->{$k});
    }
}
