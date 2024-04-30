#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsRef writeCols);
#use DpimLib qw(getLineAPMS);

# turn corum clusters listed in allComplexes.csv into something I can use
# here, we are taking the map compiled by parseCorum0 and applying it to the 
#   corum clusters
# that is, we are translating clusters listed by human entrez-ids into clusters
#   listed by FBgn

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $orthoFile = $opts{ortho};

    my @clusters;
    parseCorum(\@clusters, $in);
    
    my %orthos;
    parseOrthos(\%orthos, $orthoFile);

    for my $c (@clusters) {
	$c->{cluster} = replace($c->{cluster}, \%orthos);
	$c->{string}  = flatten($c->{cluster});
	unless(exists $opts{nested}) {
	    $c->{string} =~ s/\(//g;
	    $c->{string} =~ s/\)//g;
	}
    }
    @clusters = grep { $_->{string} } @clusters; # remove clusters with no members
    
    my @data = ([map { $_->{id} }     @clusters],
		[map { $_->{name} }   @clusters],
		[map { $_->{string} } @clusters]);
    my $header = join "\t", qw(id name proteins);
    my $format = join "\t", qw(%d %s %s);
    my $preComments = "# updated $in according to $orthoFile";
    writeCols($out, \@data, $header, $format, $preComments);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	ortho => '/home/glocke/DPiM/corum/topOrthologs.tab',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in corum.clusters -out fbgn.clusters < ".
	"$defaultString -nested >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "ortho=s", "nested");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{ortho});

    return %opts;
}

# ret[$i] = { id => $corumID, name => "name of cluster", cluster => [nested] }
sub parseCorum {
    my ($ret, $inFile) = @_;

    open my $IN, "<", $inFile or die "Can't read from $inFile. $!";
    
    my %allEntrez;

    my $idCol = 0;
    my $nameCol = 1;
    my $speciesCol = 3;
    my $entrezCol = 5;
    <$IN>; # ignore header
    while (<$IN>) {
	my @spl = split /;/;
	##next unless $spl[$speciesCol] =~ /uman/ || 
	##    $spl[$speciesCol] =~ /Mammalia/;
	#next if $spl[$entrezCol] =~ /\(/;
	my @entrez = split /,/, $spl[$entrezCol];
	
	my @cluster;
	my $ambigFlag = undef;
	# this loop is complicated because corum lists ambiguous cluster 
	# members, e.g. id1, id2, (ambig1, ambig2)
	my $tmpList;
	for my $entr (@entrez) {
	    if ($entr =~ /\(/) {
		# open an ambiguous member
		$ambigFlag = 1;
		$entr =~ s/\(//g; 
	    }
	    if ($ambigFlag) {
		my $tmpID = $entr;
		$tmpID =~ s/\)//g;
		$tmpList = [];
		push @$tmpList, $tmpID;
	    } else {
		push @cluster, $entr;
	    }
	    if ($entr =~ /\)/) {
		# close an ambiguous member
		die "double nesting" if ! $ambigFlag;
		push @cluster, $tmpList;
		$tmpList = undef;
		$ambigFlag = undef;
	    } 
	}
	push @$ret, { id => $spl[$idCol], name => $spl[$nameCol], 
		      cluster => \@cluster};
    }

    return;
}

sub parseOrthos {
    my ($ret, $orthoFile) = @_;

    my @cols = qw(entrez fbgn);
    my @read;
    readColsRef(\@read, $orthoFile, \@cols);

    for my $row (@read) {
	$ret->{$row->{entrez}}{$row->{fbgn}} = 1;
    }

    return;
}

# when there are multiple orthologs, just add all of them
# when there is no ortholog, add nothing
sub replace {
    my ($complex, $map) = @_;

    #my $null = 'FBgn0000000';

    my @ret = ();
    for my $member (@$complex) {
	if (ref($member)) {
	    push @ret, replace($member, $map); 
	} else {
	    #if (! exists $map->{$member}) {
	    #	push @ret, $null;
	    #   }
	    push @ret, $_ for keys %{ $map->{$member} };
	}
    }
    #die Dumper(\@ret) if $complex->[0] == 14682;
    
    return \@ret;
}

# basically like join "," except it must
# respect the nested structure of ambiguous clusters
sub flatten {
    my ($arr) = @_;
    
    my @flat;
    for my $a (@$arr) {
	if (ref($a)) {
	    push @flat, @$a;
	} else {
	    push @flat, $a;
	}
    }
    my %rmDupes = map { $_ => 1 } @flat;

    return join ",", sort keys %rmDupes;
}
