#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min);
use List::MoreUtils qw(each_array);
use HomeBrew::IO qw(checkExist readColsRef);
#use DpimLib qw(getLineAPMS);

# many entrez id's in the corum csv file are incorrect
# some are given as 0 when a good match exists
# some are close matches (Calm1 incorrectly identified as Calm3)
# some are totally wrong (mouse Cox1 incorrectly identified as Ptgs1)

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $uniMapFile = $opts{unimap} // die "ehh???";

    my %uniMap;
    readUniMap(\%uniMap, $uniMapFile);

    my $idCol = 0;
    my $nameCol = 1;
    my $speciesCol = 3;
    my $uniCol = 4;
    my $entrezCol = 5;

    open my $IN, "<", $in or die "Can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't read from $in. $!";
    my $line = <$IN>;
    print $OUT $line; # reproduce header
    while ($line = <$IN>) {
	my @spl = split /;/, $line;

	my @uni = split /,/, $spl[$uniCol];
	@uni = nested(\@uni);
	my @oldEntrez = split /,/, $spl[$entrezCol];
	@oldEntrez = nested(\@oldEntrez);
	my @translate = translate(\@uni, \%uniMap);
	my @newEntrez = replace(\@oldEntrez, \@translate);
	#die Dumper(\@newEntrez) if $spl[0] == 26;
	my $eString = flatten(\@newEntrez);
	warn Dumper(\@oldEntrez, \@newEntrez, $line) 
	    if $eString ne $spl[$entrezCol];
	$spl[$entrezCol] = $eString;
	print $OUT join ";", @spl;
    }
    
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	unimap => '/home/glocke/DPiM/corum/ncbi_corumUniProtID_to_entrez.tab',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "unimap=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{unimap});

    return %opts;
}

sub readUniMap {
    my ($ret, $in) = @_;

    #my @cols = qw(From To);
    my @cols = qw(uniprot entrez);
    my @read;
    readColsRef(\@read, $in, \@cols);
    
    my $duplicates = 0;
    for my $row (@read) {
	if (exists $ret->{$row->{uniprot}}) {
	    $duplicates++;
	    if (ref($ret->{$row->{uniprot}})) {
		push @{$ret->{$row->{uniprot}}}, $row->{entrez};
	    } else {
		$ret->{$row->{uniprot}} = 
		    [ $ret->{$row->{uniprot}}, $row->{entrez} ];
	    }
	} else {
	    $ret->{$row->{uniprot}} = $row->{entrez};
	}
    }
    warn "$duplicates uniprot ids out of ", (0+ keys %$ret)
	, " map to multiple Entrez ID's.  Picking the lesser value in each case"
	if $duplicates > 0;

    return;
}

# list ambiguous complex members in a nested array ref aside non-ambiguous 
#   proteins
#
# corum lists ambiguous complex members, like this: id1, id2, (ambig1, ambig2)
# output format: [ 'id1', 'id2', ['ambig1', 'ambig2'] ]
sub nested {
    my ($ids) = @_;

    my $ambigFlag = undef;
    my $tmpList;
    my @complex;
    for my $id (@$ids) {
	if ($id =~ /\(/) {
	    $id =~ s/\(//g; 

	    # open an ambiguous member
	    $tmpList = [];
	    $ambigFlag = 1;
	}
	if ($ambigFlag) {
	    my $tmpID = $id;
	    $tmpID =~ s/\)//g;
	    push @$tmpList, $tmpID;
	} else {
	    push @complex, $id;
	}
	if ($id =~ /\)/) {
	    # close an ambiguous member
	    die "double nesting" if ! $ambigFlag;
	    push @complex, $tmpList;
	    $tmpList = undef;
	    $ambigFlag = undef;
	} 
    }

    return @complex;
}

# translate uniprot to entrez 
# respect the nested structure of ambiguous clusters
sub translate {
    my ($uni, $uniMap, $depth) = @_;
    $depth //=0;

    my @translated;
    for my $u (@$uni) {
	if (ref($u)) {
	    push @translated, [translate($u, $uniMap, $depth+1)];
	} else {
	    push @translated, $uniMap->{$u} // 0;
	}
    }
    
    return @translated;
}

# replace old entrez id's with new ones if the old one is zero
# respect the nested structure of ambiguous clusters
sub replace {
    my ($old, $new, $uni, $changes) = @_;
    die "replace: bad format. Dumping", Dumper($old, $new) unless ref($old);

    my @ret;
    my $it = each_array(@$old, @$new, @$uni);	
    while (my ($o, $n, $u) = $it->()) {
	my $add;
	if (ref($o)) {
	    die "ref(o) but not ref(n). \n", Dumper($old, $new) if ! ref($n);
	    $add = [replace($o, $n, $u)];
	} else {
	    $add = ($o == 0)?$n:$o;
	    if ($o != 0 && $n != 0 && $o != $n) {
		warn "!!!!mismatch.  o = $o, n = $n\n", Dumper($old, $new);
	    }
	}
	push @ret, $add;
    }
    
    return @ret;
}

# basically like join "," except it must
# respect the nested structure of ambiguous clusters
sub flatten {
    my ($arr) = @_;
    
    for my $a (@$arr) {
	if (ref($a)) {
	    $a = join ",", @$a;
	    $a = "($a)";
	}
    }
    
    return join ",", @$arr;
}
