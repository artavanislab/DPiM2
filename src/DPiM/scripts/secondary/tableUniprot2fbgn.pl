#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(max min);
use HomeBrew::IO qw(checkExist readList readCol readColsRef readColsHashRef);
use DpimLib qw(getLineBiopAPMS networkHashFromEdgeList);

# make a table mapping from uniprot id's to FBgn (many to one)

{
    my $uniTransFile = '/home/glocke/DPiM/uniprot/uniprotOrgConversion.tsv';
    my %uniTrans;
    readUniTable(\%uniTrans, $uniTransFile);

    my $bioTransFile = '/home/glocke/DPiM/uniprot/biomartTrans.tsv';
    my %bioTrans;
    readBioTable(\%bioTrans, $bioTransFile);

    my %byHand;
    my $byHandFile = '/home/glocke/DPiM/uniprot/byHand.tsv';
    readColsHashRef(\%byHand, $byHandFile, [qw(uniprot fbgn)]);
    
    my %trans = %uniTrans;
    for my $uni (keys %bioTrans) {
	$trans{$uni}{$_} = 1 for keys %{ $bioTrans{$uni} };
    }
    $trans{$_}{$byHand{$_}}=1 for keys %byHand;
    #die Dumper(\%trans);

    my $out = '/home/glocke/DPiM/uniprot/uniprot2fbgn.many2many.tsv';
    open my $OUT, ">", $out or die "Can't write to $out: $!";
    say $OUT "# $0 merged $uniTransFile, $bioTransFile and $byHandFile";
    say $OUT join "\t", qw(uniprot fbgn);
    for my $uni (sort keys %trans) {
	for my $fbgn (sort keys %{ $trans{$uni} }) {
	    say $OUT join "\t", $uni, $fbgn;
	}
    }
    
    exit;
}

sub readUniTable {
    my ($ret, $inFile) = @_;

    open my $IN, "<", $inFile or die "can't read from $inFile: $!";

    while (<$IN>) {
	next if /^yourlist/;
	my @spl = split;
	my $uni = $spl[1];
	$uni = $spl[2] if $uni eq ',';
	$spl[0] =~ s/FBGN/FBgn/g;
	$ret->{$uni}{$_} = 1 for split /,/, $spl[0] ;
    }
    return;
}

sub readBioTable {
    my ($ret, $inFile) = @_;

    open my $IN, "<", $inFile or die "can't read from $inFile: $!";

    while (<$IN>) {
	next if /^flybase_gene_id/;
	my @spl = split;
	next if @spl < 2;
	my $uni = $spl[1];
	die $_ if 3 > length($uni);
	$ret->{$uni}{$spl[0]}=1;
    }
    return;
}
