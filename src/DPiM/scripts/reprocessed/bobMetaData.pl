#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsHashRef readColsRef);
use DpimLib qw(getLineMayAPMS);

# make search_id -> bait+retain information

##my %opts = getCommandLineOptions();
{
    my $idFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/idString2SearchID.tsv';
    my $bobFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_102216_noComments.txt';
    my $outFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.10-31-2016.tsv';
    

    my %id2id;
    readColsHashRef(\%id2id, $idFile, [qw(id_string search_id)]);

    my @bobAnn;
    my @inCols = qw(id_string Final_Updated_Bait_ID RAO_Determination);
    readColsRef(\@bobAnn, $bobFile, \@inCols, undef, "\t");
    
    open my $OUT, ">", $outFile or die "can't write to $outFile: $!";
    say $OUT "# $0 connect search_id to its bait and retain information";
    say $OUT join "\t", qw(search_id bait retain);
    for my $row (@bobAnn) {
	my @row = map { $row->{$_} } @inCols;
	$row[0] = $id2id{$row[0]} // die "Can't map $row[0] to search_id !!\n".
	    Dumper($row);;
	say $OUT join "\t", @row;
    }
    exit;
}


{
    ## map id_string to search_id
    
    my @files = readList('/home/glocke/DPiM/augRemap/apmsData/raw.list');
    my %ids;
    my %row;
    for my $f (@files) {
	open my $IN, "<", $f or die "can't read from $f: $!";
	while (getLineMayAPMS(\%row, $IN)) {
	    $ids{$row{id_string}}{$row{search_id}} = 1;
	}
    }

    my @multiKeys = grep { 1 < keys %{ $ids{$_} }} keys %ids;
    say "# $0 compiled a map from id string to search_id";
    say "# $_" for @multiKeys;
    say "# number of multi keys = ", 0+ @multiKeys;
    say join "\t", qw(search_id id_string);
    for my $id (sort keys %ids) {
	my @sid = keys %{ $ids{$id} };
	say join "\t", $_, $id for @sid;
    }
    
    exit;
}
my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    
    open my $IN, "<", $in or die "can't read $in. $!";
    close $IN; 
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 describe this";
    close $OUT; 
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	in => 
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

    checkExist('f', $opts{in});

    return %opts;
}
