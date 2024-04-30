#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsRef readColsHashRef);
#use DpimLib qw(readHS);

# select all edges connecting annotated edges

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $term = $opts{term};
    my $goFile = $opts{go};
    my $out = $opts{out};
    my $transFile = $opts{trans};
    
    my %go;
    readGO(\%go, $goFile, $term);
    die "can't find '$term' in $goFile\n" if 0 == keys %go;
    
    
    open my $IN, "<", $in or die "can't read $in. $!";
    my $OUT = *STDOUT;
    if (defined $out) {
	open $OUT, ">", $out or die "can't write to $out. $!";
    }
    my %found;
    while (<$IN>) {
	my ($p1, $p2) = split;
	$found{$p1} = 1 if exists $go{$p1};
	$found{$p2} = 1 if exists $go{$p2};
	next;
	next unless exists $go{$p1} || exists $go{$p2};
	print $OUT $_;
    }
    close $IN; 
    close $OUT if defined $out;

    my @nodes = keys %found;
    if (exists $opts{tosymb}) {
	@nodes = toSymb(\@nodes, $transFile);
    } else {
	@nodes = sort { $a <=> $b } @nodes;
    }
    say for @nodes;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	go => '/home/glocke/DPiM/human/GO/entrez.go.tsv',
	trans => '/home/glocke/DPiM/human/nsaf/entrez2symbol.concise.tsv',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in net -term GO:000000 < -out output -tosymb ".
	" $defaultString > \n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "term=s", "out=s", "tosymb", "go=s", 
	       "trans=s");
    die $usage unless exists $opts{in} && exists $opts{term};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{go});

    return %opts;
}

sub readGO {
    my ($ret, $in, $term) = @_;

    my @read;
    readColsRef(\@read, $in, [qw(entrez term)]);
    for my $row (@read) {
	$ret->{$row->{entrez}} = 1 if $row->{term} eq $term;
    }

    return;
}

sub toSymb {
    my ($nodes, $transFile) = @_;
 
    my %trans;
    readColsHashRef(\%trans, $transFile, [qw(entrez symbol)]);
    return sort map { ($trans{$_} // die "can't convert $_ to symbol") } @$nodes;
}
