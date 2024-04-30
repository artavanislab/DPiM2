#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readCol);
use Storable;
use DpimLib qw(readGoDB);

# given a network/list of FBgns, and a list of GO terms,
# find and report all proteins bearing those GO terms

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $protFile = $opts{prot};
    my $goCol = $opts{gocol};
    my $goDBFile = $opts{godb};
    my $goNameFile = $opts{goname};

    my @terms = readCol($in, $goCol);
    
    my %goDB;
    say "parsing $goDBFile";
    readGoDB(\%goDB, $goDBFile);

    my %GOProtList = map { $_ => {} } @terms; 
    # GOProtList{term} = { protein1 => cnt, protein2 => cnt, ...}
    if ($protFile ne 'all') {
	my %prot = readProt($protFile);

	for my $p (keys %prot) {
	    for my $term (@{ $goDB{$p} }) {
		next unless exists $GOProtList{$term};
		$GOProtList{$term}{$p}++;
	    }
	}
    } else {
	die "listing all proteins not currently implemented";
    }

    # maps GO:XXXX to "name" of term
    say "retrieving goNameMap";
    my $goNameMap = retrieve($goNameFile);
    $goNameMap->{UNKNOWN} = { name => '"unknown"' };
    my $preComments = "# found go terms in $in.  found proteins in $protFile";
    writer($out, \@terms, \%GOProtList, $goNameMap, $preComments);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	gocol => 'term',
	godb => '/home/glocke/DPiM/flybase/gene_association.fb',
	goname => '/home/glocke/DPiM/flybase/goMap.storable',
	prot => 'all',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in go.tab -out output < -prot file/'all' ".
	"$defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "prot=s", "gocol=s", "goname=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{prot}) unless $opts{prot} eq 'all';
    checkExist('f', $opts{godb});
    checkExist('f', $opts{goname});

    return %opts;
}

sub readProt {
    my ($in) = @_;

    open my $IN, "<", $in or die "Can't read from $in. $!";

    my %ret;
    while (my $line = <$IN>) {
	while ($line =~ /(FBgn\d+)/g) {
	    $ret{$1}++;
	}
    }

    close $IN;

    return %ret;
}

sub writer {
    my ($out, $terms, $GOProtList, $goNameMap, $preComments) = @_;

    open my $OUT, ">", $out or die "Can't write to $out. $!";

    say $OUT $preComments if defined $preComments;
    
    say $OUT join "\t", qw(term name proteins);

    for my $go (@$terms) {
	my @proteins = sort { $GOProtList->{$go}{$b} <=> $GOProtList->{$go}{$a} }
	    keys %{ $GOProtList->{$go}};
	say $OUT join "\t", $go, '"'.$goNameMap->{$go}{name}.'"', @proteins;
    }
    close $OUT;
}