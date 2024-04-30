#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(max);
use IO::Zlib;
use HomeBrew::IO qw(checkExist readList);
#use DpimLib qw(getLineAPMS);

# translate uniprot id's to entrez id's
# use ncbi's gene2accession flat file

my %opts = getCommandLineOptions();

{
    my $uniprotFile = $opts{in};
    my $out = $opts{out};
    my $dbFile = $opts{db};

    my %ids = map {$_ => 0} readList($uniprotFile);

    my $cnt = 0;
    my $matches = 0;
    
    my $entrezCol = 1;
    my $uniCol = 5;
    my $symbolCol = -1;
    #my @cols = ($uniCol, $entrezCol, $symbolCol);
    
    my $DB = new IO::Zlib;
    die "can't read from $dbFile" unless $DB->open($dbFile, "rb");
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT join "\t", qw(uniprot entrez symbol);
    while (my $line = <$DB>) {
	next if $line =~ /^#/;
	say "$cnt ($matches matches)" if 0 == ($cnt % 100000);
	$cnt++;

	my @spl = split /\t/, $line;
	my $uni = $spl[$uniCol];
	$uni =~ s/\.\d+$//;
	next unless exists $ids{$uni};

	$matches++;
	$ids{$spl[$uniCol]}++;

	print $OUT join "\t", $uni, $spl[$entrezCol], $spl[$symbolCol];
	# this relies on the fact that symbolCol includes a linebreak
    }

    my %table;
    for my $i (values %ids) {
	$table{$i}++;
    }
    my $max = max( keys %table );
    for my $i (0..$max) {
	say "count($i) = ", $table{$i}//0;
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	db => '/home/glocke/DPiM/corum/gene2accession.gz',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in uniprot.list -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "db=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
