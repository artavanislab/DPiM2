#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readHeader readColsRef);
#use DpimLib qw(getLineAPMS);

# replace dummy dates with actual dates in files describing dpim3 2015 results

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $lookupFile = $opts{lookup};
    my $dateCol = $opts{datecol};
    my $sidCol = $opts{sidcol};
    my $instrCol = $opts{inst};
    
    my @cols = ($sidCol, $dateCol);
    my @read;
    readColsRef(\@read, $in, \@cols, 'line');

    my %lookup = readLookup($lookupFile);

    my $header;
    readHeader($in, \$header);
    
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT $header;
    for my $row (@read) {
	my $newDate = $lookup{$row->{$sidCol}}{date};
	
	$row->{line} =~ s/$row->{$dateCol}/$newDate/ if defined $newDate ;
	print $OUT $row->{line};
    }
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	lookup => '/home/glocke/DPiM/dpim4/excelCsv/sid2Date.tab',
	datecol => 'sample_date',
	sidcol => 'search_id',
	instrcol => '',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "lookup=s", "datecol=s", "sidcol=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

sub readLookup {
    my ($in) = @_;

    my @cols = qw(search_id date ms_inst_run_id);
    my @read;
    readColsRef(\@read, $in, \@cols);
    
    return map { $_->{search_id} => { date => $_->{date},
				      ms_inst_run_id => $_->{ms_inst_run_id} }
    } @read;
}
