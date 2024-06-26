#!/usr/bin/env perl

use strict;
use warnings;
use v5.10; 
use Data::Dumper;
use DateTime;
use DateTime::Format::Strptime;

# parse csv files to extract a map from seach_id to date and instr_id

if (@ARGV != 1) {
    die "usage: $0 outFile\n";
}

{
    # previously located at /home/glocke/DPiM/dpim4/excelCsv/sid2DateInstrID.tsv
    my $out = $ARGV[0];

    my @csv = qw'
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_011615.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_040215.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_042015.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_050815.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_052215.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_060815.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_061815.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_071015.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_072415.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_073115.csv
'; #

    my $parseDate = DateTime::Format::Strptime->new(
        pattern   => '%m%d%y');

    my %dates; # dates{$dateTime}{$search_id}=1
    my %dateTime;
    for my $f (@csv) {
	$f =~ /_(\d+)\.csv/ or die "can't parse $f.";
	my $date = $parseDate->parse_datetime($1);
	#say "$f $date";
	$dateTime{$date} = $date;

	open my $IN, "<", $f or die "Can't read from $f. $!";
	my $header = <$IN>; # strip header
	my @spl = split /,/, $header;
	my $sidCol = (grep { $spl[$_] =~ /search_id/i } 0..$#spl)[0];
	#my $instrCol = (grep { $spl[$_] =~ /name/i } 0..$#spl)[0];
	my $instrCol = 0;
	die "can't find search_id in $f" unless defined $sidCol;
	$sidCol += 2 if $f =~ /072415.csv/; ## grumble grumble
	## these excel files are unique snowflakes!

	while (<$IN>) {
	    my @spl = split /,/;
	    my $sid = $spl[$sidCol] or die "f = $f, spl = ", join "\t", @spl;
	    my $instr = dropUnderscore($spl[$instrCol]);
	    $dates{$date}{$sid}=$instr;
	}
    }
    
    my $writeDate = DateTime::Format::Strptime->new(
        pattern   => '%Y-%m-%d');
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    
    say $OUT "search_id\tdate\tms_inst_run_id";
    for my $d ( sort {$dateTime{$a} <=> $dateTime{$b}} keys %dates) {
	my $date = $dateTime{$d};
	for my $sid (sort {$a <=> $b} keys %{$dates{$d}} ) {
	    say $OUT join "\t", $sid, $writeDate->format_datetime($date)
		, $dates{$d}{$sid};
	}
    }
    

    exit;
}

sub dropUnderscore {
    my $inst = shift;
    my @spl = split /\_/, $inst;
    return $spl[0];
}


if (0) {
    # parse csv files to extract a date/seach_id map

    my @csv = qw(
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_011615.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_040215.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_042015.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_050815.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_052215.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_060815.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_061815.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_071015.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_072415.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_073115.csv
);

    my $parseDate = DateTime::Format::Strptime->new(
        pattern   => '%m%d%y');

    my %dates; # dates{$dateTime}{$search_id}=1
    my %dateTime;
    for my $f (@csv) {
	$f =~ /_(\d+)\.csv/ or die "can't parse $f.";
	my $date = $parseDate->parse_datetime($1);
	#say "$f $date";
	$dateTime{$date} = $date;

	open my $IN, "<", $f or die "Can't read from $f. $!";
	my $header = <$IN>; # strip header
	my @spl = split /,/, $header;
	my $sidCol = (grep { $spl[$_] =~ /search_id/i } 0..$#spl)[0];
	die "can't find search_id in $f" unless defined $sidCol;
	$sidCol += 2 if $f =~ /072415.csv/; ## grumble grumble
	## these excel files are unique snowflakes!

	while (<$IN>) {
	    my @spl = split /,/;
	    my $sid = $spl[$sidCol] or die "f = $f, spl = ", join "\t", @spl;
	    $dates{$date}{$sid}=1;
	}
    }
    
    my $writeDate = DateTime::Format::Strptime->new(
        pattern   => '%Y-%m-%d');
    
    say "search_id\tdate";
    for my $d ( sort {$dateTime{$a} <=> $dateTime{$b}} keys %dates) {
	my $date = $dateTime{$d};
	for my $sid (sort {$a <=> $b} keys %{$dates{$d}} ) {
	    say join "\t", $sid, $writeDate->format_datetime($date);
	}
    }
    

    exit;
}
