#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use DateTime;
use DateTime::Format::Strptime;
use HomeBrew::IO qw(checkExist writeCols);
use DpimLib qw(getLineRawAPMS getLineAPMS);

# count up the # of peptides and # of unique proteins obtained per day
#   and # of batches

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};

    my $reader = \&getLineRawAPMS;
    $reader = \&getLineAPMS if $mode eq 'clean';

    my $parseDate = DateTime::Format::Strptime->new(
        pattern   => '%Y-%m-%d');

    my %days; # days{$date} = { tsc => $sum, prot => {$uniqFbns=> 1}, 
    #                           sid => {$uniqSearchIDs=>1}, date=>$dateTime }
    
    open my $IN, "<", $in or die "can't open $in. $!";
    my %row;
    while ($reader->(\%row, $IN)) {
	my $date = $row{sample_date};
	$days{$date}{date} //= $parseDate->parse_datetime($date);
	$days{$date}{tsc}+=$row{total_peptides};
	$days{$date}{prot}{$row{prey_ref}}=1;
	$days{$date}{sid}{$row{search_id}}=1;
    }
    close $IN;

    my @days = sort { $days{$a}{date} <=> $days{$b}{date} } keys %days;
    my @tsc  = map { $days{$_}{tsc} } @days;
    my @prot = map { 0+ keys %{ $days{$_}{prot} } } @days;
    my @sid  = map { 0+ keys %{ $days{$_}{sid} } } @days;

    my @cols = qw(date tsc prot expts);
    my $header = join "\t", @cols;
    my $format = join "\t", qw(%s %d %d %d);
    my $preComments = "# obtained daily statistics from $in";
    writeCols($out, [\@days, \@tsc, \@prot, \@sid], $header, $format, 
	      $preComments)
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my @modes = qw(raw clean);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $modeString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    checkExist('f', $opts{in});

    return %opts;
}
