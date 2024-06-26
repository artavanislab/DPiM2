#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist writeCols);
use DpimLib qw(getLineDP4APMS);

# describe the purpose for this script

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    my %baitPep; # keyed by searchID
    my %topPrey; # keyed by searchID
    my %expt; # expt{$sid} = { date, bait, topPrey }

    open my $IN, "<", $in or die "can't read from $in. $!";
    my %row;
    while (getLineDP4APMS(\%row, $IN)) {
	my $sid = $row{search_id};
	my $pep = $row{total_peptides};
	$expt{$sid} //= {date => $row{sample_date}, bait => 0, topPrey=> 0};
	if ($row{bait_ref} eq $row{prey_ref}) {
	    $expt{$sid}{bait} = $pep
	} elsif ($expt{$sid}{topPrey} < $pep) {
	    $expt{$sid}{topPrey} = $pep;
	}
    }

    my @sid = sort {$a <=> $b} keys %expt;
    my @data = (\@sid);
    for my $col (qw(date bait topPrey)) {
	push @data, [map { $expt{$_}{$col} } @sid];
    }
    
    my $header = join "\t", qw(search_id sample_date bait topPrey);
    my $format = join "\t", qw(%d %s %d %d);
    writeCols($out, \@data, $header, $format);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
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
