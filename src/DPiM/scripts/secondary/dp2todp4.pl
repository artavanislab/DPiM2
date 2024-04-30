#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(max);
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineDP4APMS);

# describe the purpose for this script

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 reformatted $in";
    parseDP2($in, $OUT);
    close $OUT; 
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


# remove FBgn0000000 
sub parseDP2 {
    my ($in, $OUT) = @_;
    
    my @cols = qw(tap_id  ms_inst_run_id  user    search_id       sample_date     total_peptides unique_peptides  bait_ref        prey_ref);
    my ($sid) =  grep { $cols[$_] =~ /search_id/ } 0..$#cols;
    my ($bait) = grep { $cols[$_] =~ /bait_ref/ } 0..$#cols;
    my ($prey) = grep { $cols[$_] =~ /prey_ref/ } 0..$#cols;
    my ($tsc) =  grep { $cols[$_] =~ /total_peptides/ } 0..$#cols;
    my ($date) =  grep { $cols[$_] =~ /sample_date/ } 0..$#cols;
    my ($instID) = grep { $cols[$_] =~ /ms_inst_run_id/ } 0..$#cols;
    my $parseLine = makeParser($sid, $bait, $prey, $tsc, $date, $instID);

    open my $IN, "<", $in or die "Can't read from $in. $!";
    
    # skip header
    while (<$IN>) {
	last if /tap_id/;
    }

    while (<$IN>) {
	next if /FBgn0000000/;
	say $OUT join "\t", $parseLine->($_);
    }
}


# return a closure that returns (search_id bait_ref prey_ref)
sub makeParser {
    my ($sid, $bait, $prey, $tsc, $date, $instID) = @_;

    my $maxI = max($sid, $bait, $prey, $tsc, $date, $instID);

    return sub {
	my $line = shift;
	chomp;
	my @spl = split /\t/;
	die "too few tokens in '$line'" unless $maxI < @spl;
	return ($spl[$sid], $spl[$bait], $spl[$prey], $spl[$tsc], $spl[$date],
		$spl[$instID]);
    }
}
