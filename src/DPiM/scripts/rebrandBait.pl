#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineDP4APMS getLineBiopAPMS);

# for each experiment, find the prey with the highest TSC, then reset the bait
#   to be that prey
# pick the lower fbgn in the case of a tie

my %opts = getCommandLineOptions();

{
    my $noBaitFile = $opts{nobait};
    my $out = $opts{out};
    my $withBaitFile = $opts{withbait};

    my %apms;
    # apms{search_id} = [ {all columns},...]
    
    my @cols = qw( search_id bait_ref prey_ref total_peptides sample_date ms_inst_run_id );
    push @cols, qw(rep plate) if $opts{mode} eq 'human';

    open my $IN, "<", $noBaitFile or die "Can't read from $noBaitFile. $!";
    my %row;
    say "parse";
    my $reader = \&getLineDP4APMS;
    $reader = \&getLineBiopAPMS if $opts{mode} eq 'human';
    while ($reader->(\%row, $IN)) {
	my $id = $row{search_id};
	$apms{$id} //= [];

	push @{ $apms{$id} }, { map { $_ => $row{$_} } @cols};
    }
    close $IN;

    for my $expt (values %apms) {
	#say "first bait = $expt->[0]{bait_ref}";
	my $maxPrey;
	my $maxTSC = 0;
	for my $row (@$expt) {
	    if ($maxTSC && $row->{total_peptides} == $maxTSC) {
		$maxPrey = (sort ($row->{prey_ref}, $maxPrey))[0];
	    } elsif ($row->{total_peptides} > $maxTSC) {
		$maxTSC = $row->{total_peptides};
		$maxPrey = $row->{prey_ref};
	    }
	}
	$_->{bait_ref} = $maxPrey for @$expt;
	#say "last bait = $expt->[0]{bait_ref}, maxTSC = $maxTSC; $maxPrey";
	#exit;
    }

    my $OUT;
    if (defined $withBaitFile) {
	open $OUT, ">", $out or die "can't write to $out. $!";
	say $OUT "# rebranded baits in $noBaitFile";
	close $OUT;
	my $cmd = "cat $withBaitFile >> $out";
	system($cmd);
	open $OUT, ">>", $out or die "can't write to $out. $!";
    } else {
	open $OUT, ">", $out or die "can't write to $out. $!";
	say $OUT "# rebranded baits in $noBaitFile";
	say $OUT join "\t", @cols;
    }

    for my $sid (sort {$a <=> $b} keys %apms) {
	for my $row (@{ $apms{$sid} }) {
	    say $OUT join "\t", map { $row->{$_} } @cols;
	}
    }
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    my @modes = qw(fly human);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -nobait baitless.apms -out output < $modeString ".
	"-withbait apms > \n";

    my %opts = ();
    GetOptions(\%opts, "nobait=s", "out=s", "mode=s", "withbait=s");
    die $usage unless exists $opts{nobait} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{mode} //= $modes[0];
    die "you must select one of the following modes: "
	, join(", ",keys(%modes)), "\n" if ! exists $modes{$opts{mode}};

    checkExist('f', $opts{nobait});
    checkExist('f', $opts{withbait}) if exists $opts{withbait};

    return %opts;
}
