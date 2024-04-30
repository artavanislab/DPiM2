#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min sum);
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineRawAPMS getLineAPMS);

# report number of
#  * experiments
#  * experiments with 0, 1, 2+ bait peptides
#  * prey
#  * baits
#  * baits with 1, 2, 3, 4+ experiments

my %opts = getCommandLineOptions();

{
    my $out = $opts{out};
    my $dp2_clean = $opts{dp2_clean};
    my $dp2_raw = $opts{dp2_raw};
    my $dp3 = $opts{dp3};

    #    my (%dp2_clean, %dp2_raw, %dp3, %dp23_clean, %dp23_raw);
    my %rawStats; # rawStats{$experiment}{protein}{$fbgn} = 
    #             #   { bait=>[$search_id, ], bait_peptides => [0, 21, 33], 
    #             #     prey=>[$search_id, ] (not including bait-bait counts)
    #             # rawStats{$experiment}{id} = {$search_id => 1}
    say "collecting dp2_clean"; 
    collectStats(\%rawStats, $dp2_clean, 'dp2_clean');
    say "collecting dp2_raw";
    collectStats(\%rawStats, $dp2_raw, 'dp2_raw');
    say "collecting dp3";
    collectStats(\%rawStats, $dp3, 'dp3');

    my @stats = qw(experiments experimentsBait0 experimentsBait1 
                   experimentsBait2Plus prey bait bait1 bait2 bait3Plus);
    my %summaryStats; # summaryStats{$exper} = map {$_ => stat} @stats;
    #                 # likewise summary{dp2_clean_dp3} when combining them
    say "summarizing dp2_clean"; 
    $summaryStats{'dp2_clean'} = summarize($rawStats{'dp2_clean'});
    say "summarizing dp2_raw"; 
    $summaryStats{'dp2_raw'} = summarize($rawStats{'dp2_raw'});
    say "summarizing dp3"; 
    $summaryStats{'dp3'} = summarize($rawStats{'dp3'});
    say "summarizing dp2_clean_dp3"; 
    $summaryStats{'dp2_clean_dp3'} = summarize($rawStats{'dp2_clean'}, 
					       $rawStats{'dp3'});
    say "summarizing dp2_raw_dp3"; 
    $summaryStats{'dp2_raw_dp3'} = summarize($rawStats{'dp2_raw'}, 
					     $rawStats{'dp3'});
    #die Dumper(\%summaryStats);
    
    my @columns = qw(dp2_clean dp2_raw dp3 dp2_clean_dp3 dp2_raw_dp3);
    writeMatrix($out, \%summaryStats, \@stats, \@columns);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	dp2_clean => '/home/glocke/DPiM/dpim3.1/dpim2_nrtap.120123',
	dp2_raw => '/home/glocke/DPiM/dpim3_raw/dpim2_all.120123.plusUniqMsInst',
	dp3 => '/home/glocke/DPiM/dpim3.1/DPiM_Data_Summary_2014_edited.reformatted.rmDupes.tab',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "dp2_clean=s", "dp2_raw=s", "dp3=s");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{dp2_clean});
    checkExist('f', $opts{dp2_raw});
    checkExist('f', $opts{dp3});

    return %opts;
}

# output: 
#
# stats{$experiment}{protein}{$fbgn} = 
#   { bait=>[$search_id, ], bait_peptides => [0, 21, 33], 
#     prey=>[$search_id, ] (not including bait-bait counts)
#     stats{id}{$experiment} = {$search_id => 1}sub collectStats {
sub collectStats {
    my ($stats, $file, $exper) = @_;

    my $reader = \&getLineRawAPMS;
    $reader = \&getLineAPMS if $exper eq 'dp2_clean';

    my %byID; # byID{id} = {bait=>$fbgn, prey { $fbgn => $TSC }}
    
    open my $IN, "<", $file or die "Can't read from $file. $!";
    my %row;
    while ($reader->(\%row, $IN)) {
	my $id = $row{search_id};
	$byID{$id}{bait} = $row{bait_ref};
	$byID{$id}{prey}{$row{prey_ref}} += $row{total_peptides};
    }

    for my $id (keys %byID) {
	$stats->{$exper}{id}{$id} = 1;
	my $bait = $byID{$id}{bait};
	$stats->{$exper}{protein}{$bait}{bait} //= [];
	push @{ $stats->{$exper}{protein}{$bait}{bait} }, $id;
	
	my $bait_peptides = 0;
	for my $prey (keys %{ $byID{$id}{prey} }) {
	    if ($prey eq $bait) {
		$bait_peptides+=$byID{$id}{prey}{$prey};
	    } else {
		$stats->{$exper}{protein}{$prey}{prey} //= [];
		push @{ $stats->{$exper}{protein}{$prey}{prey} }, $id;
	    }
	}
	
	$stats->{$exper}{protein}{$bait}{bait_peptides} //= [];
	push @{ $stats->{$exper}{protein}{$bait}{bait_peptides} }, 
	    $bait_peptides;
    }

    return;
}

sub summarize {
    my @experiments = @_;

    my @countExperiments = countExperiments(@experiments);
    my $prey = countPrey(@experiments);
    my @countBaits = countBaits(@experiments);
    
    my @statNames = qw(experiments experimentsBait0 experimentsBait1 
                       experimentsBait2Plus prey bait bait1 bait2 bait3Plus);
    my @stats = (@countExperiments, $prey, @countBaits);
    my %ret = map { $statNames[$_] => $stats[$_] } 0..$#statNames;

    return \%ret;
}

# experiment->{protein}{$fbgn} = 
#   { bait=>[$search_id, ], bait_peptides => [0, 21, 33], 
#     prey=>[$search_id, ] (not including bait-bait counts)
#     stats{id}{$experiment} = {$search_id => 1}

# count the number of experiments with X bait peptides
sub countExperiments {
    my @experiments = @_;
    
    my %expCnt; #expCnt{0} = count of experiments where there were 0 bait cnts
    # expCnt{1} is the same
    # expCnt{2} is anything above 1
    for my $exp (@experiments) {
	for my $prot (keys %{ $exp->{protein} }) {
	    next unless exists $exp->{protein}{$prot}{bait_peptides};
	    for my $pepCount (@{ $exp->{protein}{$prot}{bait_peptides} }) {
		my $p = min ($pepCount, 2);
		$expCnt{$p}++;
	    }
	}
    }

    my ($zero, $one, $two) = map {$expCnt{$_}} qw(0 1 2);
    $zero //=0;
    $one//=0;
    $two//=0;
    my $all = sum ($zero, $one, $two);

    return ($all, $zero, $one, $two);
}

# count the number of unique prey within selected experiments
sub countPrey {
    my @experiments = @_;
    
    my %uniq;
    for my $exp (@experiments) {
	for my $prot (keys %{ $exp->{protein} }) {
	    next unless exists $exp->{protein}{$prot}{prey};
	    $uniq{$prot} = 1;
	}
    }

    return 0+ keys %uniq;
}

# count the number of baits appearing X or more times
sub countBaits {
    my @experiments = @_;

    my %baitCnt; #baitCnt{1} = count of baits with one experiment
    # baitCnt{2} is the same
    # baitCnt{3} is anything above 2
    for my $exp (@experiments) {
	for my $prot (keys %{ $exp->{protein} }) {
	    next unless exists $exp->{protein}{$prot}{bait};
	    my $count = 0+ @{ $exp->{protein}{$prot}{bait} };
	    $count = min($count, 3);
	    $baitCnt{$count}++;
	}
    }

    my ($one, $two, $three) = map {$baitCnt{$_}} qw(1 2 3);
    $one//=0;
    $two//=0;
    $three//=0;
    my $all = sum ($one, $two, $three);
    return ($all, $one, $two, $three);
}

sub writeMatrix {
    my ($out, $hash, $rows, $columns) = @_;

    open my $OUT, ">", $out or die "Can't write to $out. $!";

    say $OUT join "\t", @$columns;
    for my $row (@$rows) {
	say $OUT join "\t", $row, map { $hash->{$_}{$row} } @$columns;
    }

    close $OUT;

    exit;
}
