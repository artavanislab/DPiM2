#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Temp qw(tempfile);
use List::Util qw(shuffle);
use Test::More;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineDP4_1APMS);

# find experiments that have identical prey TSC

# algorithm: for every experiment, find the prey with the most TSC
#   create a hash with abundant TSC as the key and list of experiments as value
#   for each prey, compare each experiment, find duplicates
#   once a duplicate class has been identified, select the right one
#
# The right one will:
#    not have FBgn0000000
#    will have a decimal in its ms_inst_run_id

die "this script will need to be modified to reflect the fact that duplicate ms_inst_run_id's have been changed.  start keying on bait_ref instead of rid.\n";

my %opts = getCommandLineOptions();

if (exists $opts{test}) {
    test();
}

{
    my $in = $opts{in};
    my $out = $opts{out};

    my @dupIDs = findDupIDs($in);

    my %apms; 
    readAPMS(\%apms, $in, \@dupIDs);

    my %dupRuns; # dupRuns{ms_inst_run_id}{search_id} = 1 iff this run
    # duplicates another run
    # note that this should leave one of the duplicates alone!
    findDupRuns(\%dupRuns, \%apms);
    
    filterAndReprint($in, $out, \%dupRuns);
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

    my $usage = "usage: $0 -in input -out output < -test >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "test");
    die $usage unless exists $opts{test} || exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    return %opts if exists $opts{test};
    
    checkExist('f', $opts{in});

    return %opts;
}

sub findDupIDs {
    my ($in) = @_;

    my %ids; # ids{ms_inst_run_id}{search_id}==1 iff there is an experiment with
    # these id's
    
    open my $IN, "<", $in or die "can't read $in. $!";
    my %row;
    while (getLineDP4_1APMS(\%row, $IN)) {
	$ids{$row{ms_inst_run_id}}{$row{search_id}} = 1;
    }
    close $IN; 
    return grep { 1 < keys %{$ids{$_}} } keys %ids;
}

# ret{ms_inst_run_id}{search_id}{prey_ref} = tsc
# return only experiments in @$runIDList
sub readAPMS {
    my ($ret, $inFile, $runIDList) = @_;

    my %runIDs = map {$_ => 1} @$runIDList;

    my %row;
    open my $IN, "<", $inFile or die "can't read $inFile. $!";
    while (getLineDP4_1APMS(\%row, $IN)) {
	my $rid = $row{ms_inst_run_id};
	next unless exists $runIDs{$rid};
	my $sid = $row{search_id};
	$ret->{$rid}{$sid}{$row{prey_ref}} = $row{total_peptides};
    }	
    close $IN;
    
    return;
}


# dupRuns{ms_inst_run_id}{search_id} = 1 iff this run duplicates another run
sub findDupRuns {
    my ($ret, $apms) = @_;

    for my $rid (keys %$apms) {
	my @sids = sort {$a <=> $b} keys %{$apms->{$rid}};
	do {
	    my $prevID = shift @sids;
	    next if exists $ret->{$rid}{$prevID};
	    for my $sid (@sids) {
		next if exists $ret->{$rid}{$sid};
		$ret->{$rid}{$sid} = 1 if
		    runsAreEqual($apms->{$rid}{$prevID}, $apms->{$rid}{$sid});
	    }
	} while (1 > @sids)
    }
}

# return true if every prey=>tsc pair in the two experiments match
sub runsAreEqual {
    my ($run1, $run2) = @_;

    my @prey1 = keys %$run1;
    return undef if @prey1 != keys %$run2;
    for my $prey (@prey1) {
	return undef if ! exists $run2->{$prey};
	return undef if $run1->{$prey} != $run2->{$prey};
    }

    return 1;
}

sub filterAndReprint {
    my ($in, $out, $dupRuns) = @_;

    my $nDupes = 0;
    for my $rid (keys %$dupRuns) {
	for my $sid (keys %{ $dupRuns->{$rid} }) {
	    $nDupes++;
	}
    }

    open my $IN, "<", $in or die "can't read $in. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 removed $nDupes duplicate runs";
    while (<$IN>) {
	next if /^#/;
	my @spl = split;
	my ($sid, $rid) = ($spl[0], $spl[-2]);
	next if exists $dupRuns->{$rid}{$sid};
	print $OUT $_;
    }
    close $OUT; 
    close $IN; 
}
# make a bunch of experiments, see if you remove the right ones
sub test {

    my @dupes = qw(1 2 3 4);
    my $nUniq = 5;
    my $nPrey = 10;
    my @rids = shuffle 1..10000;

    ## tests: 
    ##  * 1 test findDupIDs
    ##  * for each ms_run_id, test that you get the right number of unique runs
    ##  * test that you get the exact data back
    plan tests => 1 + @dupes + $nUniq + 1;
    

    my @trueDupRids;
    
    my %allData;
    my %nDupes;
    my @uniques; # unique experiments
    for my $n (@dupes) {
	my $rid = sprintf("f%05d", shift @rids);
	push @trueDupRids, $rid;
	
	$nDupes{$rid} = $n;
	my $expt = randExpt($nPrey);
	my %newExpt = map { $_ => 1+$expt->{$_} } keys %$expt; # almost the same
	my @expts = ($expt) x $n;
	push @expts, \%newExpt;
	push @uniques, $expt;
	push @uniques, \%newExpt;
	$allData{$rid} = \@expts;
    }

    for (1..$nUniq) {
	my $rid = sprintf("f%05d", shift @rids);
	my $expt = randExpt($nPrey);
	$allData{$rid} = [$expt];
	push @uniques, $expt;
    }

    
    my ($BASE, $baseFile) = tempfile();
    say "writing unfiltered input to $baseFile";
    writeAPMS(\%allData, $BASE);
    close $BASE;
    my @allRids = keys %allData;

    
    my @gotDupRids = findDupIDs($baseFile);
    {
	my %got = map { $_ => 1 } @gotDupRids;
	my %true = map { $_ => 1 } @trueDupRids;
	is_deeply(\%got, \%true, "find the right duplicate ids");
    }
    
    
    my %apms; 
    readAPMS(\%apms, $baseFile, \@gotDupRids);

    my %dupRuns; # dupRuns{ms_inst_run_id}{search_id} = 1 iff this run
    # duplicates another run
    # note that this should leave one of the duplicates alone!
    findDupRuns(\%dupRuns, \%apms);

    my ($FH, $outFile) = tempfile();
    close $FH;
    say "writing filtered output to $outFile";
    filterAndReprint($baseFile, $outFile, \%dupRuns);

    my %gotAPMS;
    readAPMS(\%gotAPMS, $outFile, \@allRids);
    my @gotRuns;
    for my $rid (keys %gotAPMS) {
	my @sids = keys %{ $gotAPMS{$rid} };
	my $nRun = 0+ @sids;
	if (exists $nDupes{$rid}) {
	    is($nRun, 2, "does $rid have 2 runs?");
	} else {
	    is($nRun, 1, "does $rid have 1 run?")
	}
	push @gotRuns, map { $gotAPMS{$rid}{$_} } @sids;
    }

    my @gotReordered;
    for my $u (@uniques) {
	my @tmpRuns;
	while (@gotRuns && ! runsAreEqual($u, $gotRuns[0])) {
	    push @tmpRuns, shift @gotRuns;
	}
	push @gotReordered, shift @gotRuns;
	push @tmpRuns, @gotRuns;
	@gotRuns = @tmpRuns;
    }
    is_deeply(\@uniques, \@gotReordered, "match every unique run");
    
    exit;
}

sub randExpt {
    my $nPrey = shift;
    my @tsc = map 
    { 1 + int(rand(20)) } 1..$nPrey;
    my @randInts = shuffle 1..1000;
    my @ids = map { sprintf "FBgn%07d", $_ } @randInts[0..$#tsc];
    return { map { $ids[$_] => $tsc[$_] } 0..$#tsc };
}
    
sub writeAPMS {
    my ($apms, $OUT) = @_;

    say $OUT join "\t", qw(search_id       bait_ref        prey_ref        total_peptides  sample_date    ms_inst_run_id logp);
    my $sid = 1;
    my $date = "2008-01-31";
    for my $rid (keys %$apms) {
	my $bait = sprintf "FBgn%07d", rand(10000);
	for my $expt (@{ $apms->{$rid} }) {
	    for my $prey (keys %$expt) {
		say $OUT join "\t", $sid, $bait, $prey, $expt->{$prey}, $date
		    , $rid, rand();
	    }
	    $sid++;
	}
    }
    return;
}
