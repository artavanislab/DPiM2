#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsRef );
use DpimLib qw(getLineMayAPMS);

# ms_inst_run_id is supposed to be a unique key, but it isn't.
# There are id's corresponding to multiple experiments.
# This is mainly a problem because of LC carryover correction.
# These id's are treated as sequential by the LCCC code, but the experiments
#   sharing these non-unique keys are out of order (at least some of them).
# Bob and Marina put together some spreadsheets explaining how to re-order these
#   experiments.
# This script executes these corrections.

# For each (contiguous block of) to-be-reordered experiment(s), find the nearest
#   preceding experiment that is correctly ordered.
# Use the id for that preceding experiment as the base for the TBR experiments,
#   adding decimals to it indicating the order.

# This means finding the correct order for every non-duplicate-id experiment 
#   first - keyed by date.
# It so happens that the dates for the DPIM1 era experiments are not recorded
#   accurately, but we're in luck because the only duplicates from from 2014 and
#   2015.

# The end product of this script is a table showing how id strings have been
#   changed, and, of course, a corrected version of the apms data file


my %opts = getCommandLineOptions();

my $NACode = 'dupeNotChanged';

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $correctionFile = $opts{correct};
    my $logOut = $opts{log};
    my $verboseOut = $opts{verboselog};
    
    # $correct{{YYYY-mm-dd}{search_id} = {
    #   ms_inst_run_id => '',
    #   tap_id => '', 
    #   order_per_date => '',
    # }
    
    my %correct = readCorrections($correctionFile);
	
    # byDate{YYYY-mm-dd}{search_id} = {
    #   ms_inst_run_id => '',
    #   tap_id => '', 
    #   id_string => '',
    # }
    my %byDate = readIdStringsByDate($in);

    # $change{search_id} = { from => "prev run_id", to => "new run_id" }
    #   these values will be undefined if search_id was found in the file
    #   but no change is to be made.
    my %change = whatChanges(\%byDate, \%correct, $verboseOut);
    
    makeChanges($in, $out, \%change);

    if (defined $logOut) {
	open my $OUT, ">", $logOut or die "Can't write to $logOut: $!";
	say $OUT join "\t", qw(search_id prev_rid new_rid);
	for my $sid (sort {$a <=> $b} keys %change) {
	    say $OUT join "\t", $sid ,$change{$sid}{from}, $change{$sid}{to};
	}
	close $OUT;
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	in => '/home/glocke/DPiM/augRemap/apmsData/DPiM3_r607_protein_views.out_LDAmin7_PA_PSA01_160824',
	correct => '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable3.tap.correct1.multi_RAO.txt',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString -log changesHere ".
	"-verbose evenMoreInformationHere >\n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "in=s", "correct=s", "log=s", "verboselog=s");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{correct});

    return %opts;
}

# $ret{{YYYY-mm-dd}{search_id} = {
#   ms_inst_run_id => '',
#   tap_id => '', 
#   order_per_date => '',
# }
sub readCorrections {
    my ($correctionFile) = @_;

    my %ret;
    my @correct;
    my @cols = qw(sample_date search_id ms_inst_run_id tap_id order_per_date);
    
    readColsRef(\@correct, $correctionFile, \@cols);
		
    my @capture = @cols[2..$#cols];
    for my $row (@correct) {
	next unless defined $row->{search_id};
	## stupid Excel reformats dates, so I have to re-reformat them
	## Ever see it turn the gene Sept1 into 9/1/2016?
	my @spl = split '/', $row->{sample_date};
	die Dumper($row, \@spl) if 3 != @spl || !defined $spl[0] || !defined $spl[1] || !defined $spl[2];
	my $date = join "-", $spl[2], sprintf("%02d", $spl[0])
	    , sprintf("%02d", $spl[1]);
	$ret{$date}{$row->{search_id}} = {
	    map { $_ => $row->{$_} } @capture };
    }
    return %ret;
}

# for every experiment:
# ret{YYYY-mm-dd}{search_id} = {
#   ms_inst_run_id => '',
#   tap_id => '', 
#   id_string => '',
# }
sub readIdStringsByDate {
    my ($in) = @_;
    my %ret;

    my @capture = qw(ms_inst_run_id tap_id id_string);
    open my $IN, "<", $in or die "can't read $in. $!";
    my %row;
    while (getLineMayAPMS(\%row, $IN)) {
	## convert date from MMDDyy to YYYY-MM-DD
	if (length($row{sample_date})==5) {
	    # in a few cases, the date string will omit the leading '0' for
	    # months prior to october.  Grr...
	    $row{sample_date} = '0'.$row{sample_date};
	}
	my @spl = $row{sample_date} =~ /../g;
	$spl[-1] = "20".$spl[-1];
	unshift @spl, pop @spl;
	my $date = join "-", @spl;
	$ret{$date}{$row{search_id}} = { map { $_ => $row{$_} } @capture };
    }
    close $IN; 
    return %ret;
}

# $ret{search_id} = { from => "prev id_string", to => "new id_string" }
#   these values will be undefined if search_id was found in the file
#   but no change is to be made.
sub whatChanges {
    my ($byDate, $correct, $verboseFile) = @_;

    my $VERBOSE;
    if (defined $verboseFile) {
	open $VERBOSE, ">", $verboseFile 
	    or die "Can't write to $verboseFile: $!";
    }
    my %ret;

    for my $date (keys %$correct) {
	##next unless $date eq '2015-11-19';
	say $date;
	if (! exists $byDate->{$date}) {
	    warn "debug silliness at $date!";
	    next;
	}
	say $VERBOSE "#\n# $date" if defined $VERBOSE;
	changesOnDate(\%ret, $byDate->{$date}, $correct->{$date}, $VERBOSE);
    }
    
    return %ret;
}

# $change{search_id} = { from => "prev id_string", to => "new id_string" }
sub changesOnDate {
    my ($ret, $input, $corrections, $VERBOSE) = @_;

    ## all search_id's corresponding to experiments whose run_id's are not
    ## duplicated
    my @nonDupes = grep { ! exists $corrections->{$_} } keys %$input;
    my %dupes;
    {
	my @dupes = grep { exists $corrections->{$_} } keys %$input;
	my @newDupes;
	for my $sid (@dupes) {
	    if ($corrections->{$sid}{order_per_date} =~ /NA/) {
		push @nonDupes, $sid; 
		$ret->{$sid} = { from => $input->{$sid}{ms_inst_run_id},
				 to => $NACode};
		##die Dumper($corrections);
	    } else {
		push @newDupes, $sid;
	    }
	}
	@dupes = sort { $corrections->{$a}{order_per_date} <=>
			    $corrections->{$b}{order_per_date} } @newDupes;
	for my $sid (@dupes) {
	    my $order = $corrections->{$sid}{order_per_date}-1;
	    $dupes{$order} = $sid;
	}
    }
    @nonDupes = sort {$input->{$a}->{ms_inst_run_id} cmp 
			  $input->{$b}->{ms_inst_run_id}} @nonDupes;

    my @insertAt = sort {$a <=> $b} keys %dupes;
    if ($insertAt[-1]+1 > @nonDupes+@insertAt) {
	die "There appears to be an annotation error";
    }
    ##for my $i (@insertAt) {
##	say join "\t", $i, $dupes{$i};
  ##  }
    ##die Dumper(\@nonDupes);

    
    ##
    ## how to use splice for random-access insertion
    ## splice(@a, $index, 0, $what);
    ## inserts $what prior to the element $index

    ## this code would probably be a lot cleaner if splice were used
    my @reordered;
    while (@nonDupes || @insertAt) {
	## insert into reordered until you get to a place in the order where
	## a duplicate experiment is to be inserted
	while (@nonDupes && (! @insertAt || @reordered < $insertAt[0])) {
	    push @reordered, shift @nonDupes;
	}

	my $prevRID;
	if (@reordered) {
	    my $prevSID = $reordered[-1];
	    $prevRID = $input->{$prevSID}{ms_inst_run_id};
	} else {
	    ## reordered experiments go in the front, so we put these in front
	    ## of the first non-reordred experiment
	    if (! @nonDupes) {
		## all experiments on October 26, 2015 are duplicates
		## hard-code an exception for this case
		die "what day is this??".Dumper($corrections) 
		    unless $dupes{0} == 252408;
		$prevRID = 'w35039';
	    } else {
		my $firstRID = $input->{$nonDupes[0]}{ms_inst_run_id};
		$firstRID =~ /^(\D+)(\d+)$/ 
		    or die "can't parse run_id $nonDupes[0]";
		$prevRID = $1.($2+1);
	    }
	}
	while (@insertAt && @reordered == $insertAt[0]) {
	    my $where = shift @insertAt;
	    my $which = $dupes{$where};
	    push @reordered, $which;
	    my $changeFrom = $input->{$which}{ms_inst_run_id};
	    my $changeTo = sprintf "$prevRID.%03d", $where;
	    $ret->{$which}= {from=>$changeFrom, to=>$changeTo};
	}
    }

    if (defined $VERBOSE) {
	dumpReordering($VERBOSE, \@reordered, $ret, $input, $corrections);
    }
    #my @experiments = 
    
    return;
}

sub dumpReordering {
    my ($OUT, $reordered, $change, $input, $corrections) = @_;

    say $OUT join "\t", qw(new_order sid prev_rid new_rid order_per_date);
    my $i = 1;
    for my $sid (@$reordered) {
	my $rid = $input->{$sid}{ms_inst_run_id};
	my $newRid = "same";
	my $opd = "not";
	if (exists $corrections->{$sid}) {
	    $opd = $corrections->{$sid}{order_per_date};
	    $newRid = $change->{$sid}{to} // 'poopy';
	}
	say $OUT join "\t",  $i++, $sid, $rid, $newRid, $opd;
    }
    return;
}

sub makeChanges {
    my ($in, $out, $change) = @_;

    
    open my $IN, "<", $in or die "can't read $in. $!";
	open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 reordered $opts{in} according to $opts{correct}";
    my %row;
    while (getLineMayAPMS(\%row, $IN, 'line')) {
	my $sid = $row{search_id};
	if (exists $change->{$sid} && $change->{$sid}{to} ne $NACode) {
	    $row{line} =~ s/$change->{$sid}{from}/$change->{$sid}{to}/
		or die "can't alter $row{line}\n", Dumper($change->{$sid})
	}
	print $OUT $row{line};
    }
    close $IN;
    close $OUT; 
}
