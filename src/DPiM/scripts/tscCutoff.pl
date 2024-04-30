#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempfile);
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsRef);
use DpimLib qw(getLineDP4APMS);


# remove experiments with TSC in the top/bottom X percentile

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $experimentStats = $opts{expstats};
    my $fraction = $opts{fraction};
    my $filteredOut = $opts{filtered};
    my $statsScript = $opts{statsscr};

    if (! defined $experimentStats) {
	my $FH;
	($FH, $experimentStats) = tempfile();
	my $cmd = "$statsScript -in $in -out $experimentStats";
	say $cmd;
	system($cmd);
    }
    
    my %filter = makeFilter($experimentStats, $fraction, $opts{choptop});

    open my $IN, "<", $in or die "can't read from $in. $!";

    my $header = join "\t", qw(search_id bait_ref prey_ref total_peptides sample_date ms_inst_run_id);
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# removed extreme $fraction fraction from $in";
    say $OUT $header;

    my %FILT;
    if (defined $filteredOut) {
	my $hiOut = "$filteredOut.toohigh";
	open my $HI, ">", $hiOut or die "can't write to $hiOut. $!";
	$FILT{1}=$HI;
	my $loOut = "$filteredOut.toolow";
	open my $LO, ">", $loOut or die "can't write to $loOut. $!";
	$FILT{-1}=$LO;

	say $_ $header for values %FILT;
    }

    my %row;
    while (getLineDP4APMS(\%row, $IN, 'line')) {
	my $f = $filter{$row{search_id}};
	if ($f == 0) {
	    print $OUT $row{line};
	} elsif(exists $FILT{1}) {
	    print {$FILT{$f}} $row{line};
	} 
    }
    close $IN;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	fraction => 0.05, 
	statsscr => ($ENV{DPSCR} // "full_path")."/baitPullDown.pl",
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString -expstats ".
	"statsBySID -filtered optional.out -choptop > \n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "statsscr=s", "fraction=f", 
	       "expstats=s", "filtered=s", "choptop");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{expstats}) if exists $opts{expstats};
    checkExist('f', $opts{statsscr}) if ! exists $opts{expstats};

    return %opts;
}

# %ret{search_id} = 0 if pass, 1 if too high, -1 if too low
# as of 09-20-2016, do not remove experiments for having too many TSC, only
# chop the bottom.  if $chopTop is set, then you do remove them
sub makeFilter {
    my ($file, $fraction, $chopTop) = @_;
    
    my @cols = qw(search_id tsc);
    my @read;
    readColsRef(\@read, $file, \@cols);

    @read = sort { $a->{tsc} <=> $b->{tsc} } @read;

    my %ret = map { $_->{search_id} => 0 } @read; # defaults to "pass"

    my $nDropTarget = sprintf("%.0f", $fraction * @read);

    my $loBreak = findBreakpoint(\@read, $nDropTarget, 'lo');
    for my $i (0..$loBreak) {
	$ret{$read[$i]{search_id}} = -1;
    }

    if ($chopTop) {
	my $hiBreak = findBreakpoint(\@read, $nDropTarget, 'hi');
	for my $i ($hiBreak..$#read) {
	    $ret{$read[$i]{search_id}} = 1;
	}
    }
    
    return %ret;
}

# find the index for the min(max) filtered search_id
# we want to remove $nDropTarget experiments, but we don't want to some 
#  experiments with a given TSC and not others.  So, we find the TSC at which
#  the number of experiments removed with 
sub findBreakpoint {
    my ($read, $nDropTarget, $mode) = @_;
    
    my $i;
    if ($mode eq 'lo') {
	$i=$nDropTarget-1;
    } elsif ($mode eq 'hi') {
	$i = @$read - $nDropTarget;
    } else {
	die "findBreakpoint unknown mode '$mode'";
    }

    my $ret;

    my $tsc = $read->[$i]{tsc};
    my ($bottom, $top) = ($i, $i);
    if ($mode eq 'lo') {
	# seek *highest* indices with a given tsc
	$bottom-- while $read->[$bottom]{tsc} == $tsc;
	$top++    while $read->[$top+1]{tsc} == $tsc;
	$ret = $bottom;
	$ret = $top if ($i - $bottom) > ($top - $i);
	# for example $i = 50, $top = 55, $bottom = 40;
	#   ($i - $bottom = 10) > ($top - $i = 5) so 
	#   therefore, switch to $top
    } elsif ($mode eq 'hi') {
	# seek *lowest* indices with a given tsc
	$bottom-- while $read->[$bottom-1]{tsc} == $tsc;
	$top++    while $read->[$top]{tsc} == $tsc;
	$ret = $top;
	$ret = $bottom if ($top - $i) > ($i - $bottom);
    }
    say "$mode";
    say "\tstart = $i";
    say "\tbottom - i = $bottom, sid = $read->[$bottom]{search_id}, tsc = $read->[$bottom]{tsc}";
    say "\ttop - sid = i = $top, $read->[$top]{search_id}, tsc = $read->[$top]{tsc}";
    say "\tselected $read->[$ret]{search_id}";
    
    return $ret;
}
