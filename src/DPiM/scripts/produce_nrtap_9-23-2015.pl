#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw[max];
use Data::Dumper;
use DateTime;
use DateTime::Format::Strptime;
use v5.10; 

use DpimLib qw(getLineAPMS);

sub process_cmd_args();
process_cmd_args();

my $apmsFile = $ARGV[0];

################################################################
# search for the best tap id for baits appear in different runs
################################################################

## GL update 9-23-2015
# previously, there was no protocol for what to do when multiple experiments
# pulled down the same TSC for the bait.  Now, we prefer the experiment with the
# highest total TSC for all prey in such ties.
# FURTHERMORE, and perhaps more importantly, we have identified a period of from
# April '09 to May '10 where the AP-MS runs appear very low quality.  We now
# prefer to exclude such runs where alternatives are possible.  This preference
# has *higher priority* than the TSC-based filters.



my %experiments; # $exp{bait_ref}{exp_id}={date=DateTime, bait=TSC, tot=TSC}
my $parseDate = DateTime::Format::Strptime->new(
    pattern   => '%Y-%m-%d');

open my $IN, "<", $apmsFile or die "Cannot open $apmsFile. $!";
my %row;
while(getLineAPMS(\%row, $IN)) {
    my $bait = $row{bait_ref};
    my $prey = $row{prey_ref};
    my $search_id = $row{search_id};
    my $tap_id = $row{tap_id};
    my $tsc = $row{total_peptides};
    
    my $expKey = $tap_id."_".$search_id;

    ## grouping baits
    
    # Bob: FBgn0000042, FBgn0000043, FBgn0000044, FBgn0000045, FBgn0000046 and FBgn0000047, only pick the one run with the highest total spectral count number
    # this means we need to convert all above ids into A id, pick 'FBgn0000042'
    if ($bait =~ /FBgn000004[2,3,4,5,6,7]/) {
	$bait = 'FBgn0000042';
    }
    if ($prey =~ /FBgn000004[2,3,4,5,6,7]/) {
        $prey = 'FBgn0000042';
    }
    
    ## THIS IS DONE since we are picking the nrBait
    # Bob: Annexin X is encoded by a single protein-coding gene, FBgn0000084 (CG9579), and there are three FH clones for this in the clone set, FH0025, FH0162, and FH0742.  These are all identical on the sequence level, so only one data set from them should be used

    if (! exists $experiments{$bait}{$expKey}) {
	my $date = $parseDate->parse_datetime($row{sample_date});
	$experiments{$bait}{$expKey} = { date => $date, bait => 0 };
    }

    $experiments{$bait}{$expKey}{tot} += $tsc;
    $experiments{$bait}{$expKey}{bait} += $tsc if $bait eq $prey;

    # There are cases same FBgn number show up multiple times in the same run
    # this is mainly due to multiple isoforms.  That being the case, peptides
    # in such duplicate lines represent distinct, bona fide binding events.
    # Therefore, we will sum them.

    #if($bait eq $prey) {
    #if (defined $experiments{$bait}{}) {
    #$experiments{$bait}{$expKey} = max($experiments{$bait}{$expKey},$tsc); 
    # there are cases same FBgn number show up multiple times in the same run
    #} else {
    #    $experiments{$bait}{$expKey} = $tsc;
    #}
    #}
}
close $IN;

my @tests = ( \&dateTest, \&baitTest, \&totTest );

# if one bait appear in multiple runs (tap_ids) then only keep the tap_id with the highest total peptides count
my @b = keys %experiments;
my %keep_runs;
foreach my $bait (@b) {
    my %runs = %{$experiments{$bait}};

    my @exps = keys %runs;
    for my $t (@tests) {
	@exps = $t->(\@exps, \%runs);
	last if @exps == 1;
    }
    die "can't get down to a single experiment on bait='$bait'" if @exps != 1;
    my $best_id = (sort {$runs{$b} <=> $runs{$a}} keys %runs)[0]; # pick the run with the max total peptides count
    $keep_runs{$best_id} = "";
    #if ($bait eq $debugBait) {
#	die Dumper($best_id, \%runs);
 #   }
}


################################################################
# write out the non-redundant file 
################################################################
open(IN, "<$apmsFile") or die "Cannot open $apmsFile\n";
print "#tap_id\tsearch_id\tsample_date\ttotal_peptides\tbait_ref\tprey_ref\n";
while(my $line = <IN>) {
    chomp $line;
    my @tokens = split("\t", $line);
    my $bait = $tokens[4];
    my $search_id = $tokens[1];
    #if(defined $keep_runs{$search_id}) { # remove all other runs
    #	print $line."\n";
    #}
    my $tap_id = $tokens[0];
    my $expKey = $tap_id."_".$search_id;
    if(defined $keep_runs{$expKey}) { # remove all other runs	
	$line =~ s/_TAG//; # replace ids such as X_TAG with X.
	# replace 'FBgn0266084' by 'FBgn0261259'
	# FBgn0266084 corresponds to gene called Fhos, which also have other secondary FBgn IDs listed below:
	# FBgn0010755
	# FBgn0035925
	# FBgn0035927
	# FBgn0040229
	# FBgn0052025
	# FBgn0052030
	# FBgn0261259
	# only FBgn0261259 is found in our dataset, then replace FBgn0266084 with that
	$line =~ s/FBgn0266084/FBgn0261259/;
	# remove lines contain '#' (ids such as ##FBgnxxx which are peptide matches
	# to the reverse translation of the genome (called decoys) and it is a way 
	# Gygi lab calculates the false discovery rate of peptide matches in the mass
	# spec data. It is safe it ignore them for the HGSCore analysis. 
	# GL note: this is a false positive rate, not FDR

	# fix the mis-annotation
	# Bob mentioned: FH5854 annotated as FBgn0000079 should be corrected to 'FBgn0000078'.
	if ($bait eq 'FBgn0000079' && $tap_id==5854 ) {
	    $line =~ s/FBgn0000079/FBgn0000078/; # TODO: one question left, the identified prey was FBgn0000079 should I change it to FBgn0000078?
	}

	print $line."\n" unless $line =~ /#/;
    }
}


exit(0);

sub process_cmd_args(){
    if(@ARGV != 1){
        print "Usage produce_nrtap.pl [dpim_file_after_apply_lc]\n";
	exit(0);
    }
}


# prefer runs outside of the bad year
sub dateTest {
    my ($exps, $runs) = @_;

    my $start = DateTime->new(year=> 2009, month=> 4, day=> 1);
    my $end   = DateTime->new(year=> 2010, month=> 6, day=> 1);
    
    my @passing = grep { $runs->{$_}{date} < $start && 
			     $runs->{$_}{date} >= $end } @$exps;

    return @passing if @passing > 0;
    return @$exps;
}

# prefer the runs with the highest bait TSC
sub baitTest {
    my ($exps, $runs) = @_;

    my $mx = max map {$runs->{$_}{bait} // die Dumper($exps, $runs)} @$exps;
    
    return grep { $runs->{$_}{bait} == $mx } @$exps;
}

# prefer the runs with the highest total TSC
sub totTest {
    my ($exps, $runs) = @_;

    my $mx = max map {$runs->{$_}{tot}} @$exps;
    
    return grep { $runs->{$_}{tot} == $mx } @$exps;
}
