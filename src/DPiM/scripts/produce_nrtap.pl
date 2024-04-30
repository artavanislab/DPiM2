#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw[max];
use Data::Dumper;
use DateTime;
use DateTime::Format::Strptime;

sub process_cmd_args();
process_cmd_args();

################################################################
# search for the best tap id for baits appear in different runs
################################################################

my %experiments;
open(IN, "<$ARGV[0]") or die "Cannot open $ARGV[0]\n";
while(my $line = <IN>) {
    chomp $line;
    my @tokens = split("\t", $line);
    my $bait = $tokens[4];
    my $prey = $tokens[5];
    my $search_id = $tokens[1];
    my $tap_id = $tokens[0];
    my $tsc = $tokens[3];
    
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


    if($bait eq $prey) {
	#if (defined $experiments{$bait}{$search_id}) {
	#	$experiments{$bait}{$search_id} = max($experiments{$bait}{$search_id},$tsc); 
	# there are cases same FBgn number show up multiple times in the same run
	#} else {
	#	$experiments{$bait}{$search_id} = $tsc;
	#}
	if (defined $experiments{$bait}{$tap_id."_".$search_id}) {
	    $experiments{$bait}{$tap_id."_".$search_id} = max($experiments{$bait}{$tap_id."_".$search_id},$tsc); # there are cases same FBgn number show up multiple times in the same run
        } else {
	    $experiments{$bait}{$tap_id."_".$search_id} = $tsc;
        }
    }
}
close(IN);


# if one bait appear in multiple runs (tap_ids) then only keep the tap_id with the highest total peptides count
my @b = keys %experiments;
my %keep_runs;
foreach my $bait (@b) {
    my %runs = %{$experiments{$bait}};
    my $best_id = (sort {$runs{$b} <=> $runs{$a}} keys %runs)[0]; # pick the run with the max total peptides count
    $keep_runs{$best_id} = "";
    #if ($bait eq $debugBait) {
#	die Dumper($best_id, \%runs);
 #   }
}


if (0) {
################################################################
# TODO: identify preys only have 1 peptide and only show up in one run
################################################################

    open(IN, "<$ARGV[0]") or die "Cannot open $ARGV[0]\n";
    while(my $line = <IN>) {
	chomp $line;


    }
    close(IN);
}

################################################################
# write out the non-redundant file 
################################################################
open(IN, "<$ARGV[0]") or die "Cannot open $ARGV[0]\n";
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
    if(defined $keep_runs{$tap_id."_".$search_id}) { # remove all other runs	
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
	if ($tap_id==5854 && $bait eq 'FBgn0000079') {
	    $line =~ s/FBgn0000079/FBgn0000078/; # TODO: one question left, the identified prey was FBgn0000079 should I change it to FBgn0000078?
	}

	print $line."\n" unless $line =~ /#/;
    }
}


exit(0);

sub process_cmd_args(){
    if($#ARGV != 0){
        print "Usage produce_nrtap.pl [dpim_file_after_apply_lc]\n";
		exit(0);
    }
}


