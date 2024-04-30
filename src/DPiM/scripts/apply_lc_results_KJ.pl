#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Data::Dumper;
use HomeBrew::IO qw(checkExist );
use DpimLib qw(getLineDP4APMS getLineBiopAPMS);

process_cmd_args();

my %experiments;
my $totalToRemove = 0;

my $seqFile = $ARGV[0];
my $apmsFile = $ARGV[1];

open(INPUT, "<$seqFile") or die "Cannot open $seqFile\n";
my $prey;
my $firstInSeries = 1;
while(my $buf = <INPUT>){
    chomp $buf;
    next if ($buf eq '');

    my @tokens = split("\t", $buf);
    if ($buf =~ "^\t") {
        $prey = $tokens[2];
	$firstInSeries = 1;
    } else {
	if ($firstInSeries) {
	    $firstInSeries = 0;
	} else {
	    my $instr_run_id = $tokens[0];
	    my $bait = $tokens[1];
	    die Dumper(\@tokens, $bait, $prey) 
		if ! defined $bait || ! defined $prey;
	    if ($bait ne $prey) {
		$totalToRemove++;
		$experiments{$instr_run_id}{$prey} = 1;
	    }
	}

	    
    }
}

open my $DATA, "<", $apmsFile or die "Cannot read from $apmsFile. $!";
open my $OUT, ">", "apply_lc.removed.out" or die "can't write to apply_lc.removed.out. $!";
#print $header;

print join "\t", qw(search_id bait_ref prey_ref total_peptides sample_date ms_inst_run_id);
my $reader = \&getLineDP4APMS;
if (defined $ARGV[2] && $ARGV[2] eq 'human') {
    $reader = \&getLineBiopAPMS;
    print "\trep\tplate";
}
print "\n";

open my $IN, "<", $apmsFile or die "Can't open $apmsFile. $!";

my %row;
while($reader->(\%row, $IN, 'line')){
    my $ms_inst_run_id = $row{ms_inst_run_id};
    #my $experiment = $tokens[7];
    my $prey_ref = $row{prey_ref};

    # Check to see if the experiment needs to have peptides removed.
    if(exists $experiments{$ms_inst_run_id}{$prey_ref}){
	print $OUT $row{line};
	next;
    } else {
	print $row{line};
    }
}
close($DATA);
close($OUT);

warn "$totalToRemove to remove\n";
#warn Dumper(\%experiments);

sub process_cmd_args {
    if(@ARGV < 2 || @ARGV > 3 || (@ARGV == 3 && $ARGV[2] ne 'human')){
	die "Usage: $0 seqContam in.apms [ optional: 'human' ] > out\n";
    }
    checkExist('f', $_) for @ARGV[0,1];
}
