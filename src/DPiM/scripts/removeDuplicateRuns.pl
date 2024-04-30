#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineAPMS getLineRawAPMS getLineDP4APMS getLine2015DP3APMS);

# somehow, it's happening that the same experiment appears more than once in the
#   data.  Identify such duplicates and remove them.
# Note that this script will not remove duplicates if they are adjacent
#   within the file.  It expects duplicates to appear separately in the file
#
# The input is expected to have the following format:
# tap_id  ms_inst_run_id  user    search_id       sample_date     total_peptides  unique_peptides bait_ref        pref_ref

# the algorithm works like this: 
# I have the search id for this line
# If I haven't seen this searchID:
#   Set this sID "active" so that I'll know this is a new experiment.
#   report it.
# If I have seen it and it is active:
#   report it.
# If I've seen this searchID before and it's not active:
#   Do not report it. (Do not set active.)
#
# How does this work when I'm looking at the first occurrence:
# * The first time I see it, it's new.  set active.
# * subsequent (contiguous in file) entries remain active.
# * then I find the next sid.
#
# How does this work when I'm looking at a subsequent occurrence:
# * The first entry, it's not new, but I'm not active
# * this remains the case until I find a new sid
#



my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};

    my ($reader, $header);
    if ($mode eq 'dp4') {
	$header = join "\t", qw(search_id       bait_ref        prey_ref        total_peptides  sample_date ms_inst_run_id);
	$reader = \&getLineDP4APMS;
    } elsif ($mode eq 'dp3') {
	$header = "name    bait    search_id       reference       abundance       gene";
	$reader = \&getLine2015DP3APMS;
    } elsif ($mode eq 'raw') {
	$header = "tap_id  ms_inst_run_id  user    search_id       sample_date     total_peptides  unique_peptides bait_ref        pref_ref";
	$reader = \&getLineRawAPMS if $mode eq 'raw';
    } elsif ($mode eq 'clean') {
	$header = "tap_id search_id       sample_date     total_peptides  bait_ref        prey_ref";
	$reader = \&getLineAPMS if $mode eq 'clean';
    } else {
	die "mode '$mode' not implemented";
    }

    open my $IN,  "<", $in or die "Cannot open $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";

    say $OUT $header;
    my ($prevSID, %seen, %row);
    while($reader->(\%row, $IN, 'line')) {
	my $sID = $row{search_id};
	if (exists $seen{$sID} && $prevSID ne $sID) {
	    #warn $row{line};
	    #warn "\tsID = $sID, prevSID = $prevSID\n";
	    #warn "\n";
	    next;
	} 
	$seen{$sID} = 1;
	$prevSID = $sID;
	print $OUT $row{line};
    }

    exit;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my @modes = qw(dp4 raw clean dp3);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output  < $modeString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    checkExist('f', $opts{in});

    return %opts;
}

