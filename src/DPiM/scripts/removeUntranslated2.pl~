#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsHashRef);
use DpimLib qw(getLineDP4APMS getLineBiopAPMS);

# find any FBgn's that appear in the data but not in the flybase translated
#   fasta file

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $aveLenFile = $opts{avelen};
    my $logFile = $opts{log};
    
    my %aveLen;
    my @cols = qw(fbgn avg_length);
    @cols = qw(entrez avg_length) if $opts{mode} eq 'human';
    readColsHashRef(\%aveLen, $aveLenFile, \@cols);

    my $reader = \&getLineDP4APMS;
    $reader = \&getLineBiopAPMS if $opts{mode} eq 'human';

    my (%missing, %baits, %row);
    open my $IN, "<", $in or die "can't read $in. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 removed proteins in $in that don't appear in $aveLenFile";
    say $OUT join "\t", qw(search_id       bait_ref        prey_ref        total_peptides  sample_date    ms_inst_run_id);
    while ($reader->(\%row, $IN, 'line')) {
	my $bait = $row{bait_ref};
	$baits{$bait}{$row{search_id}} = 1;
	my $prey = $row{prey_ref};
	if (! exists $aveLen{$prey}) {
	    $missing{$prey}{prey}++;
	} elsif (exists $aveLen{$bait}) {
	    print $OUT $row{line};
	}
    }
    close $IN;
    close $OUT;

    for my $b (keys %baits) {
	if (! exists $aveLen{$b}) {
	    my $n = 0+ keys %{$baits{$b}};
	    $missing{$b}{bait}+= $n;
	}
    }

    if (defined $logFile) {
	open my $LOG, ">", $logFile or die "can't write to $logFile. $!";
	say $LOG join "\t", qw(id bait prey);
	for my $fb (sort keys %missing) {
	    say $LOG join "\t", $fb, $missing{$fb}{bait} // 0
		, $missing{$fb}{prey} // 0;
	}
	close $LOG;
    }
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
	avelen => '/home/glocke/DPiM/nsaf/REVdmel-all-translation-r6.07_TAGS_sorted_trEl_vir.aveLen.tsv'
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in in.apms -out sanitized.apms < $modeString ".
	"$defaultString -log missing.tsv >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "avelen=s", "log=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{mode} //= $modes[0];
    die "you must select one of the following modes: "
	, join(", ",keys(%modes)), "\n" if ! exists $modes{$opts{mode}};

    checkExist('f', $opts{in});
    checkExist('f', $opts{avelen});

    return %opts;
}
