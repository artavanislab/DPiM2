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

    my %aveLen;
    my @cols = qw(fbgn avg_length);
    @cols = qw(entrez avg_length) if $opts{mode} eq 'human';
    readColsHashRef(\%aveLen, $aveLenFile, \@cols);

    my $reader = \&getLineDP4APMS;
    $reader = \&getLineBiopAPMS if $opts{mode} eq 'human';

    my (%missing, %baits, %row);
    open my $IN, "<", $in or die "can't read $in. $!";
    while ($reader->(\%row, $IN)) {
	my $bait = $row{bait_ref};
	$baits{$bait}{$row{search_id}} = 1;
	my $prey = $row{prey_ref};
	next if $prey eq $bait;
	if (! exists $aveLen{$prey}) {
	    $missing{$prey}{prey}++;
	}
    }
    close $IN;

    for my $b (keys %baits) {
	if (! exists $aveLen{$b}) {
	    my $n = 0+ keys %{$baits{$b}};
	    $missing{$b}{bait}+= $n;
	}
    }

    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT join "\t", qw(id bait prey);
    for my $fb (sort keys %missing) {
	say $OUT join "\t", $fb, $missing{$fb}{bait} // 0
	    , $missing{$fb}{prey} // 0;
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
	avelen => '/home/glocke/DPiM/nsaf/REVdmel-all-translation-r6.07_TAGS_sorted_trEl_vir.aveLen.tsv'
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in in.apms -out missing.txt < -$modeString ".
	"$defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "avelen=s");
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
