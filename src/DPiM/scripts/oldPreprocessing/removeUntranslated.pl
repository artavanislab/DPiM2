#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsHash);
use DpimLib qw(getLineDP4APMS getLineBiopAPMS);

# remove FBgn's referring to genes that are never translated
# (functionally, all that this script does is remove named FBgn's)
# see findUntranslatedFBgn.pl 

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $removeFile = $opts{remove};

    my %remove = readColsHash($removeFile, [qw(id id)]);
    
    open my $IN, "<", $in or die "can't read $in. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";

    my $reader = \&getLineDP4APMS;
    $reader = \&getLineBiopAPMS if $opts{mode} eq 'human';
    my %row;
    while ($reader->(\%row, $IN, 'line')) {
	next if exists $remove{$row{prey_ref}};
	next if exists $remove{$row{bait_ref}};
	print $OUT $row{line};
    }
    close $IN;
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
	remove => '/home/glocke/DPiM/dpim4/withInstr/apmsData/untranscribed.FBgn.geneSymbol.txt'
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $modeString $defaultString ".
	">\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "remove=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{mode} //= $modes[0];
    die "you must select one of the following modes: "
	, join(", ",keys(%modes)), "\n" if ! exists $modes{$opts{mode}};

    checkExist('f', $opts{in});
    checkExist('f', $opts{remove});

    return %opts;
}

