#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);
use DpimLib qw(getLineDP4APMS getLineDP4_1APMS);

# concatenate files .dp4 files

my %opts = getCommandLineOptions();

{
    my $listFile = $opts{in};
    my $out = $opts{out};

    my @files = readList($listFile);
    checkExist('f', $_) for @files;

    
    my @cols = qw(search_id       bait_ref        prey_ref        total_peptides  sample_date     ms_inst_run_id  logp);
    my $reader = \&getLineDP4_1APMS;
    if ($opts{nologp}) {
	pop @cols;
	$reader = \&getLineDP4APMS;
    }
    
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 concatenated files in $listFile";
    say $OUT join "\t", @cols;
    my %row;
    for my $in (@files) { 
	open my $IN, "<", $in or die "can't read $in. $!";
	while ($reader->(\%row, $IN, 'line')) {
	    print $OUT $row{line};
	}
	close $IN;
    }
    close $OUT; 
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

    my $usage = "usage: $0 -in in.list -out output < -nologp >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "nologp");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
