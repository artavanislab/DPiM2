#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);
#use DpimLib qw(getLineDP4APMS);

# concatenate a list of biop files
# designed concatenate previously separated files processed in parallel

my %opts = getCommandLineOptions();

{
    my $listFile = $opts{in};
    my $out = $opts{out};
    
    my @files = readList($listFile);
    checkExist('f', $_) for @files;

    system("cp $files[0] $out");

    open my $OUT, ">>", $out or die "can't write to $out. $!";
    for my $f (@files) {
	open my $IN, "<", $f or die "Can't read from $f. $!";
	while (<$IN>) {
	    next if /^#/;
	    next if /^search/;
	    next if length($_) < 4;
	    print $OUT $_;
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

    my $usage = "usage: $0 -in input -out output\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
