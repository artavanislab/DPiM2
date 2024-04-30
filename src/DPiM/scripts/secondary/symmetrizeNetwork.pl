#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineAPMS);

# take an (undirected) network and output a network where edges from A->B
# are always accompanied by an edge from B->A

# method: 
# 1. copy the input file to a new file
# 2. append lines where the first two tokens are switched

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    system("cp $in $out") and die "can't copy $in to $out. $!"; 
    # non-zero return values are fail i guess

    open my $IN, "<", $in or die "can't read from $in. $!";
    
    # skip past header
    while (<$IN>) {
	last if /node1/;
	last if /Protein1/;
    }
    
    open my $OUT, ">>", $out or die "can't write to $out. $!";
    while (<$IN>) {
	chomp;
	my @spl = split;
	($spl[1], $spl[0]) = @spl[0..1];
	say $OUT join "\t", @spl;
    }
    close $IN;
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

    my $usage = "usage: $0 -in in.net -out out.symm.net\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
