#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineDP4APMS);

# pull out GO annotations from infile

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    
    open my $IN, "<", $in or die "can't read $in. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 pulled columns 2 and 4 (1-based) from $in";
    say $OUT join "\t", qw(uniprot go);

    my %exceptions = (
	'NOT'=>1,
	'NOT|contributes_to'=>1,
	'contributes_to'=>1,
	'colocalizes_with'=>1,
	'NOT|colocalizes_with'=>1,
	);
    my %uniq;
    while (<$IN>) {
	next if /^!/;
	my @spl = split;
	my ($uni, $go) = ($spl[1], $spl[3]);
	next unless $go =~ /^GO:\d+$/;
	#$go = $spl[4] if exists $exceptions{$go};
	#die "go '$go' doesn't look like go" if 
	say $OUT join "\t", $uni, $go;
	$uniq{$uni} = 1;
    }
    close $IN; 
    close $OUT;

    say for sort keys %uniq;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	in => '/home/glocke/DPiM/human/GO/goa_human.gaf'
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "in=s");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
