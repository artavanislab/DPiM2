#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineDP4APMS);

# describe the purpose for this script

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    open my $IN, "<", $in or die "Can't read from $in. $!";

    my %fbgn2symb; # $fbgn2symb{$fbgn} = $symb
    while (my $line = <$IN>) {
	next if length($line) < 3;
	next if $line =~ /^#/;
	next if $line =~ /\\/;
	chomp $line;
	my @spl = split /\t/, $line;
	my $symb = $spl[0];
	my @fbgns = ($spl[1]);
	push @fbgns, (split /,/, $spl[2]) if @spl > 2 && length($spl[2]) > 3;
	for my $fb (@fbgns) {
	    die "can't parse '$line'" if ! defined $fb;
	    if (! exists $fbgn2symb{$fb} || $fbgn2symb{$fb} =~ /^CG\d+$/) {
		# set f2s{fb} if either
		#  - it hasn't been set
		#  - it is currently set to some CG1234
		#    + in which case anything else is no worse
		$fbgn2symb{$fb} = $symb;
	    }
	}
    }
    close $IN;

    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# fbgn2symb map gleaned from $in";
    say $OUT join "\t", qw( fbgn symbol );
    say $OUT join "\t", $_, $fbgn2symb{$_} for keys %fbgn2symb;
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	in => '/home/glocke/DPiM/flybase/fbgn_annotation_ID_fb_2016_01.tsv',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
