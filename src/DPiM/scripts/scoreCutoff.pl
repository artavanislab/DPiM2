#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(readHS);

# parse hyperspec output
# select edges above a score cutoff and output

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $cutoff = $opts{cutoff};
    my $human = $opts{mode} eq 'human';
    
    my $cmp = sub {
	my $score = shift;
	return $score >= $cutoff;
    };
    if ($cutoff < 0) {
	$cutoff = -$cutoff;
	$cmp = sub {
	    my $score = shift;
	    return $score < $cutoff;
	};
    }
    if (exists $opts{countlines}) {
	my $cnt = 0;
	$cmp = sub {
	    $cnt++;
	    return $cnt <= $cutoff;
	}
    }
    
    my $SCORE = 2;
    
    open my $IN, "<", $in or die "can't read from $in. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    unless (exists $opts{noheader}) {
	say $OUT "protein1\tprotein2\tscore";
    }
    while (my $line = readHS($IN, $human)) {
	$_ = $line;
	chomp;
	my @spl = split;
	print $OUT $line if $cmp->($spl[$SCORE]);
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
	cutoff => 73.06,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in hyper.o -out output < $modeString ".
	"$defaultString -noheader -countlines >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "cutoff=f", "noheader", 
	       "countlines");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{mode} //= $modes[0];
    die "you must select one of the following modes: "
	, join(", ",keys(%modes)), "\n" if ! exists $modes{$opts{mode}};

    checkExist('f', $opts{in});

    return %opts;
}
