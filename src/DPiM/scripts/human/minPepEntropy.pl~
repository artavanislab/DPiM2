#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineBiopAPMS);

# reject any prey unless it has at least two peptides in both replicates
#   

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    my %apms;
    parseAPMS(\%apms, $in);

    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 removed bait-prey line unless it has at least two peptides in both replicates";
    for my $bait (sort {$a <=> $b } keys %apms) {
	for my $p (sort { $apms{$bait}{A}{$b}{tsc} <=> 
			      $apms{$bait}{A}{$a}{tsc} } keys %apms)
	{
	    last if 2 > $apms{$bait}{A}{$p}{tsc};
	    next if 2 > $apms{$bait}{B}{$p}{tsc};
	    print $OUT $apms{$bait}{A}{$p}{line};
	    print $OUT $apms{$bait}{B}{$p}{line};
	}

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

# ret{bait}{A/B}{prey} = {tsc = 1234, line = $line}
sub parseAPMS {
    my ($ret, $in) = @_;

    open my $IN, "<", $in or die "can't read $in. $!";
    my %row;
    while (getLineBiopAPMS(\%row, $IN, 'line')) {
	$ret->{$row->{bait_ref}}{$row->{rep}} = {
	    tsc => $row->{total_peptides},
	    line => $row->{line};
	};
    }

    return;
}
