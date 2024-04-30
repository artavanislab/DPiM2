#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Storable;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineDP4APMS);

# in particular, I am interested to find whether edges in the max support 
#   network all derive from experiments with single-replicate-baits

# so, we want to find some summary statistics:
#   for each edge in network:
#      how many co-occurences?
#      how many co-occurences in single-bait replicates?

my %opts = getCommandLineOptions();

{
    my $net = $opts{net};
    my $apmsFile = $opts{apms};
    my $out = $opts{out};

    my %bait; # { occur, tsc = [], baitTsc = [] }
    my %prey; # { occur, occurSB, baits = [] }
    my %pbp; # pbp{prey1}{prey2} = x
    # x = { support, coOccur, coOccurSB, baitrepl=y, dates=z }
    # y = for every co-occurrence, list the number of replicates this bait had
    # z = list the dates of every co-occurrence
    {
	my $apms;
	# apms{bait}{search_id} = [ {prey_ref=>'fbgn', 'total_peptides'=>n},...]
	if (-e "$apmsFile.storable") {
	    $apms = retrieve("$apmsFile.storable");
	} else {
	    $apms = apmsRound1($apmsFile);
	}

	rawStats($apms, \%bait, \%prey);
    }    
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

    my $usage = "usage: $0 -net net -apms applyLC -out output\n";

    my %opts = ();
    GetOptions(\%opts, "net=s", "apms=s", "out=s");
    die $usage unless exists $opts{net} && exists $opts{apms} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});
    checkExist('f', $opts{apms});

    return %opts;
}

sub apmsRound1 {
    my ($apmsFile) = @_;
    die "I thought I would find a storable but I didn't.";
}

sub rawStats {
    my ($apms, $baitStats, $preyStats) = @_;
    
    for my $bait (keys %$apms) {
	my $repl = 0+ keys %{ $apms->{$bait} };
	$baitStats->{$bait}{occur} = $repl;
	my (@tsc, @baitTsc);
	for my $expt (values %{ $apms->{$bait} }) {
	    my $tsc=0;
	    push @baitTsc, 0;
	    for my $row (@$expt) {
		my $prey = $row->{prey_ref};
		if ($prey eq $bait) {
		    $baitTsc[-1]+= $row->{total_peptides};
		} else {
		    $tsc+=$row->{total_peptides};
		    $preyStats->{$prey}{baits} //= [];
		    push @{$preyStats->{$prey}{baits}}, $bait;
		    $preyStats->{$prey}{occur}++;
		    $preyStats->{$prey}{occurSB}++ if $repl == 1;
		}
	    }
	    push @tsc, $tsc;
	}
	$baitStats->{$bait}{tsc} = \@tsc;
	$baitStats->{$bait}{baitTsc} = \@baitTsc;
    }
    die Dumper($baitStats, $preyStats);
    return;
}