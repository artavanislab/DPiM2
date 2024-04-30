#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::MoreUtils qw(each_array);
use Storable qw(retrieve);
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineDP4APMS);

# find all co-appearances of proteins from set A with set B in apms data

my %opts = getCommandLineOptions();

{
    my %lookForA = map { $_ => 1} split /,/, $opts{p1};
    my %lookForB = map { $_ => 1} split /,/, $opts{p2};
    my $out = $opts{out};
    my $apmsFile = $opts{apms};
    
    my $apms = retrieve($apmsFile);

    my @co; # sid's with co-appearing A's and B's
    my @prot; # corresponding list of which A's and B's
    for my $bait (keys %$apms) {
	for my $sid (keys %{ $apms->{$bait} }) {
	    my (%foundA, %foundB);
	    $foundA{$bait} = 1 if exists $lookForA{$bait};
	    $foundB{$bait} = 1 if exists $lookForB{$bait};
	    for my $row (@{ $apms->{$bait}{$sid} }) {
		my $prey = $row->{prey_ref};
		next if $prey eq $bait;
		$foundA{$prey} = 1 if exists $lookForA{$prey};
		$foundB{$prey} = 1 if exists $lookForB{$prey};
	    }
	    if ((0 < keys %foundA) && (0 < keys %foundB)) {
		# double check that we haven't found just one protein that
		#   appears in both lists.  Given that there is at least one
		#   protein in both lists, if there are at least two different
		#   proteins, then there is at least one bona fide pair.
		my %all = map {$_ => 1} (keys %foundA, keys %foundB);
		if (1 < keys %all) {
		    push @co, $sid;
		    push @prot, join ",", (keys %foundA, keys %foundB);
		}
	    }
	}
    }

    if (@co > 0) {
	my $it = each_array(@co, @prot);
	say "sid\tprots";
	while (my ($sid, $prots) = $it->()) {
	    say "$sid\t$prots";
	}
    } else {
	say "found none!";
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	apms => '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.01-30-2016.storable',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -p1 fbgnA1,fbgnA2,... -p2 fbgnB1,fbgnB2,... < ".
	"$defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "p1=s", "p2=s", "apms=s");
    die $usage unless exists $opts{p1} && exists $opts{p2};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{apms});

    return %opts;
}
