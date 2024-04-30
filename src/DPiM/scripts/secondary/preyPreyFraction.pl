#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist readHeader);

use lib '/home/glocke/DPiM/scripts/lib';
use DpimLib qw(networkHashFromEdgeList getLineHyperspecAPMS readHS);

# what fraction of accepted edges are supported only by prey-prey interactions?
# -> given a network input, add a column saying whether it has bait-support

my %opts = getCommandLineOptions();

{
    my $netFile = $opts{net};
    my $apmsFile = $opts{apms};
    my $out = $opts{out};



    my %baitSupport; # baitSupport{$n1}{$n2} = 1 iff there is a AP-MS 
    # interaction with $n1 and $n2
    say "read network";
    networkHashFromEdgeList(\%baitSupport, $netFile);

    say "find support";
    findSupport(\%baitSupport, $apmsFile);
    my ($total, $supported) = findFraction(\%baitSupport);
    my $preComments = "# sought bait support for $netFile in $apmsFile\n";
    $preComments.= "# $supported out of $total edges have support (".
	(sprintf("%.1f", 100 * $supported/$total))."%)\n";
    
    writer($out, \%baitSupport, $netFile, $preComments);

}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	net => '/home/glocke/DPiM/data/dpim3.12232014.nrBait.network',
	apms => '/home/glocke/DPiM/data/dpim3.12232014.nrBait',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString -ternary > \n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "net=s", "apms=s", "ternary");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});
    checkExist('f', $opts{apms});

    return %opts;
}

sub findSupport {
    my ($baitSupport, $apmsFile) = @_;

    for my $n1 (keys %$baitSupport) {
	$_ = 0 for values %{$baitSupport->{$n1}};
    }

    # if ternary is defined, 
    my $map = sub { return $_[0] };
    if (exists $opts{ternary}) {
	$map = sub { return 1 };
    } 

    open my $IN, "<", $apmsFile or die "Can't read from $apmsFile. $!";
    my %row;

    while (getLineHyperspecAPMS(\%row, $IN)) {
	my ($n1, $n2) = sort($row{bait_ref}, $row{prey_ref});
	next unless exists $baitSupport->{$n1} && 
	    exists $baitSupport->{$n1}{$n2};
	$baitSupport->{$n1}{$n2} += $map->($row{total_peptides});
	if ($baitSupport->{$n1}{$n2} > 2) {
	    die Dumper("more than 2?", \%row);
	}
    }
    close $IN;

    return;
}

sub findFraction {
    my ($baitSupport) = @_;

    #                          ####
    my ($total, $supported) = (0, 0);
    #                          \o / 

    for my $n1 (keys %$baitSupport) {
	for my $n2 (keys %{$baitSupport->{$n1}}) {
	    $total++;
	    $supported++ if $baitSupport->{$n1}{$n2} > 0;
	}
    }
    return ($total, $supported);
}
sub writer {
    my ($outFile, $baitSupport, $netFile, $preComments) = @_;

    my @cols = readHeader($netFile);
    push @cols, 'baitSupport';
    
    open my $IN, "<", $netFile or die "Can't read from $netFile. $!";
    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";

    print $OUT $preComments;
    say $OUT join "\t", @cols;
    while (my $line = readHS($IN)) {
	chomp $line;
	$_ = $line;
	my @spl = split;
	my ($n1, $n2) = sort @spl[0..1];
	say $OUT $line, "\t", $baitSupport->{$n1}{$n2};
    }
    close $IN;
    close $OUT;

    return;
}
