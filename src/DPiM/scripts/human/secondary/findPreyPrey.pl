#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColRef);
use DpimLib qw(readHS);

# report how many edges connect proteins that were not baits

my %opts = getCommandLineOptions();

{
    my $netFile = $opts{net};
    my $baitFile = $opts{bait}; 
    my %outFiles = (
	0 => $opts{out0},
	1 => $opts{out1},
	2 => $opts{out2},
	);

    
    my %bait;
    {
	my @baits;
	readColRef(\@baits, $baitFile, 'bait_ref');
	%bait = map {$_ => 1} @baits;
    }


    my %OUT;
    {
	open my $NULL, ">", '/dev/null' or die "can't write to /dev/null: $!";
	my %types = ( 0 => 'prey-prey', 1 => 'bait-prey', 2 => 'bait-bait' );
	for my $n (keys %outFiles) {
	    if (defined $outFiles{$n}) { 
		open $OUT{$n}, ">", $outFiles{$n} or 
		    die "can't write to $outFiles{$n}: $!";
		
		say {$OUT{$n}} "# $0 found $types{$n} edges in $netFile (baits listed in $baitFile)";
		say {$OUT{$n}} join "\t", qw(protein1 protein2 score);
	    } else {
		$OUT{$n} = $NULL;
	    }
	}
    }
    

    open my $IN, "<", $netFile or die "Can't read from $netFile: $!";
    
    my $nEdge=0;
    my %edgesPerBaitCnt; # baitCnt{0/1/2} = 1234
    #                    # key represents the number of 
    my $doubleBait = 0;
    my $singleBait = 0;
    my $noBait = 0;
    while (my $line = readHS($IN, 'human')) {
	$nEdge++;

	my ($n1, $n2) = split /\s+/, $line;

	my $baitCnt = (exists $bait{$n1}) + exists $bait{$n2};
	$edgesPerBaitCnt{$baitCnt}++;

	print {$OUT{$baitCnt}} $line;
    }

    say join "\t", qw(nbait cnt frac);
    say join "\t", -1, $nEdge, 1;
    for my $i (0..2) {
	say join "\t", $i, $edgesPerBaitCnt{$i}, $edgesPerBaitCnt{$i}/$nEdge;
    }
    close for values %OUT;

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

    my $usage = "usage: $0 -net in.net -bait bait.list < -out0 ".
	"preyPreyEdges.tsv -out1 baitPrey.tsv -out2 baitBait.tsv >\n";

    my %opts = ();
    GetOptions(\%opts, "net=s", "bait=s", "out0=s", "out1=s", "out2=s");
    die $usage unless exists $opts{net} && exists $opts{bait};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{net});
    checkExist('f', $opts{bait});

    return %opts;
}
