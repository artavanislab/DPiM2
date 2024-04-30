#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList readCols readColsHash);
#use DpimLib qw(getLineDP4APMS);

# 361 entrez id's are not found in the ncbi protein reference i downloaded
# I plugged them into DAVID and 

my %opts = getCommandLineOptions();

{
    my $allMissFile = $opts{all};
    my $lenFile = $opts{len};
    my $davidFile = $opts{david};
    my $baseOut = $opts{out};

    my @all = readList($allMissFile);
    my @david = readCols($davidFile, [qw(To From)], undef, "\t");
    my %len = readColsHash($lenFile, [qw(n.symbol n.wt_length)], ",");

    my %davHas = map {$_->{From} => 1} @david;
    

    my %sumLength;
    my %sumN;
    my %stillMiss;
    for my $row (@david) {
	my $trz = $row->{From};
	my $sym = uc($row->{To});
	if (! exists $len{$sym}) {
	    if (defined $stillMiss{$trz}) {
		$stillMiss{$trz}.="_$sym";
	    } else {
		$stillMiss{$trz} = $sym;
	    }
	} else {
	    $sumLength{$trz} += $len{$sym} // die $trz;
	    $sumN{$trz}++;
	}
    }

    {
	my $out = "$baseOut.length";
	open my $OUT, ">", $out or die "can't write to $out. $!";
	for my $trz (sort {$a <=> $b} keys %sumN) {
	    say $OUT join "\t", $trz, -1, $sumLength{$trz} / $sumN{$trz};
	}
	close $OUT;
    }
    {    
	my $out = "$baseOut.stillMissing";
	open my $OUT, ">", $out or die "can't write to $out. $!";
	for my $trz (sort {$a <=> $b} @all) {
	    say $OUT "$trz\tnotInDavid" if ! exists $davHas{$trz};
	}
	for my $trz (sort {$a <=> $b} keys %stillMiss) {
	    say $OUT "$trz\tnotInCsv_$stillMiss{$trz}";
	}
	close $OUT; 
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	all => '/home/glocke/DPiM/human/nsaf/noSymb.log2',
	david => '/home/glocke/DPiM/human/nsaf/noSymb.david.tsv',
	len => '/home/glocke/DPiM/human/nsaf/wild_type_human_gene.csv',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out base.out < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "all=s", "david=s", "len=s");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{all});
    checkExist('f', $opts{david});
    checkExist('f', $opts{len});

    return %opts;
}
