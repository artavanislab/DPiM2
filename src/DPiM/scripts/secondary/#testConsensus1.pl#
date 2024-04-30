#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);
use Data::Dumper;
use File::Slurp;
use JSON;
use Statistics::Basic qw(variance mean median);
use HomeBrew::IO qw(checkExist readList);

#use DpimLib qw(getLineAPMS);

# collect the output results from consensus1
# do some data analysis on the results

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    my @jsonFiles = readList($in);

    say "collect data";
    # network[$i]{$fbgn}{$fbgn} = score for of i'th permutation
    my (@networks, @proteins);
    collectData(\@networks, \@proteins, \@jsonFiles);

    say "reformat"; # why didn't we do it right the first time?
    my %net2; # $net2{$fbgn}{$fbgn}[$i] == network[$i]{$fbgn}{$fbgn} // 0
    my $zeroCount = 0;
    for my $i (0..($#proteins-1)) {
	my $p1 = $proteins[$i];
	for my $p2 (@proteins[($i+1)..$#proteins]) {
	    my @arr;
	    for my $net (@networks) {
		push @arr, $net->{$p1}{$p2} // 0; 
	    }
	    if (max(@arr) == 0) {
		$zeroCount++;
		next;
	    }
	    $net2{$p1}{$p2} = \@arr;
	}
    }

    say "report";
    #writeStats0($out, \%net2, \@proteins, $zeroCount);
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

    my $usage = "usage: $0 -in json.list -out output\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

sub collectData {
    my ($collected, $proteins, $inFiles) = @_;

    my (%keys, @countKeys);
    for my $f (@$inFiles) {
	my $json = read_file($f);
	my $read = decode_json( $json );
	my $rawScore = $read->{score};
	my $ids = $read->{ids};

	$keys{$_} = 1 for values %$ids;
	push @countKeys, 0+ keys %$ids;

	my %net;
	for my $i1 (keys %$rawScore) {
	    my $fb1 = $ids->{$i1} // die "$f: can't find id{$i1}";
	    for my $i2 (keys %{ $rawScore->{$i1} }) {
		my $fb2 = $ids->{$i2} // die "$f: can't find id{$i2}";
		my ($k1, $k2) = sort ($fb1, $fb2);
		$net{$k1}{$k2} = $rawScore->{$i1}{$i2};
	    }
	}
	push @$collected, \%net;
    }
    
    push @$proteins, $_ for sort keys %keys;
}

sub writeStats0 {
    my ($out, $net, $proteins, $zeroCount) = @_;

    open my $OUT, ">", $out or die "Can't write to $out. $!";

    say $OUT "# zeroCount = $zeroCount";
    
    my @cols = qw(p1 p2 variance mean median min max);
    say $OUT join "\t", @cols;
    my $format = "%s\t%s\t".(join "\t", ('%.2e') x 5)."\n";
    for my $i1 (0..($#$proteins-1)) {
	my $p1 = $proteins->[$i1];
	for my $i2 (($i1+1)..$#$proteins) {
	    my $p2 = $proteins->[$i2];
	    my $arr = $net->{$p1}{$p2} // next;
	    printf $OUT $format, $p1, $p2, variance($arr), mean($arr)
		, median($arr), min(@$arr), max(@$arr);
	}
    }

    close $OUT;
    #$stats0{$p1}{$p2} = {
    #var => variance(\@arr),
    #mean => mean(\@arr),
    #median => median(\@arr),
    #min => min(@arr),
    #max => max(@arr),
    #}
}
