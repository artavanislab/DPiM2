#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min max);
use JSON;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(parseAPMS getLineDP4APMS getLineHyperspecAPMS readHS);

# find all co-appearances of FBgn1 with FBgn2
# do not calculate HGScore because it involves NSAF etc.
my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $netFile = $opts{net};
    my $out = $opts{out};
    my $dir = $opts{dir};
    my $n = $opts{n};

    my %edges;
    readEdges(\%edges, $netFile, $n);

    my %apms;
    say "read apms";
    my %reader = (
	"4cols" => \&getLineHyperspecAPMS,
	"6cols" => \&getLineDP4APMS,
    );
    parseAPMS(\%apms, $in, $reader{$opts{mode}} // die "don't recognize $opts{mode}");
    die "no rows in apms! choose different mode?" if 0 == keys %apms;
    
    say "count co-appearances";
    my %appear; # appear{$fbgn1}{#fbgn2}{$sid}=$tsc
    my %sidCnt; # sidCnt{$sid} = number of times this sid supports an edge
    #my $inf = 1234567890;
    my $sumSum = 0;
    for my $sid (keys %apms) {
	my %row;
	my $bait = $apms{$sid}{bait};
	my $baitCnt = 0;
	for my $p (keys %{ $apms{$sid}{prey}}) {
	    next if $p eq $bait;
	    my $tsc = $apms{$sid}{prey}{$p};
	    $row{$p} = $tsc if exists $prot{$p};
	    $sumSum+= $tsc;
	    $baitCnt = max($baitCnt, $tsc);
	}
	$row{$bait} = $baitCnt if exists $prot{$bait};
	$appear{$sid} = \%row if 0 < keys%row;
    }

    say "report";
    {
	# count-up co-appearing peptides
	my %one; # one{fbgn} = sum of TSC
	my %two; # two{fbgn1-fbgn2} = sum of minTSC where both appear
	for my $sid (sort {$a <=> $b} keys %appear) {
	    my @found = sort keys %{ $appear{$sid} };
	    for my $p (@found) {
		$one{$p} += $appear{$sid}{$p};
	    }
	    if (1 < @found) {
		for my $i (0..($#found-1)) {
		    my $p1 = $found[$i];
		    for my $j (($i+1)..$#found) {
			my $p2 = $found[$j];
			my $k = "$p1\_$p2";
			$two{$k} += min($appear{$sid}{$p1}, $appear{$sid}{$p2});
		    }
		}
	    }
	}
	
	open my $OUT, ">", $out or die "can't write to $out. $!";
	say $OUT "# $0 counted tsc for $opts{prot} in $in";
	say $OUT "all\t$sumSum";
	for my $p (sort keys %one) {
	    say $OUT join "\t", $p, $one{$p};
	}
	for my $k (sort keys %two) {
	    say $OUT join "\t", $k, $two{$k};
	}
	close $OUT;
    }

    if (exists $opts{sidout}) {
	# report which experiments include appearances
	my $out2 = $opts{sidout};
	open my $OUT, ">", $out2 or die "can't write to $out2. $!";
	say $OUT "# $0 found appearances of $opts{prot} in $in";
	for my $sid (sort {$a <=> $b} keys %appear) {
	    say $OUT "$sid:", encode_json($appear{$sid});
	}
	close $OUT;
    }


}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    # mode:
    my @modes = qw(4cols 6cols);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	dir => 'workingDir/',
	n => 1000,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in apms -net FBgn1,FBgn2  -out output < ".
	"$modeString $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "prot=s", "mode=s", "dir=s", "n=i");
    die $usage unless exists $opts{in} && exists $opts{net} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{mode} //= $modes[0];
    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" if ! exists $modes{$opts{mode}};

    checkExist('f', $opts{in});

    return %opts;
}

# if $n is non-positive, read all edges
sub readEdges {
    my ($ret, $netFile, $n) = @_;

    open my $IN, "<", $in or die "Can't read from $in: $!";
    my $cnt = 0;
    while ($_ = readHS($IN)) {
	last if $n > 0 && $cnt++ >= $n;
	my ($fb1, $fb2) = split;	
	($fb1, $fb2) = sort ($fb1, $fb2);
	$ret->{$fb1}{$fb2} = 1;
    }

    return;
}
