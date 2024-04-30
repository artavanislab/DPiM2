#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineBiopAPMS);

# 1. reject any prey unless it has at least two peptides in both replicates
# 2. reject any prey unless it passes the entropy test (description below)

my $entropyDescription = qq[
ENTROPY LC CARRY OVER CORRECTION

> Input: spectral counts for a given protein in both technical 
> replicates
s1 = # of spectral counts for protein A in replicate 1
s2 = # of spectral counts for protein A in replicate 2

> Calculate proportion of spectral counts in each run Note that an extra 
> pseudocount is split between the two replicates because zeroes will 
> cause issues in the entropy formula
p1 = (s1 + 0.5) / (s1 + s2 + 1.0)
p2 = (s2 + 0.5) / (s1 + s2 + 1.0)

> Entropy calculation (based on Shannon's entropy)
E = p1 * log2(p1) + p2 * log2(p2)

> Filter:
If E < 0.75, REJECT AS ARTIFACT (e.g. LC Carry-over)
];


my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $minEntropy = $opts{minentropy};
    $minEntropy *= log(2);

    my %apms;
    parseAPMS(\%apms, $in);
    
    my @cols = qw( search_id bait_ref prey_ref total_peptides sample_date ms_inst_run_id  rep plate );    

    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 removed bait-prey line unless it has at least two peptides in both replicates";
    say $OUT join "\t", @cols;
    for my $bait (sort {$a <=> $b } keys %apms) {
	for my $p (sort { $apms{$bait}{A}{$b}{tsc} <=> 
			      $apms{$bait}{A}{$a}{tsc} } 
		   keys %{ $apms{$bait}{A} })
	{
	    my $pepA = $apms{$bait}{A}{$p}{tsc};
	    my $pepB = $apms{$bait}{B}{$p}{tsc} // 0;
	    last if 2 > $pepA;
	    next if 2 > $pepB;
	    next unless checkEntropy($pepA, $pepB, $minEntropy);
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
	minentropy => 0.75,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "minentropy=f");
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
	$ret->{$row{bait_ref}}{$row{rep}}{$row{prey_ref}} = {
	    tsc => $row{total_peptides},
	    line => $row{line}
	};
	#die Dumper($ret);
    }

    return;
}

# note that we're not using log_2 because we assume that $min has already been
#   scaled 
sub checkEntropy {
    my ($tscA, $tscB, $min) = @_;

    my $pA = ($tscA + 0.5) / ($tscA + $tscB + 1);
    my $pB = ($tscB + 0.5) / ($tscA + $tscB + 1);

    my $E = -($pA * log($pA) + $pB * log($pB));
	#print "$E\n";
    return undef if $E < $min;
    return 1;
}
