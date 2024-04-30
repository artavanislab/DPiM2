#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(max sum);
use Math::Interpolate qw(linear_interpolate);
use HomeBrew::IO qw(checkExist readList readCols);
use Statistics::Basic qw(median);
use DpimLib qw(readHS);

# find FDR at a given HGScore from simulations

my %opts = getCommandLineOptions();

{
    my $realFile = $opts{real};
    my $simFile = $opts{sim};
    my $out = $opts{out};

    my @simFiles;
    if (exists $opts{simlist}) {
	@simFiles = readList($simFile);
	checkExist('f', $_) for @simFiles;
    } else {
	@simFiles = ($simFile);
    }

    my @realNet; # realNet[i] = {p1, p2, score}
    readNet(\@realNet, $realFile);
    
    my @interps = ();
    {
	my @realScores = map { $_->{score} } @realNet;
	# interp[$i](score) = linear interpolation for FDR at a given score
	for my $sf (@simFiles) {
	    say "interpolate $sf";
	    push @interps, mapFDR(\@realScores, $sf);
	    say "interps[-1](745) = ".$interps[-1](745);
	}
    }

    say "print";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 calculated FDR for each score from $realFile and $simFile";
    say $OUT join "\t", qw(protein1 protein2 score fdr);
    for my $row (@realNet) {
	my $fdr = sum( map { $_->($row->{score}) } @interps);
	$fdr = max(0, $fdr/@simFiles);
	say $OUT join "\t", $row->{p1}, $row->{p2}, $row->{score}, $fdr;
    }
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
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -real hypspec.o123 -sim hypspec.o123 < -out output ".
	"-simlist $modeString>\n";

    my %opts = ();
    GetOptions(\%opts, "real=s", "sim=s", "out=s", "simlist", "mode=s");
    die $usage unless exists $opts{real} && exists $opts{sim} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{mode} //= $modes[0];
    die "you must select one of the following modes: "
	, join(", ",keys(%modes)), "\n" if ! exists $modes{$opts{mode}};

    checkExist('f', $opts{real});
    checkExist('f', $opts{sim});

    return %opts;
}

# ret[i] = {p1, p2, score}
sub readNet {
    my ($ret, $realFile) = @_;

    open my $IN, "<", $realFile or die "can't read from $realFile: $!";
    while(readHS($IN)) {
	chomp;
	my @spl = split;
	push @$ret, { p1 => $spl[0], p2 => $spl[1], score => $spl[2] };
    }
}


# ret{fdr} = {score, edges}
# edges = size of network at this score
sub mapFDR {
    my ($realScores, $simFile, $fdrMax) = @_;

    $fdrMax //= 0.2; # stop mapping after FDR = 

    my @simulatedScores = ();
    open my $IN, "<", $simFile or die "can't read from $simFile. $!";
    while (my $line = <$IN>) {
	if ($line =~ "^FBgn" || $line =~ /^\d+\s+\d+\s+\d+/){
	    chomp($line);
	    my @field = split(/\t/, $line);
	    push @simulatedScores, $field[2];
	}
    }
    close($IN);
    die "can't find any simulated edges.  maybe you want to use -simlist?\n" 
	if @simulatedScores == 0;
    
    my $s=0;
    my $FDR=0;
    my (@score, @fdr);
    while ($FDR < $fdrMax && $s < $#simulatedScores) {
	my $r=0;
	while ($realScores->[$r] >= $simulatedScores[$s]) {
	    $r++;
	    die "realScores->[$r] not defined, fdr = $FDR.  \$#realScores = $#$realScores\n"
		if ! defined $realScores->[$r];
	}
	$FDR = ($s+1)/($r+1);
	$s++;
	#warn "$FDR = $fdr[-1] || $realScores->[$r+1] == $score[-1])" if 0 < @fdr;
	next if 0 < @fdr && ($FDR = $fdr[-1] || $realScores->[$r+1] == $score[-1]);
	push @fdr, $FDR;
	push @score, $realScores->[$r+1];
    }

    return sub {
	my $scoreX = shift;
	return( (linear_interpolate($scoreX, \@score, \@fdr))[0] );
    };
}

