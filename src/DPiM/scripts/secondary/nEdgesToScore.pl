#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist writeCols);

# look at the HyperSpec main output file, find the number edges necessary
# to reach FDR.  Then go to the simulation output and find the score cutoff
# implied by this number of edges

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $real = $opts{real};
    my $out = $opts{out};

    my @maxEdge = getMaxEdge($in);

    my %lines = map { $_ => 1 } @maxEdge;

    my $lineN = 1;
    my %score; # score{lineN} = HGScore for that edge
    open my $IN, "<", $real or die "can't read from $real. $!";
    while (my $line = readHS($IN)) {
	$score{$lineN} = getScore($line) if exists $lines{$lineN};
	$lineN++;
    }
    close $IN;
    
    my @score = map { $score{$_} } @maxEdge;

    my $preComments = "# found score linked to line numbers in $in, $real\n";
    my $header = "line\tscore\n";
    my $format = "%d\t%f\n";
    my @data = (\@maxEdge, \@score);
    writeCols($out, \@data, $header, $format, $preComments);
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

    my $usage = "usage: $0 -in hyperspec.out -real hyper.r -out output\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "real=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{real});

    return %opts;
}

sub getMaxEdge {
    my ($hypOut) = @_;

    open my $IN, "<", $hypOut or die "can't read from $hypOut. $!";

    my @ret;
    while(my $line = <$IN>) {
	next unless $line =~ /edges/;
	$_ = $line;
	my @spl = split;
	next unless $spl[1] =~ /^\d+$/;
	push @ret, $spl[4];
    }
    
    close $IN;

    return @ret;
}

sub readHS {
    my ($IN) = @_;

    return undef if eof($IN);

    my $goodLine = sub {
	return $_[0] =~ /^FBgn\d+\s+FBgn/;
    };
    
    my $line;
    do {
	$line = <$IN>;
    } while (!eof($IN) && ! $goodLine->($line));
    return undef if ! $goodLine->($line);

    $_ = $line;
    chomp;
    return $line;
}

sub getScore {
    my ($line) = @_;
    $_ = $line;
    chomp;
    my @spl = split;
    return $spl[-1];
}
