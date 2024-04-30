#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist);

# describe the purpose for this script

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $transFile = $opts{trans};

    my %fbgnMap;
    say "parsing $transFile...";
    makeMap(\%fbgnMap, $transFile);

    say "translating...";
    open my $IN, "<", $in or die "Can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";
	
    if (exists $opts{addcols}) {
	addCols($IN, $OUT, \%fbgnMap, $opts{addcols});
    } else {
	while (my $line = <$IN>) {
	    while ($line =~ /(FBgn\d+)/) {
		my $two = $fbgnMap{$1} // 'UNKNOWN';
		$line =~ s/$1/$two/;
	    }
	    print $OUT $line;
	}
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	trans => '/home/glocke/DPiM/flybase/fbgn_annotation_ID_fb_2015_04.tsv'
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString -addcols <number of columns> >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "trans=s", "addcols=i");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{trans});

    return %opts;
}

sub makeMap {
    my ($ret, $mapFile) = @_;

    open my $IN, "<", $mapFile or die "Can't read from $mapFile. $!";
    while (<$IN>) {
	next if /^#/;
	next if length($_) < 2;
	chomp;
	my @spl = split;
	$ret->{$spl[1]} = $spl[0];
	if (@spl > 2) {
	    my @spl2 = split ',', $spl[2];
	    $ret->{$_} = $spl[0] for @spl2;
	}
	#die Dumper($ret);
    }
    return;
}

sub addCols {
    my ($IN, $OUT, $map, $nCols) = @_;

    my $extraCols = "geneSymbol";
    if ($nCols > 1) {
	$extraCols = join "\t", map { "geneSymbol$_" } 1..$nCols;
    }

    my $line = <$IN>;
    while ($line =~ /^#/) {
	print $OUT $line;
	$line = <$IN>;
    } 
    chomp $line;
    say $OUT "$line\t$extraCols";

    while ($line = <$IN>) {
	chomp $line;
	my @genes;
	while ($line =~ /(FBgn\d+)/g) {
	    push @genes, $map->{$1} // "NOT_FOUND";
	}
	say $OUT $line, "\t", join "\t", @genes;
    }
    return;
}
