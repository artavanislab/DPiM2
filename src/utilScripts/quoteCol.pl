#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsRef readHeader);
#use DpimLib qw(getLineDP4APMS);

# put quotes around a particular column

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $col = $opts{col};
    my $out = $opts{out};
    my $sep = $opts{sep};
    if ($sep eq 'tab') {
	$sep = "\t";
    }
    
    my @allCols = readHeader($in);
    {
	my %check = map {$_ => 1} @allCols;
	die "can't find '$col' column in $in" unless exists $check{$col};
    }

    my @read;
    readColsRef(\@read, $in, \@allCols, undef, $sep);

    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# added quotes around $col column in $in";
    say $OUT join $sep, @allCols;
    for my $row (@read) {
	$row->{$col} = '"'.$row->{$col}.'"';
	say $OUT join $sep, map {$row->{$_}} @allCols;
    }
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	sep => 'tab',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -col name -out output < $defaultString ".
	">\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "col=s", "out=s", "sep=s");
    die $usage unless exists $opts{in} && exists $opts{col} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
