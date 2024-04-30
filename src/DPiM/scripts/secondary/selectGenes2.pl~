#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist readHeader);

# select edges that have certain genes in them

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $protFile = $opts{prot};
    my $out = $opts{out};

    my %lookFor;
    findProteins(\%lookFor, $protFile);

    my $header;
    readHeader($in, \$header);

    open my $OUT, ">", $out or die "Can't write to $out. $!";

    say $OUT "# seeking proteins listed in $protFile";
    say $OUT "# selecting lines listed in $in";

    open my $IN, "<", $in or die "Can't read from $in. $!";
    while (my $line = <$IN>) {
	my $flag;
	while ($line =~ /(FBgn\d+)/g) {
	    $flag = 1 if exists $lookFor{$1};
	    last if $flag;
	}
	print $OUT $line if $flag;
    }
    close $IN;
    close $OUT;
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

    my $usage = "usage: $0 -in selectFrom network -prot contains.FBgn ".
	"-out output\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "prot=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{prot});

    return %opts;
}

sub findProteins {
    my ($lookFor, $protFile) = @_;

    open my $IN, "<", $protFile or die "Can't read from $protFile. $!";
    while (my $line = <$IN>) {
	while ($line =~ /(FBgn\d+)/g) {
	    $lookFor->{$1} = 1;
	}
    }
    close $IN;

    return;
}
