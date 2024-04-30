#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsHashRef);
use DpimLib qw(getLineDP4APMS);

# read in apms data and search for fbgn's absent from the nsaf length file

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $nsafFile = $opts{nsaf};

    my %nsaf;
    readColsHashRef(\%nsaf, $nsafFile, [qw(fbgn avg_length)]);
    
    open my $IN, "<", $in or die "can't read $in. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 searching for fbgn's in $in not found in $nsafFile";
    my %row;
    my %notFound; # notFound{FBgn} = 1 if I've already identified this protein
    while (getLineDP4APMS(\%row, $IN)) {
	print $OUT check($row{bait_ref}, \%nsaf, \%notFound);
	print $OUT check($row{prey_ref}, \%nsaf, \%notFound);
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
	nsaf => '/home/glocke/DPiM/nsaf/dmel-all-translation-r6.09.aveLen.tsv',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in apms.dp4 -out output < $defaultString > \n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "nsaf=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{nsaf});

    return %opts;
}

# if there is no length for $fb (and I don't already know that) then return $fb
# otherwise, return '';
sub check {
    my ($fb, $nsaf, $notFound) = @_;

    my $ja = $fb;
    my $nein = '';
    
    return $nein if exists $notFound->{$fb};
    if (! exists $nsaf->{$fb}) {
	$notFound->{$fb} = 1;
	return $ja
    } 
    return $nein;
}
