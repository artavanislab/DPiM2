#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsMultiHashRef);
#use DpimLib qw(getLineDP4APMS);

# read in a flat file from reactome that maps from uniprot to reactome id

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    
    open my $IN, "<", $in or die "can't read $in. $!";
    close $IN; 
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 describe this";
    close $OUT; 
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	uni2fbgn => '/home/glocke/DPiM/uniprot/uniprot2fbgn.many2many.tsv',
	reactome => '/home/glocke/DPiM/reactome/R-DME-UniProt2Reactome_All_Levels.txt',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "uni2fbgn=s", );
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}
