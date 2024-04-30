#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineAPMS);

# set difference of two comma separated lists

die "usage: $0 apple1,apple2,.. orange1,apple1,...\n" unless @ARGV==2;

{
    my %l1 = map { $_ => 1 } split /,/, $ARGV[0];
    my %l2 = map { $_ => 1 } split /,/, $ARGV[1];
    
    say "missing from first";
    say join ",", findMissing(\%l2, \%l1);
    
    say "missing from second";
    say join ",", findMissing(\%l1, \%l2);
    
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub findMissing {
    my ($h1, $h2) = @_;

    my @ret;
    for my $k (keys %$h1) {
	push @ret, $k if ! exists $h2->{$k};
    }
    return @ret;
}
