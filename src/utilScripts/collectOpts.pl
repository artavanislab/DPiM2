#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;

my $optString = join " ", @ARGV;

my @args;
while ($optString =~ /"(\S+)=/g) {
    push @args, $1;
}

die "can't parse '$optString'" if @args == 0;

for my $a (@args) {
    say "    my \$$a = \$opts{$a};";
}

