#/usr/bin/env perl

use warnings;
use strict;

use POSIX;
use Test::More qw(no_plan);
use Data::Dumper;

BEGIN { use_ok('Test::Sort') };

my $numSort = sub ($$) { $_[0] <=> $_[1] };
my @nums = qw( 3 2 4 );
my @sorted1 = Test::Sort::testSort($numSort, @nums);
my @sorted2 = sort $numSort @nums;
print Dumper(\@sorted1);
print Dumper(\@sorted2);
is_deeply(\@sorted1, \@sorted2, "can I sort?");
