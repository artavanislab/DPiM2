#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;


opendir my $DIR, $ENV{PWD};
my @allFiles = readdir($DIR);
closedir $DIR;

my $ext = 'sh';
my @scr = grep { /\.$ext$/ } @allFiles;
die "found no *.$ext" unless @scr > 0;

$_ =~ s/\.sh// for @scr; 

my @outFiles = grep { /.o\d+/ } @allFiles;
my %outFiles; 
for my $o (@outFiles) {
    $o =~ s/\.o\d+//;
    $outFiles{$o} = 1;
}

my @notFound = grep { ! exists $outFiles{$_} } @scr;

say join "\n", @notFound;

warn (0+ @notFound)." out of ".(0+ @scr)." did not complete\n" if @notFound;
