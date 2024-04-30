#!/usr/bin/env perl

use strict;
use Data::Dumper;
use File::Find;
use feature qw(say);

die "can't find rename command\n" unless `which rename` !~ /^\/usr\/bin\/which/;

my $dir = $ENV{PWD};

unless (@ARGV==2 || (@ARGV==3 && $ARGV[2] eq "-r")) {
    my $usage = "usage: $0 match replace <-r>
applies s/match/replace/ to move files/directories 
  (use -r to traverse subdirectories)
WARNING: nameclashes result in lost files.  be careful!\n";
    print $usage;
    exit;
}

my $match = $ARGV[0];
my $replace = $ARGV[1];

#say $match;
#say $replace;


my $found;

sub renamer {    
    my $out = $_;
    $out =~ s/$match/$replace/;
    if ($out ne $_) {
	my $cmd = "mv $_ $out";
	#say $cmd;
	system($cmd);
	$found = 1;
    }
}

sub test { print "$_\n"; }

if (@ARGV==3) { 
    finddepth(\&renamer, $dir);
} else {
    find(\&renamer, $dir);    
}

print "failure to match (or possibly matched but changed nothing)\n"
    unless $found;
