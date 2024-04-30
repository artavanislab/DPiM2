#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Find;
use Getopt::Long;

die "can't find rename command\n" if `which rename` =~ /^\/usr\/bin\/which/;


my $usage = "usage: $0 match replace <-only cpp -r -nobak>
applies s/match/replace/ to all text in specified files
  (use -r to traverse subdirectories)
  (-only only works on filenames that /match/ this arg)
  (-nobak stops the default -i.bak switch.  better be sure!!)\n";


my ($match, $replace,%opts);
my $baseArgs = 2;
if (@ARGV==$baseArgs) {
    ($match, $replace) = @ARGV;
} elsif (@ARGV>$baseArgs) {
    $match = shift(@ARGV);
    $replace = shift(@ARGV);
    die $usage unless GetOptions(\%opts,"r","nobak","only=s");
} else {
    die $usage;
}

my $dir = $ENV{PWD};
my $bak = "-i.bak";
$bak = "" if defined $opts{nobak};

my $found = 0;

sub replacer {    
    return if -d $_;
    return if defined $opts{only} && $_ !~ /$opts{only}/;
    my $in = $_;
    $_ = "";

    # protect mtime if there is no match
    open(IN, $in) or die "can't read from $_. $!\n";
    my $lclFound;
    while (<IN>) {
	if ($_ =~ /$match/) {
	    $lclFound=1;
	    $found++;
	    last;
	}
    }
    close IN;
    return unless $lclFound;

    my $cp = "$in.bak";	
    system("cp $in $cp"); # is there another way to protect permissions?
    open(IN, $cp) or die "can't read from $cp. $!\n";
    open(OUT, ">$in") or die "can't write to $in. $!\n";
    while (<IN>) {
	s/$match/$replace/g;
	print OUT;
    }

    system("rm $cp") if defined $opts{nobak};
}

if (exists $opts{r}) { 
    finddepth(\&replacer, $dir);
} else {
    find(\&replacer, $dir);    
}

if ($found) {
    print "made replacements in $found files\n";
} else {
    print "failure to match (or possibly matched but changed nothing)\n";
}
