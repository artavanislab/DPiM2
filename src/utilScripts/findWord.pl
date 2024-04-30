#!/usr/bin/env perl

use strict;
use Data::Dumper;


# list all locations of a specific DNA word in a given fasta sequence

my %opts = &getCommandLineOptions();
my $w = $opts{word};

print "# finding $w and ".rc($w)."\n" if defined $opts{rc};

my $chr = [];
my @seqs = read_fasta_glocke($opts{fa},$chr);



my ($i,$j,@found);
for $i (0..$#$chr-1) { # check for duplicate keys;
    for $j ($i+1..$#$chr) {
	$chr->[$j].="1" if $chr->[$i] eq $chr->[$j];
    }
}

for $i (0..$#seqs) {
    $found[$i] = [];
    $j=-1;
    do {
	$j = index($seqs[$i],$w,$j+1);
	push(@{$found[$i]},$j+1) unless $j<0;
    } while $j>=0;
}

print "#found the following (one-based) locations of $w: \n";
for (0..$#found) {
    print ">".$chr->[$_].": ".scalar(@{$found[$_]})." locations.\n".
	join(", ",@{$found[$_]})."\n";
    printf "this means %.5f%% of ".$chr->[$_].
	" is covered by $w (assuming no overlapping words)\n\n",
	( 100*scalar(@{$found[$_]})*length($w)/length($seqs[$_]) );
}

exit unless defined $opts{rc};

my $old = $w;
$w = rc($w);
if ($w eq $old) {
    die "rc($w) eq $w\n";
}

@found = ();

for $i (0..$#seqs) {
    $found[$i] = [];
    $j=-1;
    do {
	$j = index($seqs[$i],$w,$j+1);
	push(@{$found[$i]},$j+1) unless $j<0;
    } while $j>=0;
}

print "#found the following (one-based) locations of $w: \n";
for (0..$#found) {
    print ">".$chr->[$_].": ".scalar(@{$found[$_]})." locations.\n".join(", ",@{$found[$_]})."\n";
    printf "this means %.5f%% of ".$chr->[$_]." is covered by $w (assuming no overlapping words)\n\n",(100*scalar(@{$found[$_]})*length($w)/length($seqs[$_]));
}


exit;

####################
###  Functions  ####
####################

sub getCommandLineOptions {
  use Getopt::Long;
  my $usage = qq{usage: $0 -fa fasta.fna -word ACGT <-rc>\n};
  # Get args
  my %opts = ();
  &GetOptions(\%opts,"fa=s","word=s","rc");
  if (!defined $opts{fa} || !defined $opts{word}) {
    print STDERR "$usage\n";
    exit -1;
  }

  checkExist('f',$opts{fa});

  return %opts;
}

sub checkExist {
  my ($type, $path) = @_;
  if ($type eq 'd') {
    if (! -d $path) {
      warn("dir not found: $path\n");
      exit -3;
    }
  }
  elsif ($type eq 'f') {
    if (! -f $path) {
      warn("file not found: $path\n");
      exit -3;
    }
    elsif (! -s $path) {
      warn("empty file: $path\n");
      exit -3;
    }
  }
}


#This sub loads fasta seqs into an array WITHOUT BioPerl
#functions. It can optionally perform various checks on
#the seqs it read.
#WARNING: this sub is NOT thoroughly tested!
sub read_fasta_glocke {
    my ($fasta_file,$chr,$mode) = @_;
    my @seqs = ();
    my $cnt = -1;
    open(II,"$fasta_file") or die "Can't open file: $!\n";
    while (<II>) {
	chomp;
	my @tt = split ' ';
	next if (@tt==0);
	next if substr($tt[0],0,1) eq '#';
	if (substr($tt[0],0,1) eq '>') {
	    if (defined $chr) {
		push(@$chr,substr($tt[0],1))
	    }
		
	    $cnt++;
	    next;
	}
	die "Unrecognized FASTA seq!\n" unless (@tt==1);
	if (!defined $seqs[$cnt]) {
	    $seqs[$cnt] = uc($tt[0]); # make everything uppercase!!!
	}else{
	    $seqs[$cnt] .= uc($tt[0]); # make everything uppercase!!!
	}
    }

    close II;
    if (defined $mode) {
	if ($mode eq "check_width") { # checks if all seqs are of the same length
	    for my $i (1 .. $#seqs) {
		die "Seqs read from $fasta_file are not of the same length!\n"
		    if ((length($seqs[$i]) - length($seqs[$i-1]) ) != 0);
	    }
	}
    }
    return @seqs;
}

sub rc {
    my $w = shift;
    my $revcomp = reverse($w);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}
